#pragma once

#include <algorithm>
#include <functional>
#include <memory>
#include <numeric>
#include <optional>
#include <random>
#include <span>
#include <sstream>
#include <cmath>
#include <string_view>
#include "reconstruction_state.hpp"
#include "recurrence_sequence.hpp"
#include <stdexcept>
#include <vector>
#include "lil_matrix.hpp"
#include "../timer_accum.hpp"
#include "block_minimal_generator.hpp"
#include "sparse_transpose.hpp"
#include "packed_projection.hpp"

namespace VectorSpace {

namespace block_wiedemann_detail {
// Each output entry is an independent dot product; no reductions or atomics.
template <Field K>
void project_block(std::span<const K> projection,std::span<const K> current,
    std::span<K> moment,std::size_t n,std::size_t b,bool transposed,int threads) {
    const bool parallel=threads>1 && n*b*b>=32768;
    (void)parallel;
#pragma omp parallel for collapse(2) schedule(static) num_threads(threads) if(parallel)
    for(std::size_t i=0;i<b;++i)
        for(std::size_t j=0;j<b;++j) {
            K sum{};
            for(std::size_t r=0;r<n;++r)sum+=projection[r*b+i]*current[r*b+j];
            moment[transposed?j*b+i:i*b+j]=sum;
        }
}
} // namespace block_wiedemann_detail

// Block vectors are row-major n-by-b arrays. The matrix callbacks overwrite
// their outputs. Sparse constructors borrow immutable matrix data and own a
// transpose; borrowed data must outlive this solver.
template <Field K>
class block_wiedemann_solver {
public:
    struct timing_t {
        timer_accum gram, sequence, generator, reconstruction;
    };
    mutable timing_t timing;
    struct sequence_statistics {
        std::size_t sequences=0, moments=0, recurrence_updates=0, validation_attempts=0;
    };
    mutable sequence_statistics sequence_stats;
    void print_timing() const {
        timing.gram.print(direct_square_ ? "Square operator applications (included in sequence/reconstruction)"
            : "Gram applications (included in sequence/reconstruction)");
        timing.sequence.print("Krylov sequence and projections");
        timing.generator.print("block minimal generator");
        timing.reconstruction.print("nullspace reconstruction and verification");
    }
    using Block = std::vector<K>;
    using Apply = std::function<void(std::span<const K>,std::span<K>,std::size_t)>;
    using ImageVec = std::unique_ptr<K[]>;
    using DomainVec = std::unique_ptr<K[]>;
    struct options {
        std::size_t block_size=8, trials=1, holdout=8;
        std::uint64_t seed=17;
        int threads=8;
        std::ostream* log=nullptr; // Borrowed stream; null keeps library use quiet.
        bool incremental_recurrence=false; // Opt-in comparison; batch remains the default.
        std::size_t recurrence_check_interval=32;
        std::size_t sequence_capacity=0; // 0 reserves the initial batch plus one extension.
        std::filesystem::path checkpoint_path;
        bool resume_checkpoint=false;
        std::filesystem::path recurrence_sequence_path; // Optional standalone comparison input.
        double reconstruction_checkpoint_seconds=60; // 0 saves every complete step.
        double checkpoint_seconds=60; // Krylov, recurrence and validation boundaries.
        std::function<void(std::string_view)> on_checkpoint; // After durable publication.
    };
    struct rank_result {
        std::size_t rank=0, nullity=0;
        std::vector<std::size_t> trial_ranks;
        bool probabilistic=true;
    };
    struct nullspace_result {
        // Each vector uses the original domain coordinates. Vectors are
        // independently verified; completeness depends on the rank estimate.
        std::vector<Block> basis;
        rank_result rank_estimate;
        bool complete=false;
    };
private:
    std::size_t rows_,cols_;
    Apply apply_,transpose_;
    options options_;
    bool direct_square_ = false;

    // Only the controlling thread reports. Flush so redirected logs stay live.
    struct progress {
        std::ostream* out;
        const char* phase;
        timer_accum::clock::time_point start=timer_accum::clock::now(),last=start;
        std::chrono::system_clock::time_point calendar_start=std::chrono::system_clock::now();
        void operator()(std::size_t done,std::size_t total) {
            if(!out)return;
            const auto now=timer_accum::clock::now();
            if(done && done!=total && now-last<std::chrono::seconds(5))return;
            *out << phase << " " << done << "/" << total << " elapsed_seconds="
                 << std::chrono::duration<double>(now-start).count()
                 << " calendar_elapsed_seconds=" << std::chrono::duration<double>(
                     std::chrono::system_clock::now()-calendar_start).count() << std::endl;
            last=now;
        }
    };

    static K nonzero(std::mt19937_64& rng) {
        K x;
        do { x=K::sample(rng); } while(x==K{});
        return x;
    }
    struct preconditioner {
        std::size_t n,other;
        Apply forward,backward;
        Block left,right;
        Block scratch;
        timer_accum* elapsed;
        bool direct = false;
        int threads = 1;
        void apply(std::span<const K> in,std::span<K> out,std::size_t b) {
            timer_accum::guard measured(*elapsed);
            if (direct) {
                forward(in,out,b);
                const bool parallel=n*b>=32768;
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
                for(std::size_t i=0;i<n;++i)
                    for(std::size_t j=0;j<b;++j) out[i*b+j]*=left[i];
                return;
            }
            scratch.resize(other*b);
            backward(in,scratch,b);
            for(std::size_t i=0;i<other;++i)
                for(std::size_t j=0;j<b;++j) scratch[i*b+j]*=right[i];
            forward(scratch,out,b);
            for(std::size_t i=0;i<n;++i)
                for(std::size_t j=0;j<b;++j) out[i*b+j]*=left[i];
        }
    };
    preconditioner precondition(std::mt19937_64& rng,bool swapped) const {
        if(direct_square_) {
            // Left scaling preserves the right kernel. Apply it AFTER the supplied
            // composition has reduced its integer accumulators modulo p.
            preconditioner p{cols_,0,apply_,{},{},{},{},&timing.gram,true,options_.threads};
            p.left.resize(p.n);
            for(auto& x:p.left)x=nonzero(rng);
            return p;
        }
        preconditioner p{swapped?cols_:rows_,swapped?rows_:cols_,
            swapped?transpose_:apply_,swapped?apply_:transpose_,{},{},{},&timing.gram};
        p.left.resize(p.n);p.right.resize(p.other);
        for(auto& x:p.left)x=nonzero(rng);
        for(auto& x:p.right)x=nonzero(rng);
        return p;
    }
    using Generator=block_wiedemann_detail::packed_generator<K>;
    struct kernel_start { preconditioner p; Block v; Generator gen; Block initial,projection; };
    std::optional<Generator> sequence_generator(preconditioner& p,const Block& initial,
        const Block& projection,std::size_t b,bool transposed_moments,bool allow_resume=true,const std::string& checkpoint_suffix="") const {
        const std::filesystem::path checkpoint_path=options_.checkpoint_path.empty() ? std::filesystem::path{}
            : std::filesystem::path(options_.checkpoint_path.string()+checkpoint_suffix);
        ++sequence_stats.sequences;
        Block current=initial,next(p.n*b);
        Block packed_projection,packed_current;
        if constexpr(std::same_as<K,Z32783>) {
            packed_projection.resize(p.n*b);packed_current.resize(p.n*b);
            block_wiedemann_detail::pack_projection<K>(projection,packed_projection,p.n,b,options_.threads);
            if(options_.log)*options_.log << "projection=packed_uint64 packed_bytes="
                << 2*p.n*b*sizeof(K) << std::endl;
        }
        const auto initial_training=2*((p.n+b-1)/b)+2;
        const auto limit=2*p.n+2;
        const auto interval=options_.recurrence_check_interval;
        const auto capacity=options_.sequence_capacity ? options_.sequence_capacity : std::min(limit,2*initial_training)+options_.holdout;
        if(capacity<=options_.holdout)throw std::invalid_argument("sequence capacity must exceed holdout");
        block_wiedemann_detail::packed_moments<K> moments(b,capacity);
        block_wiedemann_detail::minimal_generator_state<K> state(moments,b,options_.threads);
        Generator gen(b,capacity/2);
        bool accepted=false;
        auto training=options_.incremental_recurrence ? std::min(interval,initial_training) : initial_training;
        std::size_t processed=0;
        auto last_checkpoint=std::chrono::steady_clock::now();
        const auto checkpoint=[&](std::string_view phase) {
            if(checkpoint_path.empty())return;
            block_wiedemann_detail::checkpoint_output out(checkpoint_path);
            out.section("solver-state",[&]{
            out.word(p.n);out.word(p.other);out.word(b);out.word(options_.seed);
            out.word(options_.holdout);out.word(options_.incremental_recurrence);out.word(transposed_moments);
            out.word(p.direct);out.word(training);out.word(sequence_stats.sequences);
            out.word(sequence_stats.moments);out.word(sequence_stats.recurrence_updates);out.word(sequence_stats.validation_attempts);});
            out.template block<K>("initial",initial,p.n,b);
            out.template block<K>("projection",projection,p.n,b);
            out.template block<K>("left-diagonal",p.left,p.left.size(),1);
            out.template block<K>("right-diagonal",p.right,p.right.size(),1);
            out.template block<K>("current-krylov",current,p.n,b);
            moments.save(out);state.save(out);
            out.section("accepted-recurrence",[&]{out.word(accepted);if(accepted)gen.save(out);});out.finish();
            last_checkpoint=std::chrono::steady_clock::now();
            if(options_.log)*options_.log << "checkpoint saved phase=" << phase << " moments=" << moments.size()
                << " recurrence_terms=" << state.processed_terms() << " path=" << checkpoint_path << std::endl;
            if(options_.on_checkpoint)options_.on_checkpoint(phase);
        };
        const auto periodic=[&](std::string_view phase) {
            if(!checkpoint_path.empty() && std::chrono::duration<double>(
                std::chrono::steady_clock::now()-last_checkpoint).count()>=options_.checkpoint_seconds)checkpoint(phase);
        };
        if(allow_resume && options_.resume_checkpoint && !checkpoint_path.empty() && std::filesystem::exists(checkpoint_path)) {
            if(checkpoint_path.empty())throw std::invalid_argument("resume requires checkpoint_path");
            block_wiedemann_detail::checkpoint_input in(checkpoint_path);
            in.section("solver-state",[&]{
            in.expect(p.n);in.expect(p.other);in.expect(b);in.expect(options_.seed);
            in.expect(options_.holdout);in.expect(options_.incremental_recurrence);in.expect(transposed_moments);
            in.expect(p.direct);training=in.word();sequence_stats.sequences=in.word();
            sequence_stats.moments=in.word();sequence_stats.recurrence_updates=in.word();sequence_stats.validation_attempts=in.word();});
            if(training>limit)throw std::runtime_error("checkpoint training target exceeds solver limit");
            next.resize(std::max(initial.size(),p.right.size()));
            const auto compare=[&](const char* role,std::span<const K> expected,std::size_t rows,std::size_t cols) {
                auto buffer=std::span<K>(next.data(),expected.size());
                in.template block<K>(role,buffer,rows,cols);
                if(!std::equal(buffer.begin(),buffer.end(),expected.begin()))
                    throw std::runtime_error("checkpoint does not match operator, seed or projections");
            };
            compare("initial",initial,p.n,b);compare("projection",projection,p.n,b);
            compare("left-diagonal",p.left,p.left.size(),1);compare("right-diagonal",p.right,p.right.size(),1);
            in.template block<K>("current-krylov",current,p.n,b);
            moments.load(in);state.load(in);
            in.section("accepted-recurrence",[&]{
                const auto saved=in.word();if(saved>1)throw std::runtime_error("invalid recurrence acceptance flag");
                accepted=saved;if(accepted)gen.load(in);
            });in.finish();
            processed=state.processed_terms();next.resize(p.n*b);
            if(options_.log)*options_.log << "checkpoint resumed moments=" << moments.size()
                << " processed=" << processed << " capacity=" << capacity << std::endl;
        }
        if(accepted)return std::optional<Generator>(std::move(gen));
        if(options_.log)*options_.log << "packed_recurrence_bytes="
            << moments.allocated_bytes()+state.allocated_bytes()+gen.allocated_bytes()
            << " sequence_capacity=" << capacity << std::endl;
        progress report{options_.log,options_.incremental_recurrence ? "Krylov + incremental recurrence" : "Krylov sequence"};
        report(0,initial_training+options_.holdout);
        if(options_.log)*options_.log << "recurrence_mode="
            << (options_.incremental_recurrence ? "incremental" : "batch") << std::endl;
        checkpoint("sequence-start");
        for(;;) {
            while(moments.size()<training+options_.holdout) {
                if(moments.size()==moments.capacity()) {
                    checkpoint("capacity");
                    throw std::length_error("sequence capacity " + std::to_string(capacity)
                        + " exhausted; need at least " + std::to_string(training+options_.holdout)
                        + (checkpoint_path.empty()
                            ? "; set checkpoint_path to save progress before restarting"
                            : "; checkpoint saved to " + checkpoint_path.string()
                                + "; restart with larger sequence_capacity and resume_checkpoint=true"));
                }
                {
                    timer_accum::guard measured(timing.sequence);
                    auto moment=moments.append();
                    if constexpr(std::same_as<K,Z32783>) {
                        block_wiedemann_detail::pack_projection<K>(current,packed_current,p.n,b,options_.threads);
                        block_wiedemann_detail::project_packed32783(packed_projection,packed_current,moment,
                            p.n,b,transposed_moments,options_.threads);
                    } else {
                        block_wiedemann_detail::project_block<K>(projection,current,moment,
                            p.n,b,transposed_moments,options_.threads);
                    }
                    ++sequence_stats.moments;
                    p.apply(current,next,b); current.swap(next);
                }
                if(options_.incremental_recurrence && moments.size()>options_.holdout) {
                    timer_accum::guard measured(timing.generator);
                    // Keep the newest holdout terms OUT of the recurrence state.
                    const auto available=moments.size()-options_.holdout;
                    state.process_up_to(available);
                    sequence_stats.recurrence_updates+=available-processed;
                    processed=available;
                }
                periodic("krylov");
                report(moments.size(),std::max(initial_training,training)+options_.holdout);
            }
            checkpoint("krylov-complete");
            if(!options_.recurrence_sequence_path.empty())
                block_wiedemann_detail::save_recurrence_sequence(options_.recurrence_sequence_path,moments,training);
            bool found=false;
            {
                timer_accum::guard measured(timing.generator);
                if(!options_.incremental_recurrence) {
                    progress recurrence{options_.log,"block recurrence"};
                    recurrence(processed,training);
                    state.process_up_to(training,[&](auto done,auto total){
                        sequence_stats.recurrence_updates+=done-processed;processed=done;
                        periodic("recurrence");recurrence(done,total);
                    });
                    checkpoint("recurrence-complete");
                }
                if(!options_.incremental_recurrence || state.candidate_ready()) {
                    ++sequence_stats.validation_attempts;
                    progress validation{options_.log,"recurrence validation"};
                    found=state.generator(gen,[&](auto done,auto total){periodic("validation");validation(done,total);},options_.incremental_recurrence);
                }
            }
            if(found) {
                accepted=true;checkpoint("validation-complete");
                if(options_.log)*options_.log << "recurrence accepted training=" << training
                    << " moments=" << moments.size() << " holdout=" << options_.holdout << std::endl;
                return std::optional<Generator>(std::move(gen));
            }
            if(training>=limit)return std::nullopt;
            if(options_.incremental_recurrence) {
                // Always check at the old batch boundary too, so near-full-rank
                // examples never overshoot it simply because of the check interval.
                const auto boundary=training<initial_training ? initial_training : limit;
                training+=std::min(interval,boundary-training);
            } else {
                if(options_.log)*options_.log << "recurrence validation failed at training=" << training << std::endl;
                training=std::min(limit,training*2);
            }
        }
    }

    // Block analogue of the scalar solver's Horner reconstruction. The
    // accumulator contains only the requested output vectors, not b Krylov
    // vectors; in particular one RHS uses width one even for a wide generator.
    using Reconstruction=block_wiedemann_detail::reconstruction_state<K>;
    Block reconstruct(preconditioner& p,const Block& initial,const Generator& gen,
        const Block& weights,std::size_t b,std::size_t count,
        Reconstruction* restored=nullptr,
        const std::function<void(const Reconstruction&)>& checkpoint={}) const {
        std::size_t degree=0;
        for(std::size_t row=0;row<gen.size();++row)degree=std::max(degree,gen.degree(row));
        progress report{options_.log,"Horner reconstruction"};
        std::optional<Reconstruction> fresh;
        if(!restored)fresh.emplace(p.n,count,degree);
        auto& state=restored ? *restored : *fresh;
        if(state.rows!=p.n || state.columns!=count || state.degree!=degree)
            throw std::runtime_error("reconstruction state dimensions mismatch");
        auto& acc=state.accumulator;
        report(degree-state.remaining,degree);
        Block next(p.n*count),combination(b*count);
        auto saved=std::chrono::steady_clock::now();
        if(checkpoint)checkpoint(state);
        while(state.remaining) {
            const auto power=state.remaining-1;
            if(power+1<degree){p.apply(acc,next,count);acc.swap(next);}
            std::fill(combination.begin(),combination.end(),K{});
            for(std::size_t j=0;j<b;++j) {
                const auto d=gen.degree(j);
                if(power>=d)continue;
                for(std::size_t i=0;i<b;++i)
                    for(std::size_t k=0;k<count;++k)
                        combination[i*count+k]+=gen.coefficient(j,d-1-power)[i]*weights[j*count+k];
            }
            for(std::size_t r=0;r<p.n;++r)
                for(std::size_t i=0;i<b;++i)
                    for(std::size_t k=0;k<count;++k)
                        acc[r*count+k]+=initial[r*b+i]*combination[i*count+k];
            state.remaining=power;
            report(degree-power,degree);
            const auto now=std::chrono::steady_clock::now();
            if(checkpoint && (!power || std::chrono::duration<double>(now-saved).count()>=options_.reconstruction_checkpoint_seconds)) {
                checkpoint(state);saved=std::chrono::steady_clock::now();
            }
        }
        return std::move(acc);
    }
    static std::optional<Block> solve_dense(Block a,Block rhs,std::size_t b) {
        for(std::size_t c=0;c<b;++c) {
            auto r=c;while(r<b && a[r*b+c]==K{})++r;
            if(r==b)return std::nullopt;
            for(std::size_t j=0;j<b;++j)std::swap(a[r*b+j],a[c*b+j]);
            std::swap(rhs[r],rhs[c]);
            auto inverse=a[c*b+c].inv();
            for(std::size_t j=c;j<b;++j)a[c*b+j]*=inverse;
            rhs[c]*=inverse;
            for(std::size_t i=0;i<b;++i)if(i!=c) {
                auto factor=a[i*b+c];
                for(std::size_t j=c;j<b;++j)a[i*b+j]-=factor*a[c*b+j];
                rhs[i]-=factor*rhs[c];
            }
        }
        return rhs;
    }
public:
    block_wiedemann_solver(std::size_t rows,std::size_t cols,Apply apply,Apply transpose,options config={})
        :rows_(rows),cols_(cols),apply_(std::move(apply)),transpose_(std::move(transpose)),options_(config) {
        if(!config.block_size || !config.trials || !config.holdout || !config.recurrence_check_interval || config.threads<1
            || !std::isfinite(config.checkpoint_seconds) || config.checkpoint_seconds<0
            || !std::isfinite(config.reconstruction_checkpoint_seconds) || config.reconstruction_checkpoint_seconds<0)
            throw std::invalid_argument("block size, trials, holdout and threads must be positive");
        static_assert(K::characteristic()!=0,"Wiedemann requires a finite field");
    }
    // Direct rank/nullspace mode for a supplied square operator. No extra Gram
    // product. The rank estimate assumes generic left diagonal scaling makes
    // zero semisimple (as for diagonally symmetrizable graph Laplacians).
    static block_wiedemann_solver from_square_operator(std::size_t n,Apply apply,options config={}) {
        if(!apply)throw std::invalid_argument("missing square operator");
        block_wiedemann_solver solver(n,n,std::move(apply),{},config);
        solver.direct_square_=true;
        return solver;
    }
    // Borrow the original CSC; store its transpose once, as the scalar solver
    // does. Both applications gather into independent output rows.
    template <class Matrix> requires requires(const Matrix& m) { m.rows;m.offsets;m.row_indices;m.coefficients; }
    explicit block_wiedemann_solver(const Matrix& matrix,options config={})
        :block_wiedemann_solver(matrix.rows,matrix.columns(),
            [transposed=std::make_shared<signed_sparse_transpose>(matrix),workers=config.threads]
            (auto input,auto output,std::size_t b) {
                evaluate_signed_transpose<K>(*transposed,input,output,b,workers);
            },
            [&matrix,workers=config.threads](auto input,auto output,std::size_t b) {
                evaluate_signed_transpose<K>(matrix,input,output,b,workers);
            },config) {}

    explicit block_wiedemann_solver(const compressed_sparse_matrix<K>& matrix,options config={})
        :block_wiedemann_solver(matrix.image_dim(),matrix.domain_dim(),
            [transposed=std::make_shared<compressed_sparse_matrix<K>>(matrix.transpose()),workers=config.threads]
            (auto input,auto output,std::size_t b) {
                evaluate_field_transpose(*transposed,input,output,b,workers);
            },
            [&matrix,workers=config.threads](auto input,auto output,std::size_t b) {
                evaluate_field_transpose(matrix,input,output,b,workers);
            },config) {}
private:
    static void evaluate_field_transpose(const compressed_sparse_matrix<K>& matrix,
        std::span<const K> input,std::span<K> output,std::size_t b,int threads) {
        const bool parallel=matrix.rows_and_coeffs_.size()*b>=32768;
        (void)parallel; (void)threads;
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
        for(std::size_t c=0;c<matrix.domain_dim();++c) {
            auto* row=output.data()+c*b;
            std::fill_n(row,b,K{});
            for(const auto& term:matrix.get_column(c))
                for(std::size_t j=0;j<b;++j)
                    row[j]+=input[term.getValue()*b+j]*term.getCoefficient();
        }
    }
private:
    // The legacy LIL constructor also owns the original compressed copy.
    std::shared_ptr<compressed_sparse_matrix<K>> owned_;
    block_wiedemann_solver(std::shared_ptr<compressed_sparse_matrix<K>> matrix,options config)
        :block_wiedemann_solver(*matrix,config) { owned_=std::move(matrix); }
public:
    block_wiedemann_solver(const lil_matrix<K>& matrix,std::size_t block_size)
        :block_wiedemann_solver(std::make_shared<compressed_sparse_matrix<K>>(matrix.to_compressed_sparse_matrix()),
            options{block_size,1,8,17}) {}
    std::size_t image_dim() const { return rows_; }
    std::size_t domain_dim() const { return cols_; }

    rank_result rank() const {
        std::optional<kernel_start> cached;
        if(options_.resume_checkpoint && !options_.checkpoint_path.empty()
            && std::filesystem::exists(options_.checkpoint_path.string()+".rank"))
            return load_rank(cached);
        return compute_rank(options_.checkpoint_path.empty() ? nullptr : &cached);
    }
private:
    void save_rank(const rank_result& result,const kernel_start& state) const {
        block_wiedemann_detail::checkpoint_output out(options_.checkpoint_path.string()+".rank");
        const auto b=state.gen.size();
        std::size_t degree=0;for(std::size_t r=0;r<b;++r)degree=std::max(degree,state.gen.degree(r));
        out.section("rank-handoff",[&] {
            out.word(1);out.word(rows_);out.word(cols_);out.word(b);
            out.word(options_.seed);out.word(options_.trials);out.word(direct_square_);
            out.word(options_.holdout);out.word(options_.incremental_recurrence);
            out.word(options_.recurrence_check_interval);
            out.word(result.rank);out.word(result.nullity);out.word(degree);
            out.integers<std::size_t>("trial-ranks",result.trial_ranks);
        });
        out.template block<K>("left-diagonal",state.p.left,state.p.left.size(),1);
        out.template block<K>("right-diagonal",state.p.right,state.p.right.size(),1);
        out.template block<K>("starting-vectors",state.v,cols_,b);
        out.template block<K>("initial",state.initial,cols_,b);
        out.template block<K>("projection",state.projection,cols_,b);
        state.gen.save(out);out.finish();
    }
    rank_result load_rank(std::optional<kernel_start>& cached) const {
        block_wiedemann_detail::checkpoint_input in(options_.checkpoint_path.string()+".rank");
        rank_result result;const auto b=std::min(options_.block_size,cols_);
        std::size_t degree=0;
        in.section("rank-handoff",[&] {
            in.expect(1);in.expect(rows_);in.expect(cols_);in.expect(b);
            in.expect(options_.seed);in.expect(options_.trials);in.expect(direct_square_);
            in.expect(options_.holdout);in.expect(options_.incremental_recurrence);
            in.expect(options_.recurrence_check_interval);
            result.rank=in.word();result.nullity=in.word();degree=in.word();
            if(!b || result.rank>std::min(rows_,cols_) || result.nullity!=cols_-result.rank || degree>cols_+1)
                throw std::runtime_error("invalid rank handoff dimensions");
            result.trial_ranks.resize(options_.trials);
            in.integers<std::size_t>("trial-ranks",result.trial_ranks);
            if(*std::max_element(result.trial_ranks.begin(),result.trial_ranks.end())!=result.rank)
                throw std::runtime_error("invalid rank handoff trials");
        });
        std::mt19937_64 rng(options_.seed);auto p=precondition(rng,true);
        in.template block<K>("left-diagonal",p.left,p.left.size(),1);
        in.template block<K>("right-diagonal",p.right,p.right.size(),1);
        Block v(cols_*b),initial(cols_*b),projection(cols_*b);
        in.template block<K>("starting-vectors",v,cols_,b);
        in.template block<K>("initial",initial,cols_,b);
        in.template block<K>("projection",projection,cols_,b);
        Generator gen(b,degree);gen.load(in);in.finish();
        std::size_t rank=0;for(std::size_t r=0;r<b;++r)rank+=gen.degree(r);
        if(rank!=result.rank)throw std::runtime_error("rank handoff recurrence degree mismatch");
        Block check(cols_*b);p.apply(v,check,b);
        if(check!=initial)throw std::runtime_error("rank handoff operator mismatch");
        cached.emplace(kernel_start{std::move(p),std::move(v),std::move(gen),std::move(initial),std::move(projection)});
        if(options_.log)*options_.log << "rank handoff loaded rank=" << result.rank << " nullity=" << result.nullity << std::endl;
        return result;
    }
    rank_result compute_rank(std::optional<kernel_start>* cached) const {
        rank_result result;
        result.nullity=cols_;
        if(!rows_ || !cols_) {result.probabilistic=false;return result;}
        std::mt19937_64 rng(options_.seed);
        for(std::size_t trial=0;trial<options_.trials;++trial) {
            if(options_.log)*options_.log << "rank trial " << trial+1 << "/" << options_.trials
                << " rows=" << rows_ << " columns=" << cols_ << std::endl;
            // A resumable rank uses the domain operator; rank-only can use the smaller one.
            // Both independent diagonal scalings
            // are needed: AA^T alone can lose rank over finite fields.
            auto p=precondition(rng,cached || rows_>cols_);
            const auto b=std::min(options_.block_size,p.n);
            Block u(p.n*b),v(p.n*b),initial(p.n*b);
            for(auto& x:u)x=K::sample(rng);
            for(auto& x:v)x=K::sample(rng);
            // Starting at B*V removes the semisimple zero part. The generic
            // preconditioned operator has rank equal to its nonzero McMillan degree.
            p.apply(v,initial,b);
            auto gen=sequence_generator(p,initial,u,b,cached!=nullptr,true,trial ? ".trial-"+std::to_string(trial+1) : "");
            if(!gen)throw std::runtime_error("block recurrence failed held-out validation; retry with another seed");
            std::size_t rank=0;for(std::size_t row=0;row<gen->size();++row)rank+=gen->degree(row);
            if(rank>p.n)throw std::runtime_error("invalid block generator degree");
            if(options_.log)*options_.log << "trial rank=" << rank << " nullity=" << cols_-rank << std::endl;
            result.trial_ranks.push_back(rank);
            if(cached && (!*cached || rank>result.rank))
                cached->emplace(kernel_start{std::move(p),std::move(v),std::move(*gen),std::move(initial),std::move(u)});
            result.rank=std::max(result.rank,rank);
        }
        result.nullity=cols_-result.rank;
        if(cached && *cached && !options_.checkpoint_path.empty())save_rank(result,**cached);
        return result;
    }

public:
    // Right nullspace of the original rectangular operator. No dense matrix,
    // stored transpose, or Gram matrix is constructed. Output costs O(cols*h).
    nullspace_result nullspace() const {
        nullspace_result result;
        std::optional<kernel_start> cached;
        const auto reconstruction_path=options_.checkpoint_path.string()+".reconstruction";
        const auto b=std::min(options_.block_size,cols_);
        std::vector<std::size_t> pivots;
        std::mt19937_64 rng(options_.seed);
        std::size_t first_attempt=0;
        std::optional<preconditioner> restored_p;
        Block restored_v,restored_initial,restored_weights;
        std::optional<Generator> restored_gen;
        std::optional<Reconstruction> restored_state;
        if(options_.resume_checkpoint && !options_.checkpoint_path.empty()
            && std::filesystem::exists(reconstruction_path)) {
            block_wiedemann_detail::checkpoint_input in(reconstruction_path);
            std::size_t basis_size=0,degree=0,count=0;
            in.section("nullspace-reconstruction",[&] {
                in.expect(1);in.expect(rows_);in.expect(cols_);in.expect(b);
                in.expect(options_.seed);in.expect(options_.trials);in.expect(direct_square_);
                in.expect(options_.holdout);in.expect(options_.incremental_recurrence);in.expect(options_.recurrence_check_interval);
                first_attempt=in.word();result.rank_estimate.rank=in.word();
                result.rank_estimate.nullity=in.word();basis_size=in.word();
                degree=in.word();count=in.word();
                if(result.rank_estimate.rank>cols_ || result.rank_estimate.nullity!=cols_-result.rank_estimate.rank
                    || basis_size>=result.rank_estimate.nullity || !count || count!=std::min(b,result.rank_estimate.nullity-basis_size)
                    || degree>cols_ || first_attempt>=(result.rank_estimate.nullity+b-1)/b+options_.trials)
                    throw std::runtime_error("invalid nullspace reconstruction metadata");
                const auto trials=in.word();
                if(trials!=options_.trials)throw std::runtime_error("invalid saved rank trials");
                result.rank_estimate.trial_ranks.resize(trials);
                in.integers<std::size_t>("trial-ranks",result.rank_estimate.trial_ranks);
                const auto size=in.word();
                if(size>100000)throw std::runtime_error("invalid random state length");
                std::string random_state(size,' ');in.bytes(random_state.data(),size);
                std::istringstream stream(random_state);
                if(!(stream>>rng))throw std::runtime_error("invalid random generator state");
            });
            pivots.resize(basis_size);in.integers<std::size_t>("basis-pivots",pivots);
            for(auto pivot:pivots)if(pivot>=cols_)throw std::runtime_error("invalid saved pivot");
            result.basis.resize(basis_size);
            for(auto& vector:result.basis) {
                vector.resize(cols_);in.template block<K>("verified-basis-vector",vector,cols_,1);
            }
            // Reattach callbacks; saved diagonals replace the dummy generated ones.
            auto dummy_rng=rng;restored_p=precondition(dummy_rng,true);
            auto& p=*restored_p;
            in.template block<K>("left-diagonal",p.left,p.left.size(),1);
            in.template block<K>("right-diagonal",p.right,p.right.size(),1);
            restored_v.resize(cols_*b);restored_initial.resize(cols_*b);restored_weights.resize(b*count);
            in.template block<K>("starting-vectors",restored_v,cols_,b);
            in.template block<K>("initial",restored_initial,cols_,b);
            in.template block<K>("reconstruction-weights",restored_weights,b,count);
            restored_gen.emplace(b,degree);restored_gen->load(in);
            restored_state.emplace(cols_,count,degree);restored_state->load(in);in.finish();
            Block check(cols_*b);p.apply(restored_v,check,b);
            if(check!=restored_initial)throw std::runtime_error("reconstruction checkpoint operator mismatch");
            if(options_.log)*options_.log << "reconstruction resumed remaining=" << restored_state->remaining << std::endl;
        } else if(options_.resume_checkpoint && !options_.checkpoint_path.empty()
            && std::filesystem::exists(options_.checkpoint_path.string()+".rank"))
            result.rank_estimate=load_rank(cached);
        else result.rank_estimate=compute_rank(&cached);
        const auto target=result.rank_estimate.nullity;
        if(!rows_) {
            for(std::size_t i=0;i<cols_;++i) {
                Block e(cols_);e[i]=1;result.basis.push_back(std::move(e));
            }
            result.complete=true;
            return result;
        }
        if(!target) {result.complete=true;return result;}
        const auto attempts=(target+b-1)/b+options_.trials;
        for(std::size_t attempt=first_attempt;attempt<attempts && result.basis.size()<target;++attempt) {
            if(options_.log)*options_.log << "nullspace attempt " << attempt+1 << "/" << attempts
                << " verified=" << result.basis.size() << "/" << target
                << " cached_recurrence=" << bool(cached) << std::endl;
            // P = D_col A^T D_row A acts on the domain. Project random V
            // along im(P) into ker(P), then check the result against A.
            auto p=restored_p ? std::move(*restored_p) : cached ? std::move(cached->p) : precondition(rng,true);
            Block v=restored_state ? std::move(restored_v) : cached ? std::move(cached->v) : Block(cols_*b);
            Block initial(cols_*b);
            std::optional<Generator> gen;
            if(restored_state) {
                restored_p.reset();initial=std::move(restored_initial);gen=std::move(restored_gen);
            } else if(cached) {
                gen=std::move(cached->gen);
                initial=std::move(cached->initial);
                cached.reset();
            } else {
                Block u(cols_*b);
                for(auto& x:v)x=K::sample(rng);
                for(auto& x:u)x=K::sample(rng);
                p.apply(v,initial,b);
                gen=sequence_generator(p,initial,u,b,true,true,".nullspace-"+std::to_string(attempt));
            }
            if(!gen)continue;
            timer_accum::guard measured(timing.reconstruction);
            Block constant(b*b),weights(b*b);
            for(std::size_t j=0;j<b;++j) {
                for(std::size_t i=0;i<b;++i)constant[i*b+j]=gen->coefficient(j,gen->degree(j))[i];
            }
            bool invertible=true;
            for(std::size_t j=0;j<b;++j) {
                Block unit(b);unit[j]=1;
                auto column=solve_dense(constant,unit,b);
                if(!column){invertible=false;break;}
                for(std::size_t i=0;i<b;++i)weights[i*b+j]=(*column)[i];
            }
            if(!invertible)continue;
            const auto count=std::min(b,target-result.basis.size());
            Block selected(b*count);
            for(std::size_t i=0;i<b;++i)
                for(std::size_t j=0;j<count;++j)selected[i*count+j]=weights[i*b+j];
            if(restored_state && selected!=restored_weights)
                throw std::runtime_error("reconstruction weights mismatch");
            const auto save_reconstruction=[&](const Reconstruction& state) {
                if(options_.checkpoint_path.empty())return;
                block_wiedemann_detail::checkpoint_output out(reconstruction_path);
                out.section("nullspace-reconstruction",[&] {
                    out.word(1);out.word(rows_);out.word(cols_);out.word(b);
                    out.word(options_.seed);out.word(options_.trials);out.word(direct_square_);
                    out.word(options_.holdout);out.word(options_.incremental_recurrence);out.word(options_.recurrence_check_interval);
                    out.word(attempt);out.word(result.rank_estimate.rank);out.word(target);
                    out.word(result.basis.size());out.word(state.degree);out.word(count);
                    out.word(result.rank_estimate.trial_ranks.size());
                    out.integers<std::size_t>("trial-ranks",result.rank_estimate.trial_ranks);
                    std::ostringstream stream;stream<<rng;const auto random_state=stream.str();
                    out.word(random_state.size());out.bytes(random_state.data(),random_state.size());
                });
                out.integers<std::size_t>("basis-pivots",pivots);
                for(const auto& vector:result.basis)out.template block<K>("verified-basis-vector",vector,cols_,1);
                out.template block<K>("left-diagonal",p.left,p.left.size(),1);
                out.template block<K>("right-diagonal",p.right,p.right.size(),1);
                out.template block<K>("starting-vectors",v,cols_,b);
                out.template block<K>("initial",initial,cols_,b);
                out.template block<K>("reconstruction-weights",selected,b,count);
                gen->save(out);state.save(out);out.finish();
                if(options_.log)*options_.log << "checkpoint saved phase=reconstruction remaining=" << state.remaining << std::endl;
                if(options_.on_checkpoint)options_.on_checkpoint("reconstruction");
            };
            auto candidates=reconstruct(p,initial,*gen,selected,b,count,
                restored_state ? &*restored_state : nullptr,save_reconstruction);
            restored_state.reset();
            for(std::size_t r=0;r<cols_;++r)
                for(std::size_t j=0;j<count;++j)candidates[r*count+j]+=v[r*b+j];
            Block check(rows_*count);
            apply_(candidates,check,count);
            for(std::size_t j=0;j<count && result.basis.size()<target;++j) {
                bool valid=true;
                for(std::size_t r=0;r<rows_;++r)valid &= check[r*count+j]==K{};
                if(!valid)continue;
                Block x(cols_);
                for(std::size_t r=0;r<cols_;++r)x[r]=candidates[r*count+j];
                // Keep one echelon basis, serving as both output and the
                // independence filter. Linear combinations stay in ker(A).
                for(std::size_t i=0;i<pivots.size();++i) {
                    const auto factor=x[pivots[i]];
                    for(std::size_t r=0;r<cols_;++r)x[r]-=factor*result.basis[i][r];
                }
                std::size_t pivot=0;while(pivot<cols_ && x[pivot]==K{})++pivot;
                if(pivot==cols_)continue;
                const auto inverse=x[pivot].inv();
                for(auto& value:x)value*=inverse;
                pivots.push_back(pivot);
                result.basis.push_back(std::move(x));
                if(options_.log)*options_.log << "verified independent vectors=" << result.basis.size()
                    << "/" << target << std::endl;
            }
        }
        result.complete=result.basis.size()==target;
        if(options_.log)*options_.log << "nullspace complete=" << result.complete << std::endl;
        return result;
    }

    // Returns only solutions verified against the ORIGINAL rectangular matrix.
    // No solution found is inconclusive, not a proof of inconsistency.
    std::vector<DomainVec> solve_MX_equals_y(const ImageVec& rhs) const {
        if(direct_square_)throw std::logic_error("direct square mode supports rank and nullspace only");
        std::vector<DomainVec> solutions;
        bool zero=true;for(std::size_t i=0;i<rows_;++i)zero &= rhs[i]==K{};
        if(zero) {solutions.push_back(std::make_unique<K[]>(cols_));return solutions;}
        if(!cols_)return solutions;
        std::mt19937_64 rng(options_.seed);
        for(std::size_t trial=0;trial<options_.trials;++trial) {
            if(options_.log)*options_.log << "solve trial " << trial+1 << "/" << options_.trials << std::endl;
            auto p=precondition(rng,false);
            const auto b=std::min(options_.block_size,p.n);
            Block u(p.n*b),v(p.n*b),initial(p.n*b);
            for(auto& x:u)x=K::sample(rng);
            for(auto& x:v)x=K::sample(rng);
            p.apply(v,initial,b);
            for(std::size_t i=0;i<p.n;++i)initial[i*b]=p.left[i]*rhs[i];
            // Transposing the moments produces a right generator for initial.
            auto gen=sequence_generator(p,initial,u,b,true);
            if(!gen)continue;
            Block constant(b*b),unit(b);unit[0]=1;
            for(std::size_t j=0;j<b;++j) {
                for(std::size_t i=0;i<b;++i)constant[i*b+j]=gen->coefficient(j,gen->degree(j))[i];
            }
            auto weights=solve_dense(constant,unit,b);
            if(!weights)continue;
            auto z=reconstruct(p,initial,*gen,*weights,b,1);
            for(auto& value:z)value=-value;
            Block x(cols_),check(rows_);
            transpose_(z,x,1);
            for(std::size_t i=0;i<cols_;++i)x[i]*=p.right[i];
            apply_(x,check,1);
            bool valid=true;for(std::size_t i=0;i<rows_;++i)valid &= check[i]==rhs[i];
            if(valid) {
                auto out=std::make_unique<K[]>(cols_);std::copy(x.begin(),x.end(),out.get());
                solutions.push_back(std::move(out));
                break;
            }
        }
        return solutions;
    }
};
} // namespace VectorSpace
