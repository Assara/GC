#pragma once
#include <type_traits>
#include "NaturalAdjoint.hpp"

namespace GraphHomology {
// Evaluate (S_down C_down + C_up S_up)^T, without internal diagonal weights.
// Reduction happens only after both two-step branches have been accumulated.
// Borrows matrices; owns reusable workspace and is not reentrant.
template<VectorSpace::Field K, class Accumulator = GraphAccumulator>
class NaturalComposition {
    static_assert(std::is_same_v<Accumulator,std::int32_t> || std::is_same_v<Accumulator,std::int64_t>);
    const ContractionMatrix& down_;
    const ContractionMatrix& up_;
    const NaturalAdjointMatrix& split_down_;
    const NaturalAdjointMatrix& split_up_;
    std::size_t block_size_;
    int threads_;
    VectorSpace::OwnedArray<Accumulator> scratch_, result_;
    std::uint64_t path_bound_ = 0;

    static std::uint64_t magnitude(auto coefficient) {
        const auto x = std::int64_t(coefficient);
        return x < 0 ? std::uint64_t(-x) : std::uint64_t(x);
    }
    void certify() {
        // Saturate at one above the admissible path count: never overflow while checking.
        const auto limit = std::uint64_t(std::numeric_limits<Accumulator>::max()) /
            (std::uint64_t(K::characteristic())-1);
        const auto add_product = [limit](std::uint64_t sum, std::uint64_t a, std::uint64_t b) {
            if (sum > limit || (b && a > (limit-sum)/b)) return limit+1;
            return sum+a*b;
        };
        const auto weights = [&](const auto& matrix) {
            VectorSpace::OwnedArray<std::uint64_t> w(matrix.columns());
            for (std::size_t c=0;c<matrix.columns();++c) {
                for(auto i=matrix.offsets[c];i<matrix.offsets[c+1];++i)
                    w[c] = add_product(w[c],magnitude(matrix.coefficients[i]),1);
                path_bound_ = std::max(path_bound_,w[c]);
            }
            return w;
        };
        const auto lower = weights(split_down_);
        const auto upper = weights(up_);
        for(std::size_t c=0;c<down_.columns();++c) {
            std::uint64_t w=0;
            for(auto i=down_.offsets[c];i<down_.offsets[c+1];++i)
                w=add_product(w,magnitude(down_.coefficients[i]),lower[down_.row_indices[i]]);
            for(auto i=split_up_.offsets[c];i<split_up_.offsets[c+1];++i)
                w=add_product(w,magnitude(split_up_.coefficients[i]),upper[split_up_.row_indices[i]]);
            path_bound_=std::max(path_bound_,w);
        }
        if(path_bound_>limit)
            throw std::overflow_error(sizeof(Accumulator)==4
                ? "natural composition exceeds int32 certificate; use GraphAccumulator=int64_t"
                : "natural composition exceeds int64 certificate; intermediate reductions required");
    }
    template<class Matrix, class Input>
    void gather(const Matrix& matrix, const Input* input, Accumulator* output, bool add, std::size_t b) {
        const bool parallel=matrix.nonzeros()*b>=32768;
#pragma omp parallel for schedule(static) num_threads(threads_) if(parallel)
        for(std::size_t c=0;c<matrix.columns();++c) {
            auto* row=output+c*b;
            if(!add) std::fill_n(row,b,Accumulator{});
            for(auto i=matrix.offsets[c];i<matrix.offsets[c+1];++i) {
                const auto coefficient=Accumulator(matrix.coefficients[i]);
                const auto* values=input+std::size_t(matrix.row_indices[i])*b;
                for(std::size_t j=0;j<b;++j) {
                    const auto value = [&] {
                        if constexpr(std::is_same_v<Input,K>) return Accumulator(values[j].value());
                        else return values[j];
                    }();
                    row[j] += coefficient*value;
                }
            }
        }
    }
public:
    NaturalComposition(const ContractionMatrix& down, const ContractionMatrix& up,
        const NaturalAdjoints& adjoints, std::size_t block_size=8, int threads=8)
        : down_(down),up_(up),split_down_(adjoints.down),split_up_(adjoints.up),
          block_size_(block_size),threads_(threads) {
        if(!block_size || threads<1) throw std::invalid_argument("invalid block size or threads");
        if(down.columns()!=up.rows || split_down_.rows!=down.columns() ||
            split_down_.columns()!=down.rows || split_up_.rows!=up.columns() ||
            split_up_.columns()!=up.rows) throw std::invalid_argument("natural composition dimensions mismatch");
        certify();
        const auto max_rows=std::max({std::size_t(down.rows),up.columns(),down.columns()});
        if(max_rows>std::numeric_limits<std::size_t>::max()/block_size)
            throw std::length_error("natural composition workspace size overflow");
        scratch_=VectorSpace::OwnedArray<Accumulator>(std::max(std::size_t(down.rows),up.columns())*block_size);
        result_=VectorSpace::OwnedArray<Accumulator>(down.columns()*block_size);
    }
    std::uint64_t path_bound() const { return path_bound_; }
    std::size_t workspace_bytes() const { return (scratch_.size()+result_.size())*sizeof(Accumulator); }
    void apply(std::span<const K> input, std::span<K> output) {
        apply(input,output,block_size_);
    }
    // Narrow reconstruction blocks reuse the same allocation as Krylov blocks.
    void apply(std::span<const K> input, std::span<K> output, std::size_t b) {
        if(!b || b>block_size_ || input.size()!=down_.columns()*b || output.size()!=input.size())
            throw std::invalid_argument("natural composition input/output size mismatch");
        gather(split_down_,input.data(),scratch_.data(),false,b);
        gather(down_,scratch_.data(),result_.data(),false,b);
        gather(up_,input.data(),scratch_.data(),false,b);
        gather(split_up_,scratch_.data(),result_.data(),true,b);
        const bool parallel=output.size()>=32768;
#pragma omp parallel for schedule(static) num_threads(threads_) if(parallel)
        for(std::size_t i=0;i<output.size();++i) output[i]=K(result_[i]);
    }
};
}
