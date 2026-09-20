#include <chrono>
#include <iomanip>
#include <iostream>
#include "GraphHomology/NaturalComposition.hpp"
#include "GraphHomology/NaturalRepresentatives.hpp"
#include "VectorSpace/block_wiedemann.hpp"
using K=fieldType;
using Solver=VectorSpace::block_wiedemann_solver<K>;

// Canonical coordinates for comparing subspaces, independent of extraction order.
auto canonical(std::vector<std::vector<K>> basis) {
    if(basis.empty())return basis;
    std::size_t rank=0;
    for(std::size_t c=0;c<basis[0].size() && rank<basis.size();++c) {
        auto pivot=rank;
        while(pivot<basis.size() && basis[pivot][c]==K{})++pivot;
        if(pivot==basis.size())continue;
        std::swap(basis[pivot],basis[rank]);
        const auto inv=basis[rank][c].inv();
        for(auto& value:basis[rank])value*=inv;
        for(std::size_t r=0;r<basis.size();++r)if(r!=rank) {
            const auto factor=basis[r][c];
            for(std::size_t j=0;j<basis[r].size();++j)basis[r][j]-=factor*basis[rank][j];
        }
        ++rank;
    }
    if(rank!=basis.size())throw std::runtime_error("dependent output vectors");
    return basis;
}
void report(const char* name,int repeat,std::uint64_t seed,bool incremental,
    const Solver& solver,double total,std::size_t rank,std::size_t nullity) {
    std::cout << name << ',' << repeat << ',' << seed << ',' << (incremental?"incremental":"batch")
        << ',' << solver.sequence_stats.moments << ',' << solver.sequence_stats.recurrence_updates
        << ',' << solver.sequence_stats.validation_attempts << ',' << solver.timing.gram.calls
        << ',' << solver.timing.sequence.seconds << ',' << solver.timing.generator.seconds
        << ',' << solver.timing.reconstruction.seconds << ',' << total << ',' << rank << ',' << nullity
        << std::endl;
}
template<GraphHomology::Parity P>
void graphs(const char* directory,int repeats) {
    GraphHomology::ContractionWindow<8,12,P> w(directory);
    GraphHomology::NaturalAdjoints adjoints(w);
    GraphHomology::NaturalComposition<K> composition(w.down,w.up,adjoints,8,8);
    for(int repeat=0;repeat<repeats;++repeat) {
        const auto seed=std::uint64_t(17+repeat);
        std::optional<std::vector<std::vector<K>>> expected;
        for(int pass=0;pass<2;++pass) {
            const bool incremental=(pass+repeat)%2;
            Solver::options options;
            options.seed=seed;options.incremental_recurrence=incremental;
            auto solver=Solver::from_square_operator(w.middle.size(),[&](auto in,auto out,auto b){
                composition.apply(in,out,b);
            },options);
            const auto started=std::chrono::steady_clock::now();
            auto result=solver.nullspace();
            const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
            if(!result.complete)throw std::runtime_error("incomplete graph nullspace");
            GraphHomology::convert_natural_representatives<K>(result.basis,w.down,adjoints);
            auto normalized=canonical(result.basis);
            if(expected && *expected!=normalized)throw std::runtime_error("batch/incremental nullspaces differ");
            expected=std::move(normalized);
            report(P==GraphHomology::Parity::even?"L8_V12_even":"L8_V12_odd",repeat,seed,incremental,
                solver,elapsed,result.rank_estimate.rank,result.basis.size());
        }
    }
}
void diagonal(int repeats,std::size_t rank) {
    constexpr std::size_t n=2048;
    for(int repeat=0;repeat<repeats;++repeat)for(int pass=0;pass<2;++pass) {
        const bool incremental=(pass+repeat)%2;
        Solver::options options;
        options.seed=17+repeat;options.incremental_recurrence=incremental;
        auto solver=Solver::from_square_operator(n,[&](auto in,auto out,auto b){
            for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<b;++j)
                out[i*b+j]=i<rank ? K(i+1)*in[i*b+j] : K{};
        },options);
        const auto started=std::chrono::steady_clock::now();
        const auto result=solver.rank();
        const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
        if(result.rank!=rank)throw std::runtime_error("incorrect diagonal rank");
        report(rank==64?"diagonal_rank64":"diagonal_rank2047",repeat,options.seed,incremental,
            solver,elapsed,result.rank,result.nullity);
    }
}
int main(int argc,char** argv) {
    if(argc<2 || argc>3)return 2;
    try {
        const int repeats=argc==3?std::stoi(argv[2]):3;
        if(repeats<1)throw std::invalid_argument("positive repeat count required");
        std::cout << std::setprecision(9)
            << "case,repeat,seed,mode,moments,training,validation_attempts,operator_applications,sequence_seconds,recurrence_seconds,reconstruction_seconds,total_seconds,rank,nullity\n";
        diagonal(repeats,64);
        diagonal(repeats,2047);
        graphs<GraphHomology::Parity::even>(argv[1],repeats);
        graphs<GraphHomology::Parity::odd>(argv[1],repeats);
    }catch(const std::exception& e){std::cerr << e.what() << std::endl;return 1;}
}
