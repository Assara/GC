#include <chrono>
#include <iostream>
#include "types.hpp"
#include "VectorSpace/divide_conquer_generator.hpp"
#include "VectorSpace/recurrence_sequence.hpp"
using namespace VectorSpace::block_wiedemann_detail;
using K=fieldType;
using Clock=std::chrono::steady_clock;
int main(int argc,char** argv) {
    try {
        if(argc<2 || argc>6)throw std::invalid_argument("usage: bench_recurrence_finders SEQUENCE [REPEATS=3 [THREADS=8 [LEAF=16 [MULTIPLY_LEAF=8]]]]");
        const auto repeats=argc>2?std::stoul(argv[2]):3,threads=argc>3?std::stoul(argv[3]):8,
            leaf=argc>4?std::stoul(argv[4]):16,multiply_leaf=argc>5?std::stoul(argv[5]):8;
        if(!repeats || !threads || !leaf || !multiply_leaf)throw std::invalid_argument("options must be positive");
        auto input=recurrence_sequence<K>::load(argv[1]);
        const auto b=input.moments.block_size(),training=input.training;
        std::cout<<"finder,repeat,block,training,holdout,threads,setup_seconds,construction_seconds,validation_seconds,workspace_bytes,accepted,degree_sum\n";
        for(std::size_t repeat=0;repeat<repeats;++repeat) {
            std::optional<packed_generator<K>> reference_result;
            std::optional<bool> reference_accepted;
            for(std::size_t pass=0;pass<2;++pass) {
                const bool divide=(repeat+pass)%2;
                const auto setup=Clock::now();
                std::optional<minimal_generator_state<K>> baseline;
                std::optional<divide_conquer_generator<K>> candidate;
                if(divide)candidate.emplace(input.moments,training,threads,leaf,multiply_leaf);
                else baseline.emplace(input.moments,b,threads);
                packed_generator<K> result(b,training/2);
                const auto begin=Clock::now();
                if(divide)candidate->process();else baseline->process_up_to(training);
                const auto built=Clock::now();
                const bool accepted=divide?candidate->generator(result):baseline->generator(result);
                const auto checked=Clock::now();
                const auto bytes=divide?candidate->allocated_bytes():baseline->allocated_bytes();
                std::size_t degree=0;
                if(accepted)for(std::size_t r=0;r<b;++r)degree+=result.degree(r);
                if(reference_accepted) {
                    if(*reference_accepted!=accepted)throw std::runtime_error("finders disagree on acceptance");
                    if(accepted)for(std::size_t r=0;r<b;++r)
                        if(result.degree(r)!=reference_result->degree(r) || !std::ranges::equal(result.row(r),reference_result->row(r)))
                            throw std::runtime_error("finders returned different recurrences");
                } else {reference_accepted=accepted;if(accepted)reference_result.emplace(std::move(result));}
                std::cout<<(divide?"divide-conquer":"reference")<<','<<repeat<<','<<b<<','<<training<<','
                    <<input.moments.size()-training<<','<<threads<<','
                    <<std::chrono::duration<double>(begin-setup).count()<<','
                    <<std::chrono::duration<double>(built-begin).count()<<','
                    <<std::chrono::duration<double>(checked-built).count()<<','<<bytes<<','<<accepted<<','<<degree<<std::endl;
            }
        }
    }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
