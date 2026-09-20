#include <chrono>
#include <iostream>
#include "types.hpp"
#include "VectorSpace/divide_conquer_generator.hpp"
#include "VectorSpace/recurrence_sequence.hpp"
using namespace VectorSpace::block_wiedemann_detail;
using K=fieldType;
int main(int argc,char** argv) {
    try {
        if(argc<3 || argc>5)throw std::invalid_argument("usage: find_recurrence SEQUENCE OUTPUT [--resume [THREADS=8]]");
        const bool resume=argc>=4;
        if(resume && std::string(argv[3])!="--resume")throw std::invalid_argument("expected --resume");
        const int threads=argc==5?std::stoi(argv[4]):8;
        auto input=recurrence_sequence<K>::load(argv[1]);
        divide_conquer_generator<K> finder(input.moments,input.training,threads);
        const std::filesystem::path checkpoint=std::string(argv[2])+".state";
        if(resume) {
            VectorSpace::serialization::archive_input in(checkpoint);finder.load(in);in.finish();
            std::cout<<"resumed processed="<<finder.processed_terms()<<"/"<<input.training<<std::endl;
        }
        const auto save=[&] {
            VectorSpace::serialization::archive_output out(checkpoint);finder.save(out);out.finish();
        };
        save();
        std::cout<<"workspace_bytes="<<finder.allocated_bytes()<<" checkpoint="<<checkpoint<<std::endl;
        auto saved=std::chrono::steady_clock::now(),reported=saved;
        while(!finder.complete()) {
            finder.step();const auto now=std::chrono::steady_clock::now();
            if(now-reported>=std::chrono::seconds(5) || finder.complete()) {
                std::cout<<"PM-Basis leaf terms "<<finder.processed_terms()<<"/"<<input.training
                    <<" construction_complete="<<finder.complete()<<std::endl;reported=now;
            }
            if(now-saved>=std::chrono::seconds(60) || finder.complete()) {
                save();saved=std::chrono::steady_clock::now();
            }
        }
        packed_generator<K> result(input.moments.block_size(),input.training/2);
        if(!finder.generator(result))throw std::runtime_error("recurrence failed validation; completed construction checkpoint retained");
        VectorSpace::serialization::archive_output out(argv[2]);result.save(out);out.finish();
        std::cout<<"validated recurrence saved to "<<argv[2]<<std::endl;
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
