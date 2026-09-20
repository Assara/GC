#include <cassert>
#include <iostream>
#include <random>
#include "types.hpp"
#include "VectorSpace/divide_conquer_generator.hpp"
#include "VectorSpace/recurrence_sequence.hpp"
using namespace VectorSpace::block_wiedemann_detail;
using K=fieldType;
int main(int argc,char** argv) {
    std::size_t cases=0;
    for(std::size_t b:{1,2,3,4,8})for(std::size_t rank:{0,1,7,13,31})for(std::size_t leaf:{1,4,16}) {
        const auto training=2*rank+4,total=training+8;
        packed_moments<K> moments(b,total);
        std::mt19937_64 rng(17+rank+b);
        std::vector<K> u(rank*b),v(rank*b),powers(rank,K(1));
        for(auto& x:u)x=K::sample(rng);for(auto& x:v)x=K::sample(rng);
        for(std::size_t t=0;t<total;++t) {
            auto block=moments.append();
            for(std::size_t a=0;a<rank;++a) {
                for(std::size_t i=0;i<b;++i)for(std::size_t j=0;j<b;++j)block[i*b+j]+=u[a*b+i]*powers[a]*v[a*b+j];
                powers[a]*=K(a+1);
            }
        }
        minimal_generator_state<K> reference(moments,b);
        divide_conquer_generator<K> candidate(moments,training,b==8?2:1,leaf,b==8?8:2);
        reference.process_up_to(training);candidate.process();
        packed_generator<K> expected(b,training/2),actual(b,training/2);
        const auto accepted=reference.generator(expected);assert(candidate.generator(actual)==accepted);
        if(accepted) {
            std::size_t degree=0;
            for(std::size_t r=0;r<b;++r) {
                degree+=actual.degree(r);
                assert(actual.degree(r)==expected.degree(r));
                assert(std::ranges::equal(actual.row(r),expected.row(r)));
            }
            assert(degree==rank);
        }
        // Corrupt independent holdout data: neither finder may accept it.
        moments[total-1][0]+=K(1);
        assert(!reference.generator(expected));assert(!candidate.generator(actual));
        ++cases;
    }
    // Arbitrary sequences, including odd orders and singular discrepancies.
    for(std::size_t b:{1,2,3,8})for(std::size_t training:{3,9,17,33,65}) {
        packed_moments<K> moments(b,training+8);std::mt19937_64 rng(training+b);
        for(std::size_t t=0;t<training+8;++t) {
            auto block=moments.append();for(auto& x:block)x=K::sample(rng);
            if(t%3==0)std::fill(block.begin(),block.end(),K{});
        }
        minimal_generator_state<K> reference(moments,b);
        divide_conquer_generator<K> candidate(moments,training,2,4,2);
        reference.process_up_to(training);candidate.process();
        packed_generator<K> a(b,training/2),c(b,training/2);
        assert(reference.generator(a)==candidate.generator(c));
        const auto file=std::filesystem::path("build")/("recurrence_sequence_test_"+std::to_string(b)+"_"+std::to_string(training)+".bin");
        save_recurrence_sequence(file,moments,training);
        auto loaded=recurrence_sequence<K>::load(file);
        assert(loaded.training==training && loaded.moments.size()==moments.size());
        for(std::size_t t=0;t<moments.size();++t)assert(std::ranges::equal(loaded.moments[t],moments[t]));
        ++cases;
    }
    // Save/recreate at EVERY stable transition, including unfinished leaves,
    // both child branches, residual products, final composition and completion.
    for(std::size_t training:{3,17,34}) {
        packed_moments<K> moments(2,training+8);
        for(std::size_t t=0;t<training+8;++t) {
            auto term=moments.append();K power=1;
            for(std::size_t i=0;i<t;++i)power*=K(3);
            term[0]=power;term[3]=power;
        }
        using Finder=divide_conquer_generator<K>;
        Finder reference(moments,training,1,4,2);
        std::size_t expected_steps=0;
        while(!reference.complete()){reference.step();++expected_steps;}
        packed_generator<K> expected(2,training/2),actual(2,training/2);
        assert(reference.generator(expected));
        const auto path=std::filesystem::path("build")/("pm_resume_test_"+std::to_string(training)+".bin");
        {
            Finder initial(moments,training,1,4,2);
            VectorSpace::serialization::archive_output out(path);initial.save(out);out.finish();
        }
        std::size_t steps=0;
        for(;;) {
            Finder resumed(moments,training,2,4,2); // Thread count may change.
            const auto bytes=resumed.allocated_bytes();
            {VectorSpace::serialization::archive_input in(path);resumed.load(in);in.finish();}
            if(resumed.complete()) {
                assert(resumed.generator(actual));break;
            }
            resumed.step();++steps;
            assert(bytes==resumed.allocated_bytes());
            VectorSpace::serialization::archive_output out(path);resumed.save(out);out.finish();
        }
        assert(steps==expected_steps);
        for(std::size_t r=0;r<2;++r)assert(std::ranges::equal(expected.row(r),actual.row(r)));
        moments[moments.size()-1][0]+=K(1);
        bool rejected=false;
        try {
            Finder wrong(moments,training,1,4,2);
            VectorSpace::serialization::archive_input in(path);wrong.load(in);
        }catch(const std::runtime_error&){rejected=true;}
        assert(rejected);
        ++cases;
    }
    if(argc==2) {
        auto input=recurrence_sequence<K>::load(argv[1]);
        const auto b=input.moments.block_size(),training=input.training;
        minimal_generator_state<K> reference(input.moments,b);
        divide_conquer_generator<K> candidate(input.moments,training,2,1,2);
        reference.process_up_to(training);candidate.process();
        packed_generator<K> a(b,training/2),c(b,training/2);
        assert(reference.generator(a));assert(candidate.generator(c));
        for(std::size_t r=0;r<b;++r) {
            assert(a.degree(r)==c.degree(r));assert(std::ranges::equal(a.row(r),c.row(r)));
        }
        ++cases;
    }
    std::cout<<cases<<" recurrence finder comparisons and corrupted-holdout checks passed\n";
}
