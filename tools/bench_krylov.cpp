#include <chrono>
#include <iostream>
#include <random>
#include "GraphHomology/NaturalComposition.hpp"
#include "GraphHomology/NaturalRepresentatives.hpp"
#include "VectorSpace/packed_projection.hpp"
#ifndef GC_HOMOLOGY_LOOP
#define GC_HOMOLOGY_LOOP 10
#endif
#ifndef GC_HOMOLOGY_VERTICES
#define GC_HOMOLOGY_VERTICES 14
#endif
using K=fieldType;
using Clock=std::chrono::steady_clock;
int main(int argc,char** argv) {
 try {
  if(argc!=2)throw std::invalid_argument("usage: bench_krylov GRAPH_DIRECTORY");
  constexpr std::size_t b=8,steps=32;constexpr int threads=8;
  auto start=Clock::now();
  std::cerr<<"building graph matrices"<<std::endl;
  GraphHomology::ContractionWindow<GC_HOMOLOGY_LOOP,GC_HOMOLOGY_VERTICES,GraphHomology::Parity::even> window(argv[1]);
  std::cerr<<"building natural adjoints"<<std::endl;
  GraphHomology::NaturalAdjoints adjoints(window);
  GraphHomology::check_automorphism_units<K>(adjoints);
  GraphHomology::NaturalComposition<K> composition(window.down,window.up,adjoints,b,threads);
  const auto n=window.middle.size();
  std::vector<K> diagonal(n),projection(n*b),initial(n*b),current(n*b),next(n*b),packed_projection(n*b),packed_current(n*b),moment(b*b);
  std::mt19937_64 rng(17);
  for(auto& x:diagonal)do{x=K::sample(rng);}while(x==K{});
  for(auto& x:projection)x=K::sample(rng);
  for(auto& x:current)x=K::sample(rng);
  auto apply=[&](const auto& in,auto& out) {
   composition.apply(in,out,b);
#pragma omp parallel for schedule(static) num_threads(threads)
   for(std::size_t i=0;i<n;++i)for(std::size_t j=0;j<b;++j)out[i*b+j]*=diagonal[i];
  };
  apply(current,initial);current=initial;
  using namespace VectorSpace::block_wiedemann_detail;
  pack_projection<K>(projection,packed_projection,n,b,threads);
  auto step=[&]{
   pack_projection<K>(current,packed_current,n,b,threads);
   project_packed32783(packed_projection,packed_current,moment,n,b,true,threads);
   apply(current,next);current.swap(next);
  };
  std::cerr<<"setup_seconds="<<std::chrono::duration<double>(Clock::now()-start).count()
   <<" dimension="<<n<<" block="<<b<<" threads="<<threads<<" path_bound="<<composition.path_bound()<<std::endl;
  for(int i=0;i<4;++i)step();
  std::cout<<"repeat,steps,dimension,total_seconds,seconds_per_step,checksum\n";
  K checksum{};
  for(int repeat=0;repeat<3;++repeat){
   current=initial;const auto begin=Clock::now();
   for(std::size_t i=0;i<steps;++i){step();for(auto value:moment)checksum+=value;}
   const auto elapsed=std::chrono::duration<double>(Clock::now()-begin).count();
   std::cout<<repeat<<','<<steps<<','<<n<<','<<elapsed<<','<<elapsed/steps<<','<<checksum.value()<<std::endl;
  }
 }catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
