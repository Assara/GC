#include "VectorSpace/modular_dot32783.hpp"
#include <chrono>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>
#include <array>
#include <string>
#include <omp.h>
using namespace VectorSpace::experimental;
using K=Z32783;
using Clock=std::chrono::steady_clock;
void require(bool ok){if(!ok)throw std::runtime_error("modular dot mismatch");}
template<class T> void check(const std::vector<T>& a,const std::vector<T>& b,std::size_t n,std::size_t sa,std::size_t sb,K expected){
 require(dot32783<DotMethod::integer>(a.data(),b.data(),n,sa,sb)==expected);
 require(dot32783<DotMethod::floating_integer_reduce>(a.data(),b.data(),n,sa,sb)==expected);
 require(dot32783<DotMethod::floating_reciprocal>(a.data(),b.data(),n,sa,sb)==expected);
}
void correctness(){
 std::mt19937_64 rng(71);
 for(std::size_t n:{0ul,1ul,13ul,68987ul,dot_chunk-1,dot_chunk,dot_chunk+1,2*dot_chunk+3}){
  std::vector<K>a(n),b(n);for(auto&x:a)x=K::sample(rng);for(auto&x:b)x=K::sample(rng);
  auto expected=dot32783<DotMethod::field>(a.data(),b.data(),n);
  check(a,b,n,1,1,expected);
  std::vector<double> da(n),db(n);for(std::size_t i=0;i<n;++i){da[i]=a[i].value();db[i]=b[i].value();}
  check(da,db,n,1,1,expected);
  std::fill(a.begin(),a.end(),K(32782));std::fill(b.begin(),b.end(),K(32782));check(a,b,n,1,1,K(n));
 }
 for(auto sa:{1ul,3ul,8ul})for(auto sb:{1ul,5ul,8ul}){
  std::size_t n=1003;std::vector<K>a(n*sa),b(n*sb);for(auto&x:a)x=K::sample(rng);for(auto&x:b)x=K::sample(rng);
  check(a,b,n,sa,sb,dot32783<DotMethod::field>(a.data(),b.data(),n,sa,sb));
 }
 const std::uint64_t bound=(1ULL<<51)-1;
 for(auto x:{0ULL,1ULL,32782ULL,32783ULL,32784ULL,(1ULL<<51)-1})
  require(reduce_double32783(double(x))==x%dot_prime);
 for(int i=0;i<100000;++i){auto x=rng()%bound;require(reduce_double32783(double(x))==x%dot_prime);auto multiple=(x/dot_prime)*dot_prime;
  for(int delta:{-1,0,1})if(multiple>0){auto y=multiple+delta;require(reduce_double32783(double(y))==y%dot_prime);}}
 std::cout<<"correctness passed: random, maximum residues, chunk boundaries, strides, reciprocal boundaries\n";
}
template<DotMethod method,class T>
void project(const std::vector<T>& a,const std::vector<T>& b,std::vector<K>& out,std::size_t n,std::size_t width,int threads){
#pragma omp parallel for schedule(static) num_threads(threads) if(threads>1)
 for(std::size_t k=0;k<width*width;++k)out[k]=dot32783<method>(a.data()+k/width,b.data()+k%width,n,width,width);
}
template<DotMethod method,class T>
void project_packed(const std::vector<T>& a,const std::vector<T>& b,std::vector<K>& out,std::size_t n,std::size_t width,int threads){
#pragma omp parallel for schedule(static) num_threads(threads) if(threads>1)
 for(std::size_t k=0;k<width*width;++k)out[k]=dot32783<method>(a.data()+(k/width)*n,b.data()+(k%width)*n,n);
}
void benchmark(std::size_t n,std::size_t width,int threads,bool packed_only=false){
 std::mt19937_64 rng(17);std::vector<K>a(n*width),b(n*width),out(width*width),expected(width*width);
 for(auto&x:a)x=K::sample(rng);for(auto&x:b)x=K::sample(rng);
 std::vector<double> da(a.size()),db(b.size());auto start=Clock::now();
 for(std::size_t i=0;i<a.size();++i){da[i]=a[i].value();db[i]=b[i].value();}
 double conversion=std::chrono::duration<double>(Clock::now()-start).count();
 project<DotMethod::field>(a,b,expected,n,width,1);
 std::vector<K> pa(a.size()),pb(b.size());std::vector<double> pda(a.size()),pdb(b.size());
 start=Clock::now();
 for(std::size_t lane=0;lane<width;++lane)for(std::size_t r=0;r<n;++r){pa[lane*n+r]=a[r*width+lane];pb[lane*n+r]=b[r*width+lane];}
 double pack_integer=std::chrono::duration<double>(Clock::now()-start).count();
 start=Clock::now();
 for(std::size_t lane=0;lane<width;++lane)for(std::size_t r=0;r<n;++r){pda[lane*n+r]=a[r*width+lane].value();pdb[lane*n+r]=b[r*width+lane].value();}
 double pack_double=std::chrono::duration<double>(Clock::now()-start).count();
 const char* names[]={"field","uint64","double_cast_uint64_mod","double_reciprocal","preconverted_double","packed_uint64","packed_double"};
 std::array<std::vector<double>,7> timings;
 for(int round=0;round<6;++round)for(int k=0;k<7;++k){int method=(round+k)%7;if(packed_only && method<5)continue;start=Clock::now();
  switch(method){
  case 0:project<DotMethod::field>(a,b,out,n,width,threads);break;
  case 1:project<DotMethod::integer>(a,b,out,n,width,threads);break;
  case 2:project<DotMethod::floating_integer_reduce>(a,b,out,n,width,threads);break;
  case 3:project<DotMethod::floating_reciprocal>(a,b,out,n,width,threads);break;
  case 4:project<DotMethod::floating_reciprocal>(da,db,out,n,width,threads);break;
  case 5:project_packed<DotMethod::integer>(pa,pb,out,n,width,threads);break;
  case 6:project_packed<DotMethod::floating_reciprocal>(pda,pdb,out,n,width,threads);break;}
  double seconds=std::chrono::duration<double>(Clock::now()-start).count();require(out==expected);
  if(round)timings[method].push_back(seconds);
 }
 for(int k=0;k<7;++k){auto&t=timings[k];if(t.empty())continue;std::sort(t.begin(),t.end());
 std::cout<<"n="<<n<<" width="<<width<<" threads="<<threads<<" method="<<names[k]<<" median_seconds="<<t[t.size()/2]<<" min="<<t.front()<<" max="<<t.back();
 if(k==5)std::cout<<" separate_packing_seconds="<<pack_integer;
 if(k==6)std::cout<<" separate_packing_seconds="<<pack_double;
 if(k==4)std::cout<<" separate_conversion_seconds="<<conversion;
 std::cout<<std::endl;}
}
int main(int argc,char**argv){try{
 correctness();if(argc>1 && std::string(argv[1])=="--test-only")return 0;
 // Width 8 is the actual strided projection. Width 1 models a contiguous
 // recurrence dot with the same number of scalar terms as L10 degree * 8.
 if(argc==5 && std::string(argv[1])=="--case") {
  const auto n=std::stoull(argv[2]),width=std::stoull(argv[3]);const int threads=std::stoi(argv[4]);
  if(!n || !width || threads<1)throw std::invalid_argument("positive dimensions and threads required");
  benchmark(n,width,threads);return 0;
 }
 bool packed_only=argc>1 && std::string(argv[1])=="--packed-only";
 for(int threads:{1,8}){benchmark(68987,8,threads,packed_only);benchmark(1703974,8,threads,packed_only);}
 if(!packed_only)benchmark(1703976,1,1);
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}}
