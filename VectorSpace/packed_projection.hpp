#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <span>
#include "Field.hpp"

namespace VectorSpace::block_wiedemann_detail {
// Lane-major copies make each dot contiguous. The fixed projection is packed
// once per sequence; the evolving block reuses one allocation at every step.
template<Field K>
void pack_projection(std::span<const K> in,std::span<K> out,
    std::size_t n,std::size_t b,int threads) {
    const bool parallel=threads>1 && n*b>=32768;
    (void)parallel;
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
    for(std::size_t lane=0;lane<b;++lane)
        for(std::size_t r=0;r<n;++r)out[lane*n+r]=in[r*b+lane];
}

inline Z32783 integer_dot32783(const Z32783* a,const Z32783* b,std::size_t n) {
    constexpr std::size_t chunk=1u<<20;
    constexpr std::uint64_t prime=32783;
    static_assert(chunk*(prime-1)*(prime-1)<(1ULL<<51));
    Z32783 result{};
    for(std::size_t first=0;first<n;first+=chunk) {
        const auto count=std::min(chunk,n-first);
        std::uint64_t sum=0;
#pragma omp simd reduction(+:sum)
        for(std::size_t r=0;r<count;++r)
            sum+=std::uint64_t(a[first+r].value())*b[first+r].value();
        result+=Z32783(sum%prime);
    }
    return result;
}
inline void project_packed32783(std::span<const Z32783> projection,
    std::span<const Z32783> current,std::span<Z32783> moment,
    std::size_t n,std::size_t b,bool transposed,int threads) {
    const bool parallel=threads>1 && n*b*b>=32768;
    (void)parallel;
#pragma omp parallel for collapse(2) schedule(static) num_threads(threads) if(parallel)
    for(std::size_t i=0;i<b;++i)
        for(std::size_t j=0;j<b;++j) {
            // Empty spans need no pointer arithmetic.
            moment[transposed?j*b+i:i*b+j]=n
                ? integer_dot32783(projection.data()+i*n,current.data()+j*n,n)
                : Z32783{};
        }
}
} // namespace VectorSpace::block_wiedemann_detail
