#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include "../Z32783.hpp"
#include "packed_projection.hpp"

namespace VectorSpace::experimental {
// Restricted to canonical residues of this prime. Not a generic Field kernel.
inline constexpr std::uint64_t dot_prime=32783;
inline constexpr std::size_t dot_chunk=1u<<20;
static_assert(dot_chunk*(dot_prime-1)*(dot_prime-1)<(1ULL<<51));
enum class DotMethod { field, integer, floating_integer_reduce, floating_reciprocal };
inline double dot_value(double x) {return x;}
inline std::uint32_t dot_value(Z32783 x) {return x.value();}

// x is an exact nonnegative integer below 2^51. Quotient error is at most one;
// q*p and the subtraction are exact integers in double at this bound.
inline std::uint32_t reduce_double32783(double x) {
    const double q=std::floor(x*(1.0/double(dot_prime)));
    double r=x-q*double(dot_prime);
    if(r<0)r+=double(dot_prime);
    if(r>=double(dot_prime))r-=double(dot_prime);
    return static_cast<std::uint32_t>(r);
}

template<DotMethod method,class Input>
Z32783 dot32783(const Input* a,const Input* b,std::size_t n,
    std::size_t stride_a=1,std::size_t stride_b=1) {
    if constexpr(method==DotMethod::integer && std::same_as<Input,Z32783>) {
        if(stride_a==1 && stride_b==1)
            return block_wiedemann_detail::integer_dot32783(a,b,n);
    }
    Z32783 result{};
    for(std::size_t first=0;first<n;first+=dot_chunk) {
        const auto count=std::min(dot_chunk,n-first);
        const auto* x=a+first*stride_a;
        const auto* y=b+first*stride_b;
        if constexpr(method==DotMethod::field) {
            for(std::size_t i=0;i<count;++i)
                result+=x[i*stride_a]*y[i*stride_b];
        } else if constexpr(method==DotMethod::integer) {
            std::uint64_t sum=0;
#pragma omp simd reduction(+:sum)
            for(std::size_t i=0;i<count;++i)
                sum+=std::uint64_t(dot_value(x[i*stride_a]))*std::uint64_t(dot_value(y[i*stride_b]));
            result+=Z32783(sum%dot_prime);
        } else {
            double sum=0;
#pragma omp simd reduction(+:sum)
            for(std::size_t i=0;i<count;++i)
                sum+=double(dot_value(x[i*stride_a]))*double(dot_value(y[i*stride_b]));
            if constexpr(method==DotMethod::floating_integer_reduce)
                result+=Z32783(static_cast<std::uint64_t>(sum)%dot_prime);
            else result+=Z32783(reduce_double32783(sum));
        }
    }
    return result;
}
} // namespace VectorSpace::experimental
