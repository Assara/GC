#pragma once
#include <cstddef>
#include "Field.hpp"

namespace VectorSpace {
// Differential entries are usually +/-1. Branch once per sparse entry,
// avoiding field construction and modular multiplication for every lane.
template<Field K>
inline void add_signed_block(K* out,const K* in,std::size_t width,SmallSignedInt coefficient) {
    if(coefficient==1) {
        for(std::size_t j=0;j<width;++j)out[j]+=in[j];
    } else if(coefficient==-1) {
        for(std::size_t j=0;j<width;++j)out[j]-=in[j];
    } else {
        const K scalar{int(coefficient)};
        for(std::size_t j=0;j<width;++j)out[j]+=in[j]*scalar;
    }
}
}
