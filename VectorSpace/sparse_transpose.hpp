#pragma once
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <numeric>
#include <span>
#include <vector>
#include "sparse_block.hpp"
#include "OwnedArray.hpp"

namespace VectorSpace {
// CSC storage for the transpose, retaining one byte per differential coefficient.
struct signed_sparse_transpose {
    std::size_t rows;
    OwnedArray<std::size_t> offsets;
    OwnedArray<std::uint32_t> row_indices;
    OwnedArray<SmallSignedInt> coefficients;
    template<class Matrix>
    explicit signed_sparse_transpose(const Matrix& matrix)
        :rows(matrix.columns()),offsets(std::size_t(matrix.rows)+1),
         row_indices(matrix.coefficients.size()),coefficients(matrix.coefficients.size()) {
        for(auto row:matrix.row_indices)++offsets[std::size_t(row)+1];
        std::partial_sum(offsets.begin(),offsets.end(),offsets.begin());
        OwnedArray<std::size_t> next{std::span<const std::size_t>(offsets)};
        for(std::size_t c=0;c<matrix.columns();++c)
            for(auto i=matrix.offsets[c];i<matrix.offsets[c+1];++i) {
                const auto target=next[matrix.row_indices[i]]++;
                row_indices[target]=c;
                coefficients[target]=matrix.coefficients[i];
            }
    }
    std::size_t columns() const {return offsets.size()-1;}
    std::size_t allocated_bytes() const {
        return offsets.size()*sizeof(std::size_t)
            +row_indices.size()*sizeof(std::uint32_t)+coefficients.size();
    }
};

// One column of CSC is one independent output row of its transpose.
// Caller owns the output row; no atomics or reduction buffers are needed.
template<Field K,class Matrix>
inline void gather_signed_column(const Matrix& matrix,std::size_t c,
    const K* input,K* output,std::size_t b) {
    for(auto i=matrix.offsets[c];i<matrix.offsets[c+1];++i)
        add_signed_block(output,input+matrix.row_indices[i]*b,b,matrix.coefficients[i]);
}

template<Field K,class Matrix>
void evaluate_signed_transpose(const Matrix& matrix,std::span<const K> in,
    std::span<K> out,std::size_t b,int threads) {
    // Small products stay serial to avoid OpenMP launch overhead.
    const bool parallel=matrix.coefficients.size()*b>=32768;
    (void)parallel; (void)threads;
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
    for(std::size_t c=0;c<matrix.columns();++c) {
        auto* row=out.data()+c*b;
        std::fill_n(row,b,K{});
        gather_signed_column<K>(matrix,c,in.data(),row,b);
    }
}
}
