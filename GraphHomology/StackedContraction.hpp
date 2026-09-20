#pragma once
#include "ContractionMatrices.hpp"
#include "VectorSpace/sparse_transpose.hpp"

namespace GraphHomology {
// Unweighted [down; up^T], borrowing originals and owning their transposes.
// Every parallel iteration owns its entire output row, including both terms
// of the transposed stack. No output scatter or atomics.
struct StackedContraction {
    const ContractionMatrix& down;
    const ContractionMatrix& up;
    VectorSpace::signed_sparse_transpose down_transpose,up_transpose;
    int threads;
    StackedContraction(const ContractionMatrix& a,const ContractionMatrix& b,int workers=8)
        :down(a),up(b),down_transpose(a),up_transpose(b),threads(workers) {
        if(a.columns()!=b.rows)throw std::invalid_argument("incompatible stacked differentials");
        if(workers<1)throw std::invalid_argument("thread count must be positive");
    }
    std::size_t rows() const {return std::size_t(down.rows)+up.columns();}
    std::size_t columns() const {return down.columns();}
    std::size_t transpose_allocated_bytes() const {
        return down_transpose.allocated_bytes()+up_transpose.allocated_bytes();
    }
    template<VectorSpace::Field K>
    void apply(std::span<const K> in,std::span<K> out,std::size_t b) const {
        const bool parallel=(down.nonzeros()+up.nonzeros())*b>=32768;
        (void)parallel;
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
        for(std::size_t r=0;r<rows();++r) {
            auto* row=out.data()+r*b;
            std::fill_n(row,b,K{});
            if(r<down.rows)
                VectorSpace::gather_signed_column<K>(down_transpose,r,in.data(),row,b);
            else
                VectorSpace::gather_signed_column<K>(up,r-down.rows,in.data(),row,b);
        }
    }
    template<VectorSpace::Field K>
    void transpose(std::span<const K> in,std::span<K> out,std::size_t b) const {
        const bool parallel=(down.nonzeros()+up.nonzeros())*b>=32768;
        (void)parallel;
#pragma omp parallel for schedule(static) num_threads(threads) if(parallel)
        for(std::size_t c=0;c<columns();++c) {
            auto* row=out.data()+c*b;
            std::fill_n(row,b,K{});
            VectorSpace::gather_signed_column<K>(down,c,in.data(),row,b);
            // Avoid pointer arithmetic on an empty span in a zero-row stack.
            if(up.columns())VectorSpace::gather_signed_column<K>(up_transpose,c,
                in.data()+std::size_t(down.rows)*b,row,b);
        }
    }
};
}
