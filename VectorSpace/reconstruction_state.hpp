#pragma once
#include <vector>
#include "packed_recurrence.hpp"

namespace VectorSpace::block_wiedemann_detail {
// A checkpoint is taken between complete Horner steps. Scratch is disposable.
template<Field K>
struct reconstruction_state {
    std::size_t rows, columns, degree, remaining;
    std::vector<K> accumulator;
    reconstruction_state(std::size_t rows,std::size_t columns,std::size_t degree)
        :rows(rows),columns(columns),degree(degree),remaining(degree),
         accumulator(recurrence_product(rows,columns)) {}
    void save(serialization::archive_output& out) const {
        out.section("reconstruction-state",[&] {
            out.word(rows);out.word(columns);out.word(degree);out.word(remaining);
            out.block<K>("accumulated-vectors",accumulator,rows,columns);
        });
    }
    void load(serialization::archive_input& in) {
        in.section("reconstruction-state",[&] {
            in.expect(rows);in.expect(columns);in.expect(degree);remaining=in.word();
            if(remaining>degree)throw std::runtime_error("invalid reconstruction progress");
            in.block<K>("accumulated-vectors",accumulator,rows,columns);
        });
    }
};
}
