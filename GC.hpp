#pragma once

#include "GraphStandardizer.hpp"
#include "VectorSpace/LinComb.hpp"
#include "VectorSpace/BasisElement.hpp"

template <
    Int N_VERTICES, Int N_EDGES,
    Int N_OUT_HAIR, Int N_IN_HAIR,
    signedInt c, signedInt d
>
class GC {
    using ThisGC = GC<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d>;
    using GraphType = Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d>;
    using SplitGC = GC<N_VERTICES + 1, N_EDGES + 1, N_OUT_HAIR, N_IN_HAIR, c, d>;

public:
    // Default constructor
    GC() = default;

    // Construct from a single BasisElement
    explicit GC(const BasisElement<GraphType>& b) : vec(b) {}

    // Construct from a GraphType (copying it into a BasisElement)
    explicit GC(const GraphType& G) : vec(BasisElement(G)) {}

    // Construct from a full LinComb
    explicit GC(const VectorSpace::LinComb<GraphType>& v) : vec(v) {}

    // Construct from a vector of unique_ptrs to GraphType and a coefficient
    explicit GC(std::vector<std::unique_ptr<GraphType>>&& graphs, fieldType coeff) {
        for (auto& ptr : graphs) {
            BasisElement<GraphType> be(*ptr, coeff);  // dereference to get GraphType
            vec.add(VectorSpace::LinComb<GraphType>(std::move(be)));
        }
        vec.standardize_and_sort();
    }

    // Access the underlying vector
    const VectorSpace::LinComb<GraphType>& data() const {
        return vec;
    }

    // += operator for summing GCs
    ThisGC& operator+=(const ThisGC& other) {
        vec += other.vec;
        return *this;
    }

    // Standardize and sort elements
    void standardize_all() { vec.standardize_all(); }
    void sort_elements()   { vec.sort_elements(); }
    void standardize_and_sort() { vec.standardize_and_sort(); }

    // Compute delta for the whole GC
    SplitGC delta() const {
        return delta_recursive(0, vec.raw_elements().size());
    }

    // Static method to compute delta of a single basis element
    static SplitGC delta(const BasisElement<GraphType>& G) {
        return SplitGC(G.getValue().split_vertex_differential(), G.getCoefficient());
    }

    void add(BasisElement<GraphType>&& elem) {
        vec.add(VectorSpace::LinComb<GraphType>(std::move(elem)));
    }
    
    void print() const {
        const auto& elems = vec.raw_elements();
        for (const auto& elem : elems) {
            elem.getValue().print();
            std::cout << "Coefficient: " << elem.getCoefficient() << "\n\n";
        }
    }

private:
    VectorSpace::LinComb<GraphType> vec;

    // Recursive delta helper
    SplitGC delta_recursive(size_t start, size_t end) const {
        if (end - start == 0) {
            return SplitGC(); // empty
        }
        if (end - start == 1) {
            const auto& be = vec.raw_elements()[start];
            return delta(be);
        }

        size_t mid = start + (end - start) / 2;
        SplitGC left = delta_recursive(start, mid);
        SplitGC right = delta_recursive(mid, end);
        left += right;
        return left;
    }
};