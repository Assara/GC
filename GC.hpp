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
    using ContGC = GC<N_VERTICES - 1, N_EDGES - 1, N_OUT_HAIR, N_IN_HAIR, c, d>;
private:
        VectorSpace::LinComb<GraphType, fieldType> vec;

public:
    // Default constructor
    GC() = default;

    // Construct from a single BasisElement
    explicit GC(const BasisElement<GraphType, fieldType>& b) : vec(b) {}

    // Construct from a GraphType (copying it into a BasisElement)
    explicit GC(const GraphType& G) : vec(G) {}


    // Construct from a ptr to GraphType (moving it into a BasisElement)
    explicit GC(std::unique_ptr<GraphType>&& g, fieldType coeff = 1.0f)
    : vec(BasisElement<GraphType, fieldType>(std::move(g), coeff)) {}

    // Construct from a full LinComb
    explicit GC(const VectorSpace::LinComb<GraphType, fieldType>& v) : vec(v) {}

    // Construct from a vector of unique_ptrs to GraphType and a coefficient
    explicit GC(std::vector<std::unique_ptr<GraphType>>&& graphs, fieldType coeff)
                : vec(std::move(graphs), coeff) {} 

    explicit GC(std::vector<BasisElement<GraphType, fieldType>>&& elems)
                : vec(std::move(elems)) {} 

    // Access the underlying vector
    const VectorSpace::LinComb<GraphType, fieldType>& data() const {
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

    BasisElement<GraphType, fieldType>& back() {
            return vec.back();
    }

    // Compute delta for the whole GC
    SplitGC delta() const {
            return delta_recursive(0, vec.raw_elements().size());
    }

    // Static method to compute delta of a single basis element
    static SplitGC delta(const BasisElement<GraphType, fieldType>& G) {
            return SplitGC(G.getValue().split_vertex_differential(), G.getCoefficient());
    }

    vector<BasisElement<typename GraphType::ContGraph, fieldType>> d_contraction_without_sort() const {
            using ContGraph = typename GraphType::ContGraph;
            vector<BasisElement<typename GraphType::ContGraph, fieldType>> result;
            result.reserve(vec.raw_elements().size() * GraphType::N_EDGES_);

            for (const auto& elem : vec.raw_elements()) {
                    const auto& be = elem;
                    for (Int i = 0; i < GraphType::N_EDGES_; ++i) {
                            // Contract the i-th edge
                            BasisElement<ContGraph, fieldType> contracted = GraphType::contract_edge(be, i);
                            
                            if (contracted.getCoefficient() != 0) {
                                result.push_back(std::move(contracted));;
                            }
                    }
            }
            return result;
    }

    ContGC d_contraction() {
            std::vector<BasisElement<typename GraphType::ContGraph, fieldType>> elems = d_contraction_without_sort();
            ContGC dThis(std::move(elems));
            return dThis;
    }

    void add(BasisElement<GraphType, fieldType>&& elem) {
            vec.add(VectorSpace::LinComb<GraphType, fieldType>(std::move(elem)));
    }
    
    void print() const {
            const auto& elems = vec.raw_elements();
            for (const auto& elem : elems) {
                    elem.getValue().print();
                    std::cout << "Coefficient: " << elem.getCoefficient() << "\n\n";
            }
    }


    //This should not be here
    void candidates() {
            if (vec.raw_elements().empty()) return;

            const GraphType& G = vec.back().getValue();

            std::vector<std::unique_ptr<typename GraphType::SplitGraph>> splits;
            G.split_vertex(0, G.adjacent(0), splits);

            vector<ThisGC> boundaries;
            for (auto& split : splits) {
                    boundaries.push_back(ThisGC::SplitGC(std::move(split)).d_contraction());
            }


            // scale boundaries 

        /*    for (ThisGC bouundary : boundaries) {
                    if (boundary.back().getValueRef().compare(G) == 0) {

                        
                    }

            }
        */

          //  cout << "aaa = " << boundaries.size() << endl;

    }

    GC& scalar_multiply(fieldType scalar) {
            vec.scalar_multiply(scalar);
            return *this;
    }


private:

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