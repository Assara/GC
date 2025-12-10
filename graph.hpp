#pragma once

#include <array>
#include <vector>
#include <memory>
#include <utility>
#include <cstring>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#include "types.hpp"
#include "permutation.hpp"
#include "CombinatorialUtils.hpp"
#include "VectorSpace/LinComb.hpp"
#include "graph_hash.hpp"
#include "GraphStandardizer.hpp"
#include <ranges>

using namespace std;

template<typename T, typename k>
class BasisElement;

template <
    Int N_VERTICES,
    Int N_EDGES,
    Int N_OUT_HAIR,
    Int N_IN_HAIR,
    signedInt c,
    signedInt d,
    typename fieldType
>
class Graph {
public:
    static constexpr Int SIZE = N_IN_HAIR + N_OUT_HAIR + (Int)2 * N_EDGES;
    static constexpr Int N_HAIR = N_IN_HAIR + N_OUT_HAIR;
    static constexpr Int N_EDGES_ = N_EDGES;
    static constexpr signedInt FLIP_EDGE_SIGN = ((c % 2) != 0 && (d % 2) != 0) ? -1 : 1;
    static constexpr signedInt SWAP_EDGE_SIGN = (((c + d) % 2) != 0) ? -1 : 1;
    static constexpr signedInt SWAP_VERTICES_SIGN = -1 * SWAP_EDGE_SIGN;

    using SplitGraph = Graph<N_VERTICES + 1, N_EDGES + 1, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>;
   
    using ContGraph = Graph<N_VERTICES - 1, N_EDGES - 1, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>;

    using ThisGraph = Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>;
    
    using ExtraEdgeGraph = Graph<N_VERTICES, N_EDGES + 1, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>;

    using Basis = BasisElement<ThisGraph, fieldType>;


    Graph() = default;

    explicit Graph(const array<Int, SIZE>& arr) : half_edges(arr) {}

    inline signedInt custom_filter() const {
        return std::ranges::count_if(valence_array(), [](Int v){ return v % 2 != 0; });
    }

    inline signedInt total_order_comp(const ThisGraph& other) const {
       return combutils::compareHalfEdges(half_edges, other.half_edges);
    }

    signedInt compare(const ThisGraph& other) const {
        
        // for homological filtration
        /*
        signedInt r = other.custom_filter() - custom_filter();
        if (r != 0) {
            return r;
        }
        */

        // For absolute ordering
        return total_order_comp(other);
    }

    bool operator<(const Graph& o) const { return compare(o) < 0; }

    static void std(Basis& b) {
			GraphStandardizer<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType> s;
			b = s.standardize(b);
	}

    static Basis canonized(Basis& b) {
		    GraphStandardizer<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType> s;
			return s.standardize(b);
	}
	
	ThisGraph canonical_represesentation() {
			Basis be = Basis(*this);
			return ThisGraph::canonized(be).getValue();
	}
	

    pair<Int, Int> getEdge(Int i) const {
        return { half_edges[N_HAIR + 2 * i], half_edges[N_HAIR + 2 * i + 1] };
    }

    void setEdge(Int i, Int v, Int w) {
        half_edges[N_HAIR + 2 * i] = v;
        half_edges[N_HAIR + 2 * i + 1] = w;
    }

    vector<Int> adjacent(Int v) const {
        vector<Int> result;
        result.reserve(SIZE);
        for (Int i = 0; i < SIZE; ++i) {
            if (half_edges[i] == v) {
                result.push_back(i);
            }
        }
        return result;
    }

    array<Int, N_VERTICES> valence_array() const {
        array<Int, N_VERTICES> valence{};
        for (Int v : half_edges) {
            ++valence[v];
        }
        return valence;
    }

    VectorSpace::LinComb<SplitGraph, fieldType> split_vertex_differential(fieldType coef) const {
            VectorSpace::LinComb<SplitGraph, fieldType> splits = unsorted_splits(coef);
            splits.standardize_and_sort();
            return splits;
    }



    VectorSpace::LinComb<SplitGraph, fieldType> unsorted_splits(fieldType coef) const {
            vector<vector<Int>> adjRepresentation;
            adjRepresentation.reserve(N_VERTICES);

            bigInt resultSize = 0;
            for (Int v = 0; v < N_VERTICES; v++) {
                    adjRepresentation.push_back(adjacent(v));
                    resultSize += combutils::n_splits(adjRepresentation.back().size());
            }

            VectorSpace::LinComb<SplitGraph, fieldType> result;
            result.reserve(resultSize);

            for (Int v = 0; v < N_VERTICES; v++) {
                    split_vertex(v, adjRepresentation[v], result, coef);
            }
            return result;
    }


    void split_vertex(Int split_vertex, const vector<Int>& adjacent, VectorSpace::LinComb<SplitGraph, fieldType>& result, fieldType coef) const {
        if (adjacent.size() < 4) {
    
            if (adjacent.size() == 2 ) {

                result.append_in_basis_order(splitGraph(split_vertex, adjacent, vector<Int>(adjacent.back())), -coef);
            }
            else if (adjacent.size() < 2) {
                result.append_in_basis_order(splitGraph(split_vertex, adjacent, vector<Int>()), coef);
            }
            //do nothing for adjacent.size() == 3
            return; 
        }
 
        Int max_index = adjacent.size() - 1;
        for (Int i = 2; i < max_index; i++) {
            vector<Int> S = combutils::firstSubset(1, i);

            do {
                result.append_in_basis_order(splitGraph(split_vertex, adjacent, S), coef);
            } while (combutils::nextSubset(S, max_index));
        }
    }

    SplitGraph splitGraph(Int split_vertex, const vector<Int>& adjacent, const vector<Int>& S) const {
            auto sg = SplitGraph();
            
            std::copy_n(half_edges.begin(), SIZE, sg.half_edges.begin());
            for (auto s : S) {
                sg.half_edges[adjacent[s]] = N_VERTICES;
            }

            sg.half_edges[SplitGraph::SIZE - 2] = split_vertex;
            sg.half_edges[SplitGraph::SIZE - 1] = N_VERTICES;

            return sg;
    }

    void add_splits_to_set(unordered_set<SplitGraph>& result) {
            VectorSpace::LinComb<SplitGraph, fieldType> splits = unsorted_splits();
            splits.standardize_all();

            
            constexpr fieldType zero{}; 
            for (auto const& be : splits) {
                if (be.getCoefficient() != zero) {
                    result.emplace(be.getValue());
                }
            }
    }


     VectorSpace::LinComb<SplitGraph, fieldType> unsorted_even_splits(fieldType coef) const {
            vector<vector<Int>> adjRepresentation;
            adjRepresentation.reserve(N_VERTICES);

            bigInt resultSize = 0;
            for (Int v = 0; v < N_VERTICES; v++) {
                    adjRepresentation.push_back(adjacent(v));
                    resultSize += combutils::n_splits(adjRepresentation.back().size());
            }

            VectorSpace::LinComb<SplitGraph, fieldType> result;
            result.reserve(resultSize);

            for (Int v = 0; v < N_VERTICES; v++) {
                    split_vertex_even(v, adjRepresentation[v], result, coef);
            }
            return result;
    }

    void split_vertex_even(Int vertex_to_split, const vector<Int>& adjacent, VectorSpace::LinComb<SplitGraph,fieldType>& result, fieldType coef) const {
        //odd vertices and bivalent vertices are handled the same way

        if (adjacent.size()%2 == 1 || adjacent.size() == 2) {
                split_vertex(vertex_to_split, adjacent, result, coef);
                return;
        }

        if (adjacent.size() < 5) {
            // do nothing for 4
            return; 
        }
 
        Int max_index = adjacent.size() - 1;
        for (Int i = 2; i < max_index; i++) {
            vector<Int> S = combutils::firstSubset(1, i);

            do {
                result.append_in_basis_order(splitGraph(vertex_to_split, adjacent, S), coef);
            } while (combutils::nextSubset(S, max_index));
        }
    }

    void add_even_splits_to_set(unordered_set<SplitGraph>& result) const {
            VectorSpace::LinComb<SplitGraph, fieldType> splits = unsorted_even_splits(fieldType{1});
            splits.standardize_all();

            constexpr fieldType zero{}; 
            for (auto const& be : splits) {
                if (be.getCoefficient() != zero) {
                    result.emplace(be.getValue());
                }
            }
    }

    VectorSpace::LinComb<ExtraEdgeGraph, fieldType> add_edge_differential(fieldType coef) const {
            VectorSpace::LinComb<ExtraEdgeGraph, fieldType> result;
            result.reserve(N_VERTICES * (N_VERTICES - 1) - N_EDGES);

                
            ExtraEdgeGraph base_graph;
            std::copy_n(half_edges.begin(), SIZE, base_graph.half_edges.begin() + 2);

            for (Int u = 0; u < N_VERTICES - 1; u++) {
                for (Int v = u+1; v < N_VERTICES; v++) {
                    base_graph.half_edges[0] = u;
                    base_graph.half_edges[1] = v;
                    result.append_in_basis_order(base_graph, coef);
                }
            }
            result.standardize_and_sort();

            return result;
    }


    //only suitable after sorting and directing edges
    bool has_double_edge() {
        for (Int i = 0; i < ThisGraph::N_EDGES_ - 1; i++) {
            if (getEdge(i) == getEdge(i+1)) {
                return true;
            }
        }
        return false;

    }
    
    VectorSpace::LinComb<ContGraph, fieldType> contraction_differential(fieldType k) const {
		
		VectorSpace::LinComb<ContGraph, fieldType> result;
		
		for(Int i = 0; i < N_EDGES_; i++) {
				result.append_in_basis_order(contract_edge(i, k));
		}
		result.standardize_and_sort();
		return result;
	}
	
	unordered_set<ContGraph> contraction_set() const {
		
		unordered_set<ContGraph> result;
		result.reserve(N_EDGES_);
		
		for(Int i = 0; i < N_EDGES_; i++) {
				auto be = contract_edge(i, fieldType{1});
				ContGraph::std(be);
				
				if (be.getCoefficient() != fieldType{0}) {
						result.emplace(be.getValue());
				}		 
		}
		return result;
	}
	
    
    
    BasisElement<ContGraph, fieldType> contract_edge(Int i, fieldType k) const {
            const Int edge_index = N_HAIR + 2 * i;
            const Int contraction_vertex = min(half_edges[edge_index], half_edges[edge_index + 1]);
            const Int deletion_vertex = max(half_edges[edge_index], half_edges[edge_index + 1]);
     
            // skip for tadpoles
            if (contraction_vertex == deletion_vertex) {
                    return BasisElement<ContGraph, fieldType>(ContGraph{}, static_cast<fieldType>(0));
            }

            BasisElement<ContGraph, fieldType> contracted(ContGraph(), k);
 
            for (Int j = 0; j < edge_index; ++j) {
                    contracted.getValue().half_edges[j] =  
                            contraction_value(half_edges[j], contraction_vertex, deletion_vertex);
            } 

            for (Int j = edge_index + 2 ; j < ThisGraph::SIZE; ++j) {
                    contracted.getValue().half_edges[j - 2] =  // shift down by 2 to account for removed edge
                            contraction_value(half_edges[j], contraction_vertex, deletion_vertex);
            }
            
            if ( (N_EDGES - i) % 2 == 0) {
                    contracted.multiplyCoefficient(SWAP_EDGE_SIGN);
            } 

            if (deletion_vertex != N_VERTICES - 1) {
                    contracted.multiplyCoefficient(SWAP_VERTICES_SIGN);
            } 

            if (half_edges[edge_index] > half_edges[edge_index + 1]) {
                    contracted.multiplyCoefficient(FLIP_EDGE_SIGN);
            }

            return contracted;
	}
    

    static BasisElement<ContGraph, fieldType> contract_edge(const BasisElement<ThisGraph, fieldType>& be, Int i) {
            return be.getValue().contract_edge(i, be.getCoefficient());
    }

    static Int contraction_value(Int v, Int contraction_vertex, Int deletion_vertex) {
        if (v == deletion_vertex) {
            return contraction_vertex;
        }
        if (v == N_VERTICES - 1) {
            return deletion_vertex;
        }
        return v;

    }

    signedInt flipEdge(Int i) {
        const Int base = N_HAIR + 2 * i;
        std::swap(half_edges[base], half_edges[base + 1]);
        return FLIP_EDGE_SIGN;
    }

    signedInt swapEdges(Int i, Int j) {
        const Int base_i = N_HAIR + 2 * i;
        const Int base_j = N_HAIR + 2 * j;
        std::swap(half_edges[base_i], half_edges[base_j]);
        std::swap(half_edges[base_i + 1], half_edges[base_j + 1]);
        return SWAP_EDGE_SIGN;
    }

    signedInt directEdges() {
        signedInt sign = 1;
        for (Int edgeIndex = 0; edgeIndex < N_EDGES; ++edgeIndex) {
            Int base = N_HAIR + 2 * edgeIndex;
            if (half_edges[base] > half_edges[base + 1]) {
                sign *= flipEdge(edgeIndex);
            }
        }
        return sign;
    }

    signedInt compareEdge(Int e1, Int e2) const {
        Int base1 = N_HAIR + 2 * e1;
        Int base2 = N_HAIR + 2 * e2;
        return std::memcmp(&half_edges[base1], &half_edges[base2], 2 * sizeof(Int));
    }

    signedInt sortEdgesInsertion() {
        signedInt overallSign = 1;
        for (Int i = 1; i < N_EDGES; ++i) {
            Int j = i;
            while (j > 0 && compareEdge(j - 1, j) > 0) {
                overallSign *= swapEdges(j - 1, j);
                j--;
            }
        }
        return overallSign;
    }

    signedInt directAndSortEdges() {
        return directEdges() * sortEdgesInsertion();
    }

    signedInt permuteVertices(Permutation<N_VERTICES> perm) {
        for (auto& h : half_edges) {
            h = perm[h];
        }
        return (SWAP_VERTICES_SIGN == -1) ? perm.sign() : 1;
    }

    signedInt swapVertices(Int v, Int w) {
        if (v == w) return 1;
        for (Int& h : half_edges) {
            if (h == v) h = w;
            else if (h == w) h = v;
        }
        return SWAP_VERTICES_SIGN;
    }

	void print(std::ostream& out = std::cout) const {
		out << "printing graph:\n";

		if (N_OUT_HAIR > 0) {
			out << "out_hair: ";
			for (Int i = 0; i < N_OUT_HAIR; ++i) {
				out << +half_edges[i];
				if (i < N_OUT_HAIR - 1) out << ", ";
			}
			out << "\n";
		}

		if (N_IN_HAIR > 0) {
			out << "in_hair: ";
			for (Int i = 0; i < N_IN_HAIR; ++i) {
				out << +half_edges[N_OUT_HAIR + i];
				if (i < N_IN_HAIR - 1) out << ", ";
			}
			out << "\n";
		}

		out << "edges: ";
		for (Int e = 0; e < N_EDGES; ++e) {
			Int base = N_HAIR + 2 * e;
			out << "(" << +half_edges[base] << "," << +half_edges[base + 1] << ")";
			if (e < N_EDGES - 1) out << ", ";
		}
		out << "\n";

		out << "grade = " << custom_filter() << "\n";
	}
	
	
    bool operator==(Graph const& other) const {
        return half_edges == other.half_edges;
    }

    bool operator!=(Graph const& other) const {
        return !(*this == other);
    }

    std::array<Int, SIZE> half_edges{};


    Int half_edge(Int i) const {
        return half_edges[i];
    }

};
