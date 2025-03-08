


#include <array>
#include <algorithm>
#include <cstring>
#include<vector>
#include <iostream>
#include <memory>
#include <utility>
#include "CombinatorialUtils.cpp"
#include "permutation.hpp"
#include "Types.hpp"


using namespace std;



template <
    Int N_VERTICES,
    Int N_EDGES,
    Int N_OUT_HAIR,
    Int N_IN_HAIR,
    signedInt c,
    signedInt d
>
class Graph
{
public:
    using SplitGraph = Graph<
        (N_VERTICES + 1),
        (N_EDGES + 1),
        N_OUT_HAIR,
        N_IN_HAIR,
        c,
        d
    >;

    static constexpr Int SIZE = N_IN_HAIR + N_OUT_HAIR + 2 * N_EDGES;
    static constexpr Int N_HAIR = N_IN_HAIR + N_OUT_HAIR;

    static constexpr signedInt FLIP_EDGE_SIGN = ((c % 2) != 0 && (d % 2) != 0) ? -1 : 1;
    static constexpr signedInt SWAP_EDGE_SIGN = (((c + d) % 2) != 0) ? -1 : 1;
    static constexpr signedInt SWAP_VERTICES_SIGN = -1*SWAP_EDGE_SIGN;


    Graph() = default;
//private:
    // Two arrays, each with fixed size. The first N_EDGES entries will be our edges.
    // outArray[i] is the "from" vertex for edge i
    // inArray[i]  is the "to"   vertex for edge i
    array<Int, SIZE> half_edges{};

    explicit Graph(const array<Int, SIZE>& arr)
        : half_edges(arr)
    {}

    pair<Int,Int> getEdge(Int i) {
        return {half_edges[ N_HAIR + 2*i],  N_HAIR + 2*i + 1 };
    }
    
    vector<Int> adjacent(Int v) const
    {
        vector<Int> result;
        result.reserve(SIZE);

        for (Int i = 0; i < SIZE; ++i) {
            if (half_edges[i] == v) {
                result.push_back(i);
            }
        }
        return result;
    }


    vector<unique_ptr<SplitGraph>> split_vertex_differential() const
    {
        vector<vector<Int>> adjRepresentation;
        adjRepresentation.reserve(N_VERTICES);

        bigInt resultSize = 0;

        for (Int v = 0; v < N_VERTICES; v++) {
            adjRepresentation.push_back(adjacent(v));
            resultSize += combutils::n_splits(adjRepresentation.back().size());
        }

        vector<unique_ptr<SplitGraph>>  result;
        result.reserve(resultSize);

        for (Int v = 0; v < N_VERTICES; v++) {
            split_vertex(v, adjRepresentation[v], result);
        }

        return result;
    }


    void split_vertex(Int split_vertex, vector<Int> &adjacent, vector<unique_ptr<SplitGraph>> &result) const
    {
        if (adjacent.size() < 4) {
            return;
        }
        Int max_index = adjacent.size()-1;
        for (Int i = 2; i < max_index; i++) {
            vector<Int> S = combutils::firstSubset(1, i);

            do {
                result.push_back(splitGraph(split_vertex, adjacent[0], S));
            } while (combutils::nextSubset(S, max_index));

        }
    }

    unique_ptr<SplitGraph> splitGraph(Int split_vertex, Int firstIndex, vector<Int> &S) const {
        auto sg = make_unique<SplitGraph>();

        // Copy half_edges from firstIndex onward.
        for (Int i = firstIndex; i < SIZE; ++i) {
            sg->half_edges[i] = half_edges[i];
        }

        // For each index in S, set that position in half_edges to N_VERTICES.
        for (auto s : S) {
            sg->half_edges[s] = N_VERTICES;
        }

        // Set the last two entries in half_edges:
        // - Second-to-last: split_vertex
        // - Last: N_VERTICES
        sg->half_edges[SplitGraph::SIZE - 2] = split_vertex;
        sg->half_edges[SplitGraph::SIZE - 1] = N_VERTICES;

        return sg;
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
        return FLIP_EDGE_SIGN;
    }

    signedInt directEdges() {
        signedInt sign = 1;
        for (Int i = 0; i < 2 * N_EDGES; i += 2) {
            // If the "from" vertex is greater than the "to" vertex,
            // swap them so that the lower value is in the "from" slot.
            if (half_edges[N_EDGES + i] > half_edges[N_EDGES + i +1]) {
                sign *= flipEdge(i);
            }
        }
        return sign;
    }


    signedInt compareEdge(Int e1, Int e2) const {
        Int base1 = N_EDGES + 2 * e1;
        Int base2 = N_EDGES + 2 * e2;
        return std::memcmp(&half_edges[base1], &half_edges[base2], 2 * sizeof(Int));
    }



    signedInt sortEdgesInsertion() {
        signedInt overallSign = 1;
        // Loop over edges starting at the second edge.
        for (Int i = 1; i < N_EDGES; ++i) {
            Int j = i;
            // While the current edge is less than the previous one, swap them.
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
        for (auto it = half_edges.begin(); it < half_edges.end(); ++it) {
            (*it) = perm[(*it)];
        }

        if (SWAP_VERTICES_SIGN == -1) {
            return perm.sign();
        }
        return 1;
    }

    // maybe return a unique_ptr instead. 
    signedInt swapVertices(Int v, Int w) const {
    // Iterate over all half_edges and swap every occurrence of i and j.
    for (Int k = 0; k < SIZE; ++k) {
        if (half_edges[k] == v) {
            half_edges[k] = w;
        } else if (half_edges[k] == w) {
            half_edges[k] = v;
        }
    }
        return SWAP_VERTICES_SIGN;
    }



};