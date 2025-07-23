#pragma once

#include <array>
#include <vector>
#include <memory>
#include <utility>
#include <cstring>
#include <iostream>
#include "Types.hpp"
#include "permutation.hpp"
#include "CombinatorialUtils.hpp"
#include "VectorSpace/HasCompare.hpp"

using namespace std;

template<typename T>
class BasisElement;

template <
    Int N_VERTICES,
    Int N_EDGES,
    Int N_OUT_HAIR,
    Int N_IN_HAIR,
    signedInt c,
    signedInt d
>
class Graph {
public:
    static constexpr Int SIZE = N_IN_HAIR + N_OUT_HAIR + 2 * N_EDGES;
    static constexpr Int N_HAIR = N_IN_HAIR + N_OUT_HAIR;
    static constexpr Int EDGES = N_EDGES;
    static constexpr signedInt FLIP_EDGE_SIGN = ((c % 2) != 0 && (d % 2) != 0) ? -1 : 1;
    static constexpr signedInt SWAP_EDGE_SIGN = (((c + d) % 2) != 0) ? -1 : 1;
    static constexpr signedInt SWAP_VERTICES_SIGN = -1 * SWAP_EDGE_SIGN;

    using SplitGraph = std::conditional_t<
        N_VERTICES + 1 < SIZE,
        Graph<N_VERTICES + 1, N_EDGES + 1, N_OUT_HAIR, N_IN_HAIR, c, d>,
        void
    >;

    using ThisGraph = Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d>;

    Graph() = default;

    explicit Graph(const array<Int, SIZE>& arr) : half_edges(arr) {}

    signedInt compare(const ThisGraph& other) const {
        return combutils::compareHalfEdges(half_edges, other.half_edges);
    }

    static void std(BasisElement<ThisGraph>& b);

    pair<Int, Int> getEdge(Int i) {
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

    vector<unique_ptr<SplitGraph>> split_vertex_differential() const {
        if constexpr (std::is_same_v<SplitGraph, void>) {
            static_assert(!std::is_same_v<SplitGraph, void>,
                "split_vertex_differential() called on a Graph without valid SplitGraph");
        } else {
            vector<vector<Int>> adjRepresentation;
            adjRepresentation.reserve(N_VERTICES);

            bigInt resultSize = 0;
            for (Int v = 0; v < N_VERTICES; v++) {
                adjRepresentation.push_back(adjacent(v));
                resultSize += combutils::n_splits(adjRepresentation.back().size());
            }

            vector<unique_ptr<SplitGraph>> result;
            result.reserve(resultSize);

            for (Int v = 0; v < N_VERTICES; v++) {
                split_vertex(v, adjRepresentation[v], result);
            }

            return result;
        }
    }

    void split_vertex(Int split_vertex, vector<Int>& adjacent, vector<unique_ptr<SplitGraph>>& result) const {
        if constexpr (std::is_same_v<SplitGraph, void>) return;
        if (adjacent.size() < 4) return;

        Int max_index = adjacent.size() - 1;
        for (Int i = 2; i < max_index; i++) {
            vector<Int> S = combutils::firstSubset(1, i);
            do {
                result.push_back(splitGraph(split_vertex, adjacent[0], S));
            } while (combutils::nextSubset(S, max_index));
        }
    }

    auto splitGraph(Int split_vertex, Int firstIndex, vector<Int>& S) const {
        if constexpr (std::is_same_v<SplitGraph, void>) {
            static_assert(!std::is_same_v<SplitGraph, void>,
                "splitGraph() called on a Graph without valid SplitGraph");
        } else {
            auto sg = make_unique<SplitGraph>();
            for (Int i = firstIndex; i < SIZE; ++i) {
                sg->half_edges[i] = half_edges[i];
            }

            for (auto s : S) {
                sg->half_edges[s] = N_VERTICES;
            }

            sg->half_edges[SplitGraph::SIZE - 2] = split_vertex;
            sg->half_edges[SplitGraph::SIZE - 1] = N_VERTICES;

            return sg;
        }
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

    void print() const {
        if (N_OUT_HAIR > 0) {
            cout << "out_hair: ";
            for (Int i = 0; i < N_OUT_HAIR; ++i) {
                cout << half_edges[i];
                if (i < N_OUT_HAIR - 1) cout << ", ";
            }
            cout << "\n";
        }

        if (N_IN_HAIR > 0) {
            cout << "in_hair: ";
            for (Int i = 0; i < N_IN_HAIR; ++i) {
                cout << half_edges[N_OUT_HAIR + i];
                if (i < N_IN_HAIR - 1) cout << ", ";
            }
            cout << "\n";
        }

        cout << "edges: ";
        for (Int e = 0; e < N_EDGES; ++e) {
            Int base = N_HAIR + 2 * e;
            cout << "(" << half_edges[base] << "," << half_edges[base + 1] << ")";
            if (e < N_EDGES - 1) cout << ", ";
        }
        cout << endl;
    }

    bool operator==(Graph const& other) const {
        return half_edges == other.half_edges;
    }

    bool operator!=(Graph const& other) const {
        return !(*this == other);
    }

    std::array<Int, SIZE> half_edges{};
};
