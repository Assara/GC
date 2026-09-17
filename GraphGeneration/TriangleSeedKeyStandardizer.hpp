#pragma once

#include <algorithm>
#include <span>
#include <stdexcept>
#include <vector>
#include "graph.hpp"
#include "GraphStandardizer.hpp"
#include "GraphGeneration/ValenceGraphStandardizer.hpp"

namespace GraphGeneration {

// Baseline: canonicalize the full graph independently of construction history.
// The caller keeps working edges and stack masks in their original labels.
template <int MaxLoop>
class TriangleSeedKeyStandardizer {
    template <int V = 3, int E = (3 * (V - 1) + 1) / 2>
    static std::vector<Int> dispatch(int vertices, std::span<const Int> edges) {
        if (vertices == V) {
            if (edges.size() == 2 * E) {
                using G = Graph<V, E, 0, 0, 0, 0, fieldType>;
                G graph;
                std::copy(edges.begin(), edges.end(), graph.half_edges.begin());
                GraphStandardizer<V, E, 0, 0, 0, 0, fieldType> standardizer;
                const auto canonical = standardize_sorted_valences(graph, standardizer);
                const auto& data = canonical.half_edges;
                return {data.begin(), data.end()};
            }
            if constexpr (E < std::min(V * (V - 1) / 2, V - 1 + MaxLoop))
                return dispatch<V, E + 1>(vertices, edges);
        } else if constexpr (V < 2 * MaxLoop + 1) {
            return dispatch<V + 1>(vertices, edges);
        }
        throw std::invalid_argument("unreachable triangle graph dimensions");
    }
public:
    std::vector<Int> canonical_key(int vertices, std::span<const Int> edges) const {
        return dispatch(vertices, edges);
    }
};

} // namespace GraphGeneration
