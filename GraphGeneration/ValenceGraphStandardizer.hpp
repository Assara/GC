#pragma once

#include <algorithm>
#include "permutation.hpp"

namespace GraphGeneration {

// Order an existing canonical representative by (degree, canonical index).
// This is deterministic on isomorphism classes, including ties in degree.
template <class Graph, class Standardizer>
Graph standardize_sorted_valences(const Graph& graph, const Standardizer& standardizer) {
    auto canonical = standardizer.standardize_no_sign(graph);
    const auto degrees = canonical.valence_array();
    Permutation<Graph::N_VERTICES_> order;
    std::sort(order.p.begin(), order.p.end(), [&](Int a, Int b) {
        return degrees[a] != degrees[b] ? degrees[a] < degrees[b] : a < b;
    });
    Graph result;
    result.assignPermutedDirectedSortedEdgesNoSign(canonical, order.inverse());
    return result;
}

} // namespace GraphGeneration
