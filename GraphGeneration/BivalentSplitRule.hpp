#pragma once

#include <algorithm>
#include <array>
#include <cassert>

#include "graph.hpp"

namespace GraphGeneration {

// Inputs are standardized: leaves precede bivalent vertices, then higher valences.
// Return the common neighbours of all bivalent vertices, or every vertex if none.
// A bivalent vertex cannot neighbour itself, so it is never selected.
template <typename G>
auto bivalent_split_vertices(const G& graph,
    const std::array<Int, G::N_VERTICES_>& valences) {
    assert(std::is_sorted(valences.begin(), valences.end()));
    std::array<bool, G::N_VERTICES_> eligible;
    eligible.fill(true);
    std::size_t vertex = 0;
    while (vertex < valences.size() && valences[vertex] < 2) ++vertex;
    for (; vertex < valences.size() && valences[vertex] == 2; ++vertex) {
        std::array<bool, G::N_VERTICES_> neighbours{};
        for (const auto position : graph.adjacent(static_cast<Int>(vertex))) {
            const auto other = position % 2 == 0 ? position + 1 : position - 1;
            neighbours[graph.half_edges[other]] = true;
        }
        for (std::size_t candidate = 0; candidate < eligible.size(); ++candidate)
            eligible[candidate] = eligible[candidate] && neighbours[candidate];
    }
    return eligible;
}

} // namespace GraphGeneration
