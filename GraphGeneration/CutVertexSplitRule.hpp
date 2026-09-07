#pragma once

#include <array>
#include <bit>
#include <cstdint>

namespace GraphGeneration {

// Connected, simple, hairless graphs only. Select one cut vertex and record
// its neighbours in each component of G-v. No canonicalization or splitting.
template <int V>
struct CutVertexSplitRule {
    static_assert(V >= 3 && V < 64);
    int vertex = -1;
    int branch_count = 0;
    std::array<std::uint64_t, V> branches{};

    template <class Graph>
    explicit CutVertexSplitRule(const Graph& graph) {
        std::array<std::uint64_t, V> adjacent{};
        for (int e = 0; e < graph.N_EDGES_; ++e) {
            const auto [a, b] = graph.getEdge(e);
            adjacent[a] |= std::uint64_t{1} << b;
            adjacent[b] |= std::uint64_t{1} << a;
        }
        constexpr auto all = (std::uint64_t{1} << V) - 1;
        for (int removed = 0; removed < V; ++removed) {
            auto remaining = all & ~(std::uint64_t{1} << removed);
            branch_count = 0;
            while (remaining) {
                auto pending = remaining & -remaining;
                std::uint64_t component = 0;
                remaining &= ~pending;
                while (pending) {
                    const int v = std::countr_zero(pending);
                    pending &= pending - 1;
                    component |= std::uint64_t{1} << v;
                    const auto added = adjacent[v] & remaining;
                    remaining &= ~added;
                    pending |= added;
                }
                branches[branch_count++] = component & adjacent[removed];
            }
            if (branch_count > 1) {
                vertex = removed;
                return;
            }
        }
        branch_count = 0;
    }

    // Both replacement vertices must meet every old branch. Otherwise one
    // replacement remains a cut vertex. Complementary partitions are equivalent.
    bool resolves(std::uint64_t moved_neighbours) const noexcept {
        for (int i = 0; i < branch_count; ++i)
            if (!(branches[i] & moved_neighbours)
                || !(branches[i] & ~moved_neighbours)) return false;
        return true;
    }
};

} // namespace GraphGeneration
