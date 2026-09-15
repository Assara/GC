#pragma once

#include <array>
#include <bit>
#include <cstdint>

namespace GraphGeneration {

// Cut-free, simple parents only. Test the proposed incidence partition before
// constructing a child. Prefer another triangle-free edge with larger endpoint
// degree sum whose contraction is cut-free; equal scores remain eligible.
template <class Graph>
class SplitReductionRule {
    static constexpr int V = Graph::N_VERTICES_;
    static_assert(V + 1 < 64);
    using Mask = std::uint64_t;
    using Adjacency = std::array<Mask, V + 1>;
    const Graph& graph_;
    Adjacency adjacent_{};
    std::array<int, V> degrees_{};
    Mask skipped_ = 0;

    static bool connected_without(const Adjacency& adjacent, Mask all, int a, int b) {
        auto remaining = all & ~(Mask{1} << a) & ~(Mask{1} << b);
        auto pending = remaining & -remaining;
        remaining &= ~pending;
        while (pending) {
            const int v = std::countr_zero(pending);
            pending &= pending - 1;
            const auto added = adjacent[v] & remaining;
            remaining &= ~added;
            pending |= added;
        }
        return !remaining;
    }

public:
    explicit SplitReductionRule(const Graph& graph) : graph_(graph) {
        for (int e = 0; e < Graph::N_EDGES_; ++e) {
            const auto [a, b] = graph.getEdge(e);
            adjacent_[a] |= Mask{1} << b;
            adjacent_[b] |= Mask{1} << a;
        }
        for (int v = 0; v < V; ++v) degrees_[v] = std::popcount(adjacent_[v]);
        // A split's new-edge score is always degree(v)+2. For an eligible
        // parent edge disjoint from v, contracting and splitting commute.
        // Such an edge therefore rejects every partition of v.
        for (int e = 0; e < Graph::N_EDGES_; ++e) {
            const auto [a, b] = graph.getEdge(e);
            if (adjacent_[a] & adjacent_[b]) continue;
            Mask victims = 0;
            for (int v = 0; v < V; ++v)
                if (v != a && v != b && degrees_[v] > 3
                    && degrees_[a] + degrees_[b] > degrees_[v] + 2)
                    victims |= Mask{1} << v;
            if ((victims & ~skipped_)
                && connected_without(adjacent_, (Mask{1} << V) - 1, a, b))
                skipped_ |= victims;
        }
    }

    bool skip_vertex(int vertex) const { return skipped_ & (Mask{1} << vertex); }

    bool redundant(int vertex, Mask moved) const {
        // Update only the moved incidences in a temporary adjacency view.
        // No child Graph, edge copying or child valence scan is needed.
        auto adjacent = adjacent_;
        adjacent[vertex] = (adjacent[vertex] & ~moved) | (Mask{1} << V);
        adjacent[V] = moved | (Mask{1} << vertex);
        for (auto pending = moved; pending; pending &= pending - 1) {
            const int v = std::countr_zero(pending);
            adjacent[v] = (adjacent[v] & ~(Mask{1} << vertex)) | (Mask{1} << V);
        }
        const int count = std::popcount(moved);
        const int score = degrees_[vertex] + 2;
        const auto degree = [&](int v) {
            return v == V ? count + 1
                : v == vertex ? degrees_[vertex] - count + 1 : degrees_[v];
        };
        for (int e = 0; e < Graph::N_EDGES_; ++e) {
            auto [a, b] = graph_.getEdge(e);
            if (a == vertex && (moved & (Mask{1} << b))) a = V;
            else if (b == vertex && (moved & (Mask{1} << a))) b = V;
            if (degree(a) + degree(b) <= score || (adjacent[a] & adjacent[b])) continue;
            if (connected_without(adjacent, (Mask{1} << (V + 1)) - 1, a, b)) return true;
        }
        return false;
    }
};

} // namespace GraphGeneration
