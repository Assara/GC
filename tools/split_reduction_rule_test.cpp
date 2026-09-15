#include <cassert>
#include <iostream>
#include "graph.hpp"
#include "GraphGeneration/CutVertexSplitRule.hpp"
#include "GraphGeneration/SplitReductionRule.hpp"

std::size_t checked = 0, rejected = 0, separating = 0;

template <class Child>
bool check_child(const Child& child, int vertex) {
    using Contracted = Graph<Child::N_VERTICES_ - 1, Child::N_EDGES_ - 1,
                             0, 0, 0, 0, fieldType>;
    const auto child_degrees = child.valence_array();
    bool expected = false;
    // Independent oracle: actually contract every better edge,
    // reject parallel edges, then run the existing cut finder.
    for (int e = 0; e < Child::N_EDGES_; ++e) {
        const auto [a, b] = child.getEdge(e);
        if (child_degrees[a] + child_degrees[b]
            <= child_degrees[vertex] + child_degrees[Child::N_VERTICES_ - 1]) continue;
        Contracted contracted;
        std::array<std::array<bool, Child::N_VERTICES_ - 1>, Child::N_VERTICES_ - 1> seen{};
        bool simple = true;
        int output_edge = 0;
        const auto map = [a, b](int v) {
            if (v == b) v = a;
            return v - (v > b);
        };
        for (int f = 0; f < Child::N_EDGES_; ++f) {
            if (f == e) continue;
            const auto [x, y] = child.getEdge(f);
            const int u = map(x), v = map(y);
            if (u == v || seen[u][v]) { simple = false; break; }
            seen[u][v] = seen[v][u] = true;
            contracted.setEdge(output_edge++, u, v);
        }
        if (!simple) continue;
        if (GraphGeneration::CutVertexSplitRule<Child::N_VERTICES_ - 1>(contracted).vertex >= 0) {
            ++separating;
            continue;
        }
        expected = true;
    }
    ++checked;
    rejected += expected;
    return expected;
}

template <int E>
void check() {
    using Parent = Graph<6, E, 0, 0, 0, 0, fieldType>;
    for (unsigned mask = 0; mask < (1u << 15); ++mask) {
        if (std::popcount(mask) != E) continue;
        Parent parent;
        int bit = 0, edge = 0;
        for (int b = 1; b < 6; ++b)
            for (int a = 0; a < b; ++a, ++bit)
                if (mask & (1u << bit)) parent.setEdge(edge++, a, b);
        const auto degrees = parent.valence_array();
        if (std::ranges::any_of(degrees, [](Int d) { return d < 3; })
            || GraphGeneration::CutVertexSplitRule<6>(parent).vertex >= 0) continue;
        const GraphGeneration::SplitReductionRule reduction(parent);
        for (int vertex = 0; vertex < 6; ++vertex) {
            const auto adjacent = parent.adjacent(vertex);
            const int last = adjacent.size() - 1;
            for (int moved = 2; moved < last; ++moved) {
                auto subset = combutils::firstSubset(1, moved);
                do {
                    const auto child = parent.splitGraph(vertex, adjacent, subset);
                    std::uint64_t moved_neighbours = 0;
                    for (auto index : subset)
                        moved_neighbours |= std::uint64_t{1} << parent.half_edges[adjacent[index] ^ 1];
                    const bool expected = check_child(child, vertex);
                    assert(reduction.redundant(vertex, moved_neighbours) == expected);
                    assert(!reduction.skip_vertex(vertex) || expected);
                } while (combutils::nextSubset(subset, last));
            }
        }
    }
    if constexpr (E < 15) check<E + 1>();
}

int main() {
    check<9>();
    // Two K4s joined through adjacent separating vertices 8 and 9.
    Graph<10, 21, 0, 0, 0, 0, fieldType> parent;
    int e = 0;
    parent.setEdge(e++, 8, 9);
    for (int base : {0, 4}) {
        for (int b = base + 1; b < base + 4; ++b)
            for (int a = base; a < b; ++a) parent.setEdge(e++, a, b);
        for (int v = base; v < base + 4; ++v)
            parent.setEdge(e++, v, v < base + 2 ? 8 : 9);
    }
    assert(GraphGeneration::CutVertexSplitRule<10>(parent).vertex < 0);
    const auto child = parent.splitGraph(0, parent.adjacent(0), vector<Int>{1, 2});
    const bool expected = check_child(child, 0);
    // The only better triangle-free edge is 8--9, whose contraction has a cut.
    const GraphGeneration::SplitReductionRule reduction(parent);
    assert(!expected && !reduction.skip_vertex(0));
    assert(!reduction.redundant(0, (std::uint64_t{1} << 2) | (std::uint64_t{1} << 3)));
    // K4,4 has a better eligible edge disjoint from every degree-four
    // vertex, so every vertex can be skipped without enumerating partitions.
    Graph<8, 16, 0, 0, 0, 0, fieldType> bipartite;
    e = 0;
    for (int a = 0; a < 4; ++a)
        for (int b = 4; b < 8; ++b) bipartite.setEdge(e++, a, b);
    const GraphGeneration::SplitReductionRule skip_all(bipartite);
    for (int v = 0; v < 8; ++v) {
        assert(skip_all.skip_vertex(v));
        const auto adjacent = bipartite.adjacent(v);
        auto subset = combutils::firstSubset(1, 2);
        do {
            std::uint64_t moved = 0;
            for (auto i : subset) moved |= std::uint64_t{1} << bipartite.half_edges[adjacent[i] ^ 1];
            assert(check_child(bipartite.splitGraph(v, adjacent, subset), v));
            assert(skip_all.redundant(v, moved));
        } while (combutils::nextSubset(subset, 3));
    }
    assert(checked && rejected && separating);
    std::cout << "Checked " << checked << " splits; " << rejected
              << " redundant; " << separating << " separating-edge cases\n";
}
