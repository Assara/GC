#include <cassert>
#include "graph.hpp"
#include "GraphGeneration/CutVertexSplitRule.hpp"

int main() {
    Graph<7,12,0,0,0,0,fieldType> graph;
    int edge = 0;
    for (auto block : {std::array{0,1,2,3}, std::array{0,4,5,6}})
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j) graph.setEdge(edge++, block[i], block[j]);
    GraphGeneration::CutVertexSplitRule<7> rule(graph);
    assert(rule.vertex == 0 && rule.branch_count == 2);
    const auto adjacent = graph.adjacent(0);
    int accepted = 0, rejected = 0;
    for (unsigned mask = 0; mask < 64; ++mask) {
        if (mask & 1 || std::popcount(mask) < 2 || std::popcount(mask) > 4) continue;
        std::vector<Int> subset;
        std::uint64_t neighbours = 0;
        for (int i = 0; i < 6; ++i) if (mask >> i & 1) {
            subset.push_back(i);
            neighbours |= std::uint64_t{1} << graph.half_edges[adjacent[i] ^ 1];
        }
        const auto child = graph.splitGraph(0, adjacent, subset);
        const bool resolved = child.is_connected_after_removing_vertex(0)
            && child.is_connected_after_removing_vertex(7);
        assert(rule.resolves(neighbours) == resolved);
        assert((GraphGeneration::CutVertexSplitRule<8>(child).vertex < 0) == child.is_biconnected());
        if (resolved) ++accepted; else ++rejected;
    }
    assert(accepted == 18 && rejected > 0);

    // A chain of three K4s: repairing one articulation leaves the other.
    Graph<10,18,0,0,0,0,fieldType> chain;
    edge = 0;
    for (auto block : {std::array{0,1,2,3}, std::array{0,4,5,6}, std::array{3,7,8,9}})
        for (int i = 0; i < 4; ++i)
            for (int j = i + 1; j < 4; ++j) chain.setEdge(edge++, block[i], block[j]);
    const GraphGeneration::CutVertexSplitRule<10> first(chain);
    assert(first.vertex == 0);
    const auto incidences = chain.adjacent(0);
    std::vector<Int> subset;
    for (Int i = 0; i < incidences.size(); ++i) {
        const auto neighbour = chain.half_edges[incidences[i] ^ 1];
        if (neighbour == 1 || neighbour == 4) subset.push_back(i);
    }
    assert(first.resolves((1ULL << 1) | (1ULL << 4)));
    const auto child = chain.splitGraph(0, incidences, subset);
    assert(child.is_connected_after_removing_vertex(0));
    assert(child.is_connected_after_removing_vertex(10));
    assert(GraphGeneration::CutVertexSplitRule<11>(child).vertex == 3);
}
