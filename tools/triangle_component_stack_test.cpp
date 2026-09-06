#include <cassert>
#include <cstdint>
#include "GraphGeneration/TriangleComponentStack.hpp"

int main() {
    using Stack = GraphGeneration::TriangleComponentStack<4>;
    const auto bit = [](int v) { return std::uint64_t{1} << v; };
    Stack graph;
    assert(graph.size() == 1 && graph.live() == 7);

    // Attach at vertex 2: the shared anchor stays live, older vertices freeze.
    graph.add_triangle(bit(2), bit(3) | bit(4), 2);
    assert(graph.size() == 2);
    assert(graph.live() == (bit(2) | bit(3) | bit(4)));
    assert(graph.touches(bit(2)));
    assert(!graph.touches(bit(0) | bit(1)));

    // Ordinary growth belongs to the live component without pushing.
    graph.add_triangle(bit(3) | bit(4), bit(5), 1);
    assert(graph.size() == 2 && graph.touches(bit(5)));

    // A sibling candidate must retain the parent's unchanged history.
    auto child = graph;
    child.add_triangle(bit(5), bit(6) | bit(7), 2);
    assert(child.size() == 3 && !child.touches(bit(3)));
    assert(graph.size() == 2 && graph.touches(bit(3)));

    // Reach the immediate parent: pop once and retain all merged vertices.
    auto parent_merge = child;
    parent_merge.add_triangle(bit(6) | bit(3), bit(8), 1);
    assert(parent_merge.size() == 2);
    assert(parent_merge.touches(bit(8)) && parent_merge.touches(bit(7)));
    assert(!parent_merge.touches(bit(0)));

    // Reach K3: pop through both frames, merging all labels.
    child.add_triangle(bit(6) | bit(0) | bit(1), 0, 0);
    assert(child.size() == 1 && child.live() == 255);

    // The declared capacity includes the initial K3 frame.
    Stack deepest;
    deepest.add_triangle(bit(2), bit(3) | bit(4), 2);
    deepest.add_triangle(bit(4), bit(5) | bit(6), 2);
    deepest.add_triangle(bit(6), bit(7) | bit(8), 2);
    assert(deepest.size() == 4);
}
