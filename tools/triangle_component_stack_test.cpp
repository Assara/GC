#include <cassert>
#include <cstdint>
#include "GraphGeneration/TriangleComponentStack.hpp"
#include "GraphGeneration/TriangleCompletionBudget.hpp"
#include "GraphGeneration/TriangleSeedEntry.hpp"
#include "LinearProbeSet.hpp"

int main() {
    using GraphGeneration::triangle_completion_fits;
    static_assert(triangle_completion_fits(3, 3, 1, 3)); // K3 can reach K4.
    static_assert(!triangle_completion_fits(7, 3, 5, 6));
    static_assert(triangle_completion_fits(7, 2, 5, 6)); // Tight deficit budget.
    static_assert(!triangle_completion_fits(7, 1, 6, 6));
    static_assert(triangle_completion_fits(7, 0, 6, 6)); // Cut vertices allowed.
    static_assert(triangle_completion_fits(5, 0, 6, 6)); // K5 can split to cubic.
    static_assert(!triangle_completion_fits(17, 0, 8, 9)); // Cubic vertex cap.
    using Stack = GraphGeneration::TriangleComponentStack<4>;
    const auto bit = [](int v) { return std::uint64_t{1} << v; };
    Stack graph;
    assert(graph.size() == 1 && graph.live() == 7);
    assert(!graph.last_exceeds_first());

    // Attach at vertex 2: the shared anchor stays live, older vertices freeze.
    graph.add_triangle(bit(2), bit(3) | bit(4), 2);
    assert(graph.size() == 2);
    assert(graph.live() == (bit(2) | bit(3) | bit(4)));
    assert(graph.touches(bit(2)));
    assert(!graph.touches(bit(0) | bit(1)));
    assert(!graph.last_exceeds_first()); // Equal triangles share their anchor.

    // Ordinary growth belongs to the live component without pushing.
    graph.add_triangle(bit(3) | bit(4), bit(5), 1);
    assert(graph.size() == 2 && graph.touches(bit(5)));
    assert(graph.last_exceeds_first()); // Diamond is larger than the old triangle.

    // Regression: two diamonds sharing a vertex have seven vertices, not eight.
    Stack diamonds;
    diamonds.add_triangle(bit(0) | bit(1), bit(3), 1);
    diamonds.add_triangle(bit(2), bit(4) | bit(5), 2);
    diamonds.add_triangle(bit(4) | bit(5), bit(6), 1);
    assert(diamonds.size() == 2 && !diamonds.last_exceeds_first());

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
    assert(!child.last_exceeds_first());

    // The declared capacity includes the initial K3 frame.
    Stack deepest;
    deepest.add_triangle(bit(2), bit(3) | bit(4), 2);
    deepest.add_triangle(bit(4), bit(5) | bit(6), 2);
    deepest.add_triangle(bit(6), bit(7) | bit(8), 2);
    assert(deepest.size() == 4);

    // Compare the ends, not the last against the total graph size.
    Stack ends;
    const std::array<std::uint64_t, 3> unequal_ends{7, 60, 480};
    ends.assign(unequal_ends); // Sizes 3, 4, 4; last is under half of all nine vertices.
    assert(ends.last_exceeds_first());
    const std::array<std::uint64_t, 3> equal_ends{15, 56, 480};
    ends.assign(equal_ends); // Sizes 4, 3, 4: equal ends pass, middle is irrelevant.
    assert(!ends.last_exceeds_first());
    assert(!ends.touches(bit(4))); // Only the last component is live.

    // Compare mirrored sizes inward; the first unequal pair decides.
    GraphGeneration::TriangleComponentStack<6> sizes;
    const auto check_sizes = [&](std::initializer_list<int> counts, bool discard) {
        std::array<std::uint64_t, 6> masks{};
        std::size_t n = 0;
        for (int count : counts) masks[n++] = (bit(count) - 1);
        sizes.assign({masks.data(), n});
        assert(sizes.reversed_sizes_are_greater() == discard);
    };
    check_sizes({3}, false);
    check_sizes({3, 4}, true);
    check_sizes({4, 3}, false);
    check_sizes({3, 3, 4, 3}, true);
    check_sizes({3, 4, 3, 3}, false);
    check_sizes({3, 3, 3, 4, 3, 3}, true);
    check_sizes({3, 3, 4, 3, 3, 3}, false);
    check_sizes({3, 4, 5, 4, 3}, false); // Odd middle cannot decide reversal.
    check_sizes({3, 4, 4, 3}, false);
    check_sizes({4, 3, 5, 3}, false); // Outer decision overrides inner sizes.
    check_sizes({3, 5, 3, 4}, true);

    using Entry = GraphGeneration::TriangleSeedEntry<4>;
    linear_probe_set<Entry> entries;
    Entry first;
    first.endpoint_count = 6;
    first.canonical = {0, 1, 0, 2, 1, 2};
    assert(entries.insert(first));
    auto duplicate = first;
    duplicate.components.add_triangle(bit(2), bit(3) | bit(4), 2);
    duplicate.canonical.back() = 9; // Unused capacity is not part of the key.
    assert(duplicate.hash() == first.hash() && duplicate == first);
    assert(!entries.insert(duplicate) && entries.size() == 1);

    // Force reallocation, then verify the original payload survived unchanged.
    for (int i = 1; i <= 40; ++i) {
        auto other = first;
        other.canonical[0] = static_cast<Int>(i);
        assert(entries.insert(other));
    }
    assert(entries.size() == 41 && entries.contains(duplicate));
    bool found = false;
    for (std::size_t i = 0; i < entries.capacity(); ++i) {
        const auto& stored = entries.data()[i];
        if (!stored.empty() && stored == first) {
            found = true;
            assert(stored.components.size() == 1);
        }
    }
    assert(found);

    // Completed graphs forget the stack and permit attachment at every vertex.
    auto completed = deepest;
    completed.forget();
    assert(completed.size() == 1 && completed.live() == 511);
    completed.forget(); // Idempotent.
    for (int v = 0; v < 9; ++v) assert(completed.touches(bit(v)));
    completed.add_triangle(bit(0), bit(9) | bit(10), 2);
    assert(completed.size() == 2 && completed.live() == (bit(0) | bit(9) | bit(10)));
    assert(!completed.touches(bit(1))); // Frozen old vertex.
    assert(completed.touches(bit(9)));
    completed.forget();
    assert(completed.size() == 1 && completed.live() == 2047);
}
