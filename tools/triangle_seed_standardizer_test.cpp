#include <algorithm>
#include <array>
#include <cassert>
#include <numeric>
#include <random>
#include <vector>
#include "GraphGeneration/TriangleSeedStandardizer.hpp"

int main() {
    using Stack = GraphGeneration::TriangleComponentStack<7>;
    using Entry = GraphGeneration::TriangleSeedEntry<7>;
    using Orientation = GraphGeneration::TriangleSeedStandardizer<7>::Orientation;
    GraphGeneration::TriangleSeedStandardizer<7> standardizer;
    const auto edges_of = [](const Entry& e) {
        return std::vector<Int>(e.canonical.begin(), e.canonical.begin() + e.endpoint_count);
    };
    const auto verify = [&](int vertices, std::vector<Int> edges, std::vector<std::uint64_t> masks) {
        Stack components;
        components.assign(masks);
        Orientation direction;
        const auto key = standardizer.standardize(vertices, edges, components, &direction);
        auto reversed_masks = masks;
        std::reverse(reversed_masks.begin(), reversed_masks.end());
        Stack reversed;
        reversed.assign(reversed_masks);
        Orientation reversed_direction;
        const auto reversed_key = standardizer.standardize(vertices, edges, reversed, &reversed_direction);
        if (masks.size() == 1) {
            assert(direction == Orientation::Keep && reversed_direction == Orientation::Keep);
            assert(!key.empty() && key == reversed_key);
        } else {
            if (direction == Orientation::Both) {
                assert(reversed_direction == Orientation::Both);
                assert(!key.empty() && key == reversed_key);
            } else {
                assert((direction == Orientation::Keep && reversed_direction == Orientation::Discard)
                    || (direction == Orientation::Discard && reversed_direction == Orientation::Keep));
                assert(key.empty() != reversed_key.empty());
            }
        }
        if (!key.empty()) {
            auto again = standardizer.standardize(vertices, edges_of(key), key.components);
            assert(key == again);
            assert(std::ranges::equal(key.components.masks(), again.components.masks()));
            for (std::size_t i = 0; i < masks.size(); ++i)
                assert(std::popcount(key.components.masks()[i]) == std::popcount(masks[i]));
        }
        std::vector<Int> labels(vertices);
        std::iota(labels.begin(), labels.end(), 0);
        std::mt19937 random(7189);
        for (int trial = 0; trial < 100; ++trial) {
            std::shuffle(labels.begin(), labels.end(), random);
            auto input = edges;
            for (auto& v : input) v = labels[v];
            std::reverse(input.begin(), input.end());
            const auto original = input;
            auto renamed = masks;
            for (auto& mask : renamed) {
                std::uint64_t result = 0;
                for (int v = 0; v < vertices; ++v)
                    if (mask & (std::uint64_t{1} << v)) result |= std::uint64_t{1} << labels[v];
                mask = result;
            }
            Stack reordered;
            reordered.assign(renamed);
            Orientation renamed_direction;
            const auto result = standardizer.standardize(vertices, input, reordered, &renamed_direction);
            assert(renamed_direction == direction);
            assert(result == key && input == original);
            assert(std::ranges::equal(result.components.masks(), key.components.masks()));
        }
    };
    verify(4, {0,1,0,2,1,2,1,3,2,3}, {15});
    // Equal diamonds sharing an articulation retain their input stack orientation.
    verify(7, {0,1,0,2,1,2,0,3,1,3,2,4,2,5,4,5,4,6,5,6}, {15,116});
    // Three equal branches sharing one vertex.
    verify(7, {0,1,0,2,1,2,0,3,0,4,3,4,0,5,0,6,5,6}, {7,25,97});
    // Both orientations must preserve their last component, with no size sorting.
    verify(6, {0,1,0,2,1,2,1,3,2,3,3,4,3,5,4,5}, {15,56});
    verify(6, {0,1,0,2,1,2,1,3,2,3,3,4,3,5,4,5}, {56,15});
    Orientation forward, backward;
    Stack stack;
    const std::vector<Int> unequal{0,1,0,2,1,2,1,3,2,3,3,4,3,5,4,5};
    stack.assign(std::array<std::uint64_t, 2>{15,56});
    assert(!standardizer.standardize(6, unequal, stack, &forward).empty());
    stack.assign(std::array<std::uint64_t, 2>{56,15});
    assert(standardizer.standardize(6, unequal, stack, &backward).empty());
    assert(forward == Orientation::Keep && backward == Orientation::Discard);

    const std::vector<Int> equal{0,1,0,2,1,2,0,3,1,3,2,4,2,5,4,5,4,6,5,6};
    stack.assign(std::array<std::uint64_t, 2>{15,116});
    const auto forward_key = standardizer.standardize(7, equal, stack, &forward);
    stack.assign(std::array<std::uint64_t, 2>{116,15});
    const auto backward_key = standardizer.standardize(7, equal, stack, &backward);
    assert(forward == Orientation::Both && backward == Orientation::Both);
    assert(!forward_key.empty() && forward_key == backward_key);

    // Equal-sized ends distinguished by how they attach to the middle diamond.
    const std::vector<Int> asymmetric{0,1,0,2,1,2,0,3,1,3,
        3,4,3,5,4,5,4,6,5,6,4,7,4,8,7,8,7,9,8,9};
    stack.assign(std::array<std::uint64_t, 3>{15,120,912});
    const auto a = standardizer.standardize(10, asymmetric, stack, &forward);
    stack.assign(std::array<std::uint64_t, 3>{912,120,15});
    const auto b = standardizer.standardize(10, asymmetric, stack, &backward);
    assert(forward != Orientation::Both && backward != Orientation::Both);
    assert(a.empty() != b.empty());
    verify(10, asymmetric, {15,120,912});
    verify(10, asymmetric, {912,120,15});
    // Equal outer sizes must not hide the first asymmetry farther inward.
    for (const auto& sizes : {std::vector<int>{3,4,3,3},
                              std::vector<int>{3,3,4,3,3,3}}) {
        int vertices = 1;
        int attachment = 0;
        std::vector<Int> edges;
        std::vector<std::uint64_t> masks;
        for (int size : sizes) {
            const int x = vertices++, y = vertices++;
            for (int v : {attachment,x,attachment,y,x,y}) edges.push_back(v);
            auto mask = (std::uint64_t{1} << attachment)
                      | (std::uint64_t{1} << x) | (std::uint64_t{1} << y);
            if (size == 4) {
                const int z = vertices++;
                for (int v : {x,z,y,z}) edges.push_back(v);
                mask |= std::uint64_t{1} << z;
                attachment = z;
            } else attachment = y;
            masks.push_back(mask);
        }
        stack.assign(masks);
        assert(!standardizer.standardize(vertices, edges, stack, &forward).empty());
        assert(forward == Orientation::Keep);
        verify(vertices, edges, masks);
        std::reverse(masks.begin(), masks.end());
        stack.assign(masks);
        assert(standardizer.standardize(vertices, edges, stack, &backward).empty());
        assert(backward == Orientation::Discard);
        verify(vertices, edges, masks);
    }
    // Exercise reversal through refinement and individualization on actual
    // triangle-addition histories, including stacks merged by reconnection.
    std::mt19937 histories(29183);
    for (int trial = 0; trial < 100; ++trial) {
        int vertices = 3;
        std::vector<Int> edges{0,1,0,2,1,2};
        Stack components;
        for (int step = 0; step < 6; ++step) {
            struct Addition { int a, b, c, vertices; std::vector<Int> edges; };
            std::vector<Addition> additions;
            for (int x = 0; x < vertices; ++x)
                for (int y = x + 1; y < vertices + 2; ++y)
                    for (int z = y + 1; z < vertices + 2; ++z) {
                        if (z == vertices + 1 && y != vertices) continue;
                        std::uint64_t existing = 0;
                        for (int v : {x,y,z}) if (v < vertices) existing |= std::uint64_t{1} << v;
                        if (!components.touches(existing)) continue;
                        auto next = edges;
                        for (auto [u,v] : {std::pair{x,y}, std::pair{x,z}, std::pair{y,z}}) {
                            bool found = false;
                            for (std::size_t e = 0; e < edges.size(); e += 2)
                                found |= edges[e] == u && edges[e + 1] == v;
                            if (!found) { next.push_back(u); next.push_back(v); }
                        }
                        const int count = std::max(vertices, z + 1);
                        if (next.size() == edges.size() || int(next.size()/2) - count + 1 > 7) continue;
                        additions.push_back({x,y,z,count,std::move(next)});
                    }
            if (additions.empty()) break;
            auto next = additions[histories() % additions.size()];
            std::uint64_t existing = 0, introduced = 0;
            for (int v : {next.a,next.b,next.c})
                (v < vertices ? existing : introduced) |= std::uint64_t{1} << v;
            components.add_triangle(existing, introduced, next.vertices - vertices);
            vertices = next.vertices;
            edges = std::move(next.edges);
        }
        std::vector<std::uint64_t> masks(components.masks().begin(), components.masks().end());
        verify(vertices, edges, masks);
        std::reverse(masks.begin(), masks.end());
        verify(vertices, edges, masks);
    }
    bool rejected = false;
    try { standardizer.standardize(4, std::vector<Int>{0,1,0}, Stack{}); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
}
