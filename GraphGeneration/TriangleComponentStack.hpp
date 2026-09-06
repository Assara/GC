#pragma once

#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>

namespace GraphGeneration {

// Each integer is a vertex mask in the graph's original, unstandardized labels.
// An attachment vertex belongs to both the parent and the new component.
template <std::size_t Capacity>
class TriangleComponentStack {
    static_assert(Capacity > 0);
    std::array<std::uint64_t, Capacity> components_{};
    std::size_t size_ = 1;

public:
    TriangleComponentStack() { components_[0] = 7; } // K3.

    std::size_t size() const { return size_; }
    std::uint64_t live() const { return components_[size_ - 1]; }
    bool touches(std::uint64_t existing_vertices) const {
        return (live() & existing_vertices) != 0;
    }

    void add_triangle(std::uint64_t existing_vertices, std::uint64_t new_vertices,
                      int new_vertex_count) {
        assert(touches(existing_vertices));
        // Reaching an older component merges the intervening stack frames.
        while ((existing_vertices & ~live()) != 0) {
            assert(size_ > 1);
            const auto popped = components_[--size_];
            components_[size_ - 1] |= popped;
        }
        if (new_vertex_count == 2) {
            assert(size_ < Capacity);
            components_[size_++] = existing_vertices | new_vertices;
        } else {
            components_[size_ - 1] |= new_vertices;
        }
    }
};

} // namespace GraphGeneration
