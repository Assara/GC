#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <span>

namespace GraphGeneration {

// Ordered construction stack. Growth stays at the last end; reconnecting to an
// earlier component pops and merges the intervening suffix. Never sort entries.
template <std::size_t Capacity>
class TriangleComponentStack {
    static_assert(Capacity > 0);
    std::array<std::uint64_t, Capacity> components_{};
    std::size_t size_ = 1;

public:
    TriangleComponentStack() noexcept { components_[0] = 7; }
    std::size_t size() const { return size_; }
    std::uint64_t live() const { return components_[size_ - 1]; }
    std::span<const std::uint64_t> masks() const { return {components_.data(), size_}; }
    void assign(std::span<const std::uint64_t> masks) noexcept {
        assert(!masks.empty() && masks.size() <= Capacity);
        size_ = masks.size();
        std::copy(masks.begin(), masks.end(), components_.begin());
    }
    bool last_exceeds_first() const noexcept {
        return std::popcount(live()) > std::popcount(components_[0]);
    }
    bool reversed_sizes_are_greater() const noexcept {
        for (std::size_t i = 0; i < size_ / 2; ++i) {
            const int first = std::popcount(components_[i]);
            const int last = std::popcount(components_[size_ - 1 - i]);
            if (first != last) return last > first;
        }
        return false;
    }
    bool touches(std::uint64_t existing_vertices) const {
        return (live() & existing_vertices) != 0;
    }

    // Once no bivalent vertices remain, all vertices become available for
    // growth again. This forgets construction history, not actual cut vertices.
    void forget() noexcept {
        for (std::size_t i = 1; i < size_; ++i) components_[0] |= components_[i];
        size_ = 1;
    }

    void add_triangle(std::uint64_t existing, std::uint64_t introduced, int new_count) {
        assert(touches(existing));
        while ((existing & ~live()) != 0) {
            assert(size_ > 1);
            const auto popped = components_[--size_];
            components_[size_ - 1] |= popped;
        }
        if (new_count == 2) {
            assert(size_ < Capacity);
            components_[size_++] = existing | introduced;
        } else {
            components_[size_ - 1] |= introduced;
        }
    }
};

} // namespace GraphGeneration
