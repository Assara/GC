#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include "types.hpp"
#include "graph_hash.hpp"
#include "GraphGeneration/TriangleComponentStack.hpp"

namespace GraphGeneration {

template <std::size_t MaxLoop>
struct TriangleSeedEntry {
    // V <= 2L+1, E = V-1+L <= 3L. Fixed storage satisfies linear_probe_set's
    // noexcept value requirements and keeps payloads directly in its slots.
    std::array<Int, 6 * MaxLoop> canonical{};
    TriangleComponentStack<MaxLoop> components;
    std::size_t endpoint_count = 0;

    bool empty() const noexcept { return endpoint_count == 0; }
    std::size_t hash() const noexcept {
        return graph_hash_detail::hash_bytes(
            reinterpret_cast<const unsigned char*>(canonical.data()),
            endpoint_count * sizeof(Int));
    }
    bool operator==(const TriangleSeedEntry& other) const noexcept {
        return endpoint_count == other.endpoint_count
            && std::equal(canonical.begin(), canonical.begin() + endpoint_count,
                          other.canonical.begin());
    }
};

} // namespace GraphGeneration
