#pragma once

#include <algorithm>
#include <array>
#include <cstddef>
#include "types.hpp"
#include "graph_hash.hpp"

namespace GraphGeneration {

// A graph alone determines every possible continuation. Store no construction
// history or second labelling; expand directly from these canonical edges.
template <std::size_t MaxLoop>
struct TriangleGraphEntry {
    std::array<Int, 6 * MaxLoop> canonical{};
    std::size_t endpoint_count = 0;

    bool empty() const noexcept { return endpoint_count == 0; }
    std::size_t hash() const noexcept {
        return graph_hash_detail::hash_bytes(
            reinterpret_cast<const unsigned char*>(canonical.data()),
            endpoint_count * sizeof(Int));
    }
    bool operator==(const TriangleGraphEntry& other) const noexcept {
        return endpoint_count == other.endpoint_count
            && std::equal(canonical.begin(), canonical.begin() + endpoint_count,
                          other.canonical.begin());
    }
};

} // namespace GraphGeneration
