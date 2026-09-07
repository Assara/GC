#pragma once

namespace GraphGeneration {

// Necessary conditions for completing by triangle additions and vertex splits
// to a connected cubic graph. Splits preserve loop order and cannot repair a
// bivalent vertex. One added loop repairs at most two bivalent vertices.
constexpr bool triangle_completion_fits(int vertices, int bivalent_vertices,
                                        int loop, int max_loop) noexcept {
    return loop <= max_loop
        && vertices <= 2 * (max_loop - 1)
        && bivalent_vertices <= 2 * (max_loop - loop);
}

} // namespace GraphGeneration
