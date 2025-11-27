#pragma once
#include "Types.hpp"
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <bit>        
#include <functional>
static_assert(sizeof(std::size_t) == 8, "This hash expects 64-bit size_t.");

template <Int N_VERTICES, Int N_EDGES, Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename fieldType>
class Graph;

namespace std {


template <Int N_VERTICES, Int N_EDGES, Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename fieldType>
struct hash<Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>> {
    std::size_t operator()(const Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>& g) const noexcept {
        // Treat half_edges (std::array<UInt, M>) as a byte buffer
        const unsigned char* p = reinterpret_cast<const unsigned char*>(g.half_edges.data());
        std::size_t len = sizeof(g.half_edges);

        // 64-bit accumulator; seed mixes length
        std::uint64_t h = 0x9E3779B185EBCA87ull ^ static_cast<std::uint64_t>(len);

        // Process 8 bytes at a time
        while (len >= 8) {
            std::uint64_t v;
            std::memcpy(&v, p, 8);                 // unaligned load, UB-free
            h ^= v * 0x9E3779B185EBCA87ull;
            h  = std::rotl(h, 27) * 0xC2B2AE3D27D4EB4Full;
            p  += 8;
            len -= 8;
        }

        // Tail (0..7 bytes)
        if (len) {
            std::uint64_t t = 0;
            switch (len) {
                case 7: t |= (std::uint64_t)p[6] << 48; [[fallthrough]];
                case 6: t |= (std::uint64_t)p[5] << 40; [[fallthrough]];
                case 5: t |= (std::uint64_t)p[4] << 32; [[fallthrough]];
                case 4: t |= (std::uint64_t)p[3] << 24; [[fallthrough]];
                case 3: t |= (std::uint64_t)p[2] << 16; [[fallthrough]];
                case 2: t |= (std::uint64_t)p[1] << 8;  [[fallthrough]];
                case 1: t |= (std::uint64_t)p[0];
                        h ^= t * 0x9E3779B185EBCA87ull;
                        h  = std::rotl(h, 23) * 0xC2B2AE3D27D4EB4Full;
                        break;
            }
        }

        // Final avalanche
        h ^= h >> 33; h *= 0xFF51AFD7ED558CCDull;
        h ^= h >> 33; h *= 0xC4CEB9FE1A85EC53ull;
        h ^= h >> 33;

        return static_cast<std::size_t>(h);
    }
};

} // namespace std
