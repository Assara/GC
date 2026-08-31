#pragma once
#include "types.hpp"
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <bit>        
#include <functional>
static_assert(sizeof(std::size_t) == 8, "This hash expects 64-bit size_t.");

namespace graph_hash_detail {

inline std::size_t hash_bytes(const unsigned char* p, std::size_t len) noexcept {
	std::uint64_t h = 0x9E3779B185EBCA87ull ^ static_cast<std::uint64_t>(len);
	while (len >= 8) {
		std::uint64_t v;
		std::memcpy(&v, p, 8);
		h ^= v * 0x9E3779B185EBCA87ull;
		h = std::rotl(h, 27) * 0xC2B2AE3D27D4EB4Full;
		p += 8;
		len -= 8;
	}
	if (len) {
		std::uint64_t t = 0;
		for (std::size_t i = 0; i < len; ++i) {
			t |= static_cast<std::uint64_t>(p[i]) << (8 * i);
		}
		h ^= t * 0x9E3779B185EBCA87ull;
		h = std::rotl(h, 23) * 0xC2B2AE3D27D4EB4Full;
	}
	h ^= h >> 33; h *= 0xFF51AFD7ED558CCDull;
	h ^= h >> 33; h *= 0xC4CEB9FE1A85EC53ull;
	h ^= h >> 33;
	return static_cast<std::size_t>(h);
}

} // namespace graph_hash_detail

template <Int N_VERTICES, Int N_EDGES, Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename fieldType>
class Graph;

namespace std {


	template <Int N_VERTICES, Int N_EDGES, Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename fieldType>
		struct hash<Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>> {
			std::size_t operator()(const Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>& g) const noexcept {
				return g.hash();
			}
		};

} // namespace std
