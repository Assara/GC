#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <cstring>
#include <span>
#include <stdexcept>
#include <unordered_set>
#include <utility>
#include <vector>
#include "types.hpp"

using namespace std;

namespace combutils {
	inline bool nextSubset(vector<Int>& S, Int A_max_index) {
		Int i = 0;
		for (; i < S.size(); ++i) {
			if (S[i] < A_max_index - i) {
				S[i]++;
				break;
			}
		}
		if (i == S.size()) return false;
		while(i > 0) {
			--i;
			S[i] = S[i+1] + 1;

		}
		return true;
	}

	inline vector<Int> firstSubset(Int startIndex, Int size) {
		vector<Int> subset(size);
		for (Int i = 0; i < size; ++i) {
			subset[i] = startIndex + (size - 1 - i);
		}
		return subset;
	}

	// Enumerate the Cartesian product
	//
	//     0 <= counts[i] <= upper_bounds[i].
	//
	// The count vector is reused between callbacks and is therefore exposed as
	// a read-only span.  An empty list of bounds has one element: the empty
	// count vector.
	template <typename Emit>
	std::uint64_t for_each_bounded_count_vector(
		std::span<const Int> upper_bounds,
		Emit&& emit
	) {
		std::vector<Int> counts(upper_bounds.size(), Int{0});
		std::uint64_t emitted = 0;
		while (true) {
			std::invoke(
				emit,
				std::span<const Int>(counts.data(), counts.size())
			);
			++emitted;

			std::size_t index = 0;
			while (index < counts.size()
				&& counts[index] == upper_bounds[index]) {
				counts[index] = 0;
				++index;
			}
			if (index == counts.size()) {
				break;
			}
			++counts[index];
		}
		return emitted;
	}

	// Counts describe how identical objects in each group are assigned to one
	// side of an unordered bipartition; group_sizes-counts describes the other
	// side.  Keep the lexicographically first of the complementary descriptions.
	// Equality is retained because a self-complementary description occurs only
	// once in the bounded count product.
	inline bool is_first_grouped_bipartition(
		std::span<const Int> counts,
		std::span<const Int> group_sizes
	) {
		if (counts.size() != group_sizes.size()) {
			throw std::invalid_argument("grouped bipartition size mismatch");
		}
		for (std::size_t i = 0; i < counts.size(); ++i) {
			if (counts[i] > group_sizes[i]) {
				throw std::invalid_argument("grouped bipartition count exceeds group size");
			}
			const Int complement = static_cast<Int>(group_sizes[i] - counts[i]);
			if (counts[i] < complement) {
				return true;
			}
			if (counts[i] > complement) {
				return false;
			}
		}
		return true;
	}


	inline bigInt n_splits(bigInt n_adjacent) {
		if (n_adjacent < 4) {
			return 0;
		}

		return (bigInt(1) << (n_adjacent - 1)) - n_adjacent - 1;
	}


	template<size_t N>
		inline signedInt compareHalfEdges(const std::array<Int, N>& a, const std::array<Int, N>& b) {
			return std::memcmp(a.data(), b.data(), N * sizeof(Int));

		}


	template <std::size_t N>
		[[nodiscard]] inline constexpr bool
		lessThan(const std::array<Int, N>& a, const std::array<Int, N>& b) noexcept {
			for (std::size_t i = 0; i < N; ++i) {
				if (a[i] < b[i]) return true;   // first differing element decides
				if (a[i] > b[i]) return false;
			}
			return false; // equal
		}

	template<typename T>
		vector<T> intersection(unordered_set<T> const& a,
				unordered_set<T> const& b) {
			vector<T> result;

			// iterate over the smaller set for efficiency
			auto const& small = (a.size() < b.size()) ? a : b;
			auto const& large = (a.size() < b.size()) ? b : a;

			for (auto const& x : small)
				if (large.contains(x))
					result.push_back(x);

			return result;
		}


	template<typename T>
		std::vector<T> ordered_intersection(std::vector<T> const& a,
				std::vector<T> const& b) {
			std::vector<T> result;
			result.reserve(std::min(a.size(), b.size()));

			auto it_a = a.begin();
			auto it_b = b.begin();

			while (it_a != a.end() && it_b != b.end()) {
				if (*it_a < *it_b) {
					++it_a;
				} else if (*it_b < *it_a) {
					++it_b;
				} else {
					// *it_a == *it_b
					result.push_back(*it_a);
					++it_a;
					++it_b;
				}
			}

			return result;
		}

}
