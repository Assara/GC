#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <utility>
#include <vector>

#include "GraphGeneration/UnrootedSupportTransientGraph.hpp"

namespace GraphGeneration {

// Standardize4-style individualization/refinement for an unrooted simple
// support graph.  The first colors are the literal support valences: they are
// sorted and grouped before any hashing.  Consequently the final initial cell
// contains exactly the maximum-valence vertices, and its size is returned with
// the canonical graph for storage partitioning and split control.
template <Int N_VERTICES, Int N_EDGES>
class unrooted_support_transient_standardizer {
	public:
		using graph_type = unrooted_support_transient_graph<N_VERTICES, N_EDGES>;
		using adjacency_rows = std::array<std::uint64_t, N_VERTICES>;
		using permutation = std::array<Int, N_VERTICES>;
		using hash_type = std::uint64_t;

		struct canonicalization_result {
			graph_type canonical_graph;
			Int maximum_valence = 0;
			Int maximum_valence_count = 0;
			std::uint64_t automorphism_count = 1;
			std::size_t final_attempt_count = 1;
		};

		graph_type standardize_no_sign(const graph_type& graph) const {
			return standardize_with_info(graph).canonical_graph;
		}

		canonicalization_result standardize_with_info(
			const graph_type& input
		) const {
			const adjacency_rows adjacency = input.support().support_adjacency();
			const auto valences = input.support_valences();
			const Int maximum = *std::max_element(valences.begin(), valences.end());
			const Int maximum_count = static_cast<Int>(
				std::count(valences.begin(), valences.end(), maximum)
			);
			if constexpr (N_VERTICES == 1) {
				return {input, maximum, maximum_count, 1, 1};
			} else {
				auto [attempts, valid] = create_final_attempts(adjacency, valences);
				graph_type best = input.permuted(
					attempts[valid.front()].create_permutation()
				);
				std::uint64_t automorphisms = 1;
				for (std::size_t i = 1; i < valid.size(); ++i) {
					graph_type candidate = input.permuted(
						attempts[valid[i]].create_permutation()
					);
					const int comparison = best.support().compare(candidate.support());
					if (comparison < 0) {
						best = std::move(candidate);
						automorphisms = 1;
					} else if (comparison == 0) {
						++automorphisms;
					}
				}
				return {
					std::move(best), maximum, maximum_count,
					automorphisms, valid.size()
				};
			}
		}

	private:
		struct vertex_bucket {
			std::array<Int, N_VERTICES> data{};
			std::size_t size = 0;

			Int operator[](std::size_t index) const noexcept { return data[index]; }
			void push_back(Int value) noexcept { data[size++] = value; }

			int compare(const vertex_bucket& other) const noexcept {
				if (size < other.size) return -1;
				if (size > other.size) return 1;
				return size == 0 ? 0 : std::memcmp(
					data.data(), other.data.data(), size * sizeof(Int)
				);
			}
		};

		struct canon_builder {
			std::array<hash_type, N_VERTICES> colors{};
			std::array<Int, N_VERTICES> vertex_groups{};
			vertex_bucket group_separators;

			canon_builder() {
				std::iota(vertex_groups.begin(), vertex_groups.end(), Int{0});
			}

			void initialize_valence_groups(
				const std::array<Int, N_VERTICES>& valences
			) {
				for (std::size_t vertex = 0; vertex < N_VERTICES; ++vertex) {
					colors[vertex] = static_cast<hash_type>(valences[vertex]);
				}
				group_separators.size = 0;
				update_vertex_groups();
			}

			void refine(const adjacency_rows& adjacency) {
				std::array<hash_type, N_VERTICES> next{};
				for (std::size_t vertex = 0; vertex < N_VERTICES; ++vertex) {
					next[vertex] = hash(colors[vertex]);
					std::uint64_t neighbors = adjacency[vertex];
					while (neighbors != 0) {
						const std::size_t neighbor = std::countr_zero(neighbors);
						next[vertex] += colors[neighbor];
						neighbors &= neighbors - 1;
					}
				}
				colors = next;
			}

			void update_vertex_groups() {
				vertex_bucket next_separators;
				for (std::size_t group = 0;
					 group <= group_separators.size;
					 ++group) {
					const std::size_t begin = group == 0
						? 0 : static_cast<std::size_t>(group_separators[group - 1]);
					const std::size_t end = group < group_separators.size
						? static_cast<std::size_t>(group_separators[group])
						: N_VERTICES;
					if (group > 0) next_separators.push_back(static_cast<Int>(begin));
					if (begin + 1 >= end) continue;
					std::sort(
						vertex_groups.begin() + static_cast<std::ptrdiff_t>(begin),
						vertex_groups.begin() + static_cast<std::ptrdiff_t>(end),
						[this](Int first, Int second) {
							return colors[first] < colors[second];
						}
					);
					for (std::size_t separator = begin + 1;
						 separator < end;
						 ++separator) {
						if (colors[vertex_groups[separator]]
							!= colors[vertex_groups[separator - 1]]) {
							next_separators.push_back(static_cast<Int>(separator));
						}
					}
				}
				group_separators = next_separators;
			}

			int compare(const canon_builder& other) const noexcept {
				const int separator_comparison
					= group_separators.compare(other.group_separators);
				if (separator_comparison != 0) return separator_comparison;
				std::size_t begin = 0;
				for (std::size_t group = 0;
					 group <= group_separators.size;
					 ++group) {
					const std::size_t end = group < group_separators.size
						? static_cast<std::size_t>(group_separators[group])
						: N_VERTICES;
					if (begin < end) {
						const hash_type ours = colors[vertex_groups[begin]];
						const hash_type theirs = other.colors[other.vertex_groups[begin]];
						if (ours < theirs) return -1;
						if (ours > theirs) return 1;
					}
					begin = end;
				}
				return 0;
			}

			std::pair<std::size_t, std::size_t> branching_group() const {
				std::size_t begin = 0;
				for (std::size_t group = 0;
					 group <= group_separators.size;
					 ++group) {
					const std::size_t end = group < group_separators.size
						? static_cast<std::size_t>(group_separators[group])
						: N_VERTICES;
					if (end - begin > 1) return {begin, end};
					begin = end;
				}
				return {N_VERTICES, N_VERTICES};
			}

			void branch(
				std::pair<std::size_t, std::size_t> range,
				std::vector<canon_builder>& output
			) const {
				for (std::size_t chosen = range.first; chosen < range.second; ++chosen) {
					output.emplace_back(*this);
					canon_builder& child = output.back();
					++child.colors[child.vertex_groups[chosen]];
				}
			}

			permutation create_permutation() const noexcept {
				permutation result{};
				for (std::size_t new_vertex = 0;
					 new_vertex < N_VERTICES;
					 ++new_vertex) {
					result[static_cast<std::size_t>(vertex_groups[new_vertex])]
						= static_cast<Int>(new_vertex);
				}
				return result;
			}

			static hash_type hash(hash_type value) noexcept {
				value += 0x9e3779b97f4a7c15ULL;
				value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
				value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
				return value ^ (value >> 31);
			}
		};

		using final_attempt_set = std::pair<
			std::vector<canon_builder>, std::vector<std::size_t>
		>;

		final_attempt_set create_final_attempts(
			const adjacency_rows& adjacency,
			const std::array<Int, N_VERTICES>& valences
		) const {
			static constexpr std::size_t RELOAD_ITERATIONS
				= N_VERTICES > 2 ? N_VERTICES / 3 : 1;
			std::vector<canon_builder> attempts(1);
			std::vector<canon_builder> next_attempts;
			std::vector<std::size_t> valid{0};
			std::vector<std::size_t> next_valid;
			attempts.front().initialize_valence_groups(valences);

			while (attempts[valid.front()].group_separators.size
				< N_VERTICES - 1) {
				for (std::size_t reload = 0; reload < RELOAD_ITERATIONS; ++reload) {
					next_valid.clear();
					if (valid.size() == 1) {
						auto& attempt = attempts[valid.front()];
						attempt.refine(adjacency);
						attempt.update_vertex_groups();
						continue;
					}
					for (const std::size_t index : valid) {
						auto& attempt = attempts[index];
						attempt.refine(adjacency);
						attempt.update_vertex_groups();
						if (next_valid.empty()) {
							next_valid.push_back(index);
						} else {
							const int comparison = attempt.compare(
								attempts[next_valid.front()]
							);
							if (comparison > 0) next_valid.assign(1, index);
							else if (comparison == 0) next_valid.push_back(index);
						}
					}
					valid.swap(next_valid);
				}
				if (attempts[valid.front()].group_separators.size
					>= N_VERTICES - 1) break;
				const auto range = attempts[valid.front()].branching_group();
				next_attempts.clear();
				next_attempts.reserve((range.second - range.first) * valid.size());
				for (const std::size_t index : valid) {
					attempts[index].branch(range, next_attempts);
				}
				attempts.swap(next_attempts);
				valid.resize(attempts.size());
				std::iota(valid.begin(), valid.end(), std::size_t{0});
			}
			return {std::move(attempts), std::move(valid)};
		}
};

} // namespace GraphGeneration
