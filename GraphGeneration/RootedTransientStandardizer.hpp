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

#include "GraphGeneration/RootedTransientGraph.hpp"

namespace GraphGeneration {

// Canonicalizer specialized to rooted_transient_graph.  Vertex zero is fixed;
// all refinement, individualization, and final permutations act only on the
// remaining vertices.  The implementation deliberately lives beside the
// general GraphStandardizer: the compact rooted representation has different
// invariants and should not complicate the sign-aware graph code.
template <
	Int N_VERTICES,
	Int N_EDGES,
	typename FieldType = fieldType
>
class rooted_transient_standardizer {
	public:
		using graph_type = rooted_transient_graph<N_VERTICES, N_EDGES, FieldType>;
		using hash_type = std::uint64_t;
		using vertex_type = Int;

		static constexpr std::size_t N_NONROOT
			= static_cast<std::size_t>(N_VERTICES - 1);
		static_assert(N_NONROOT <= 64,
			"rooted transient bit-mask refinement supports at most 64 nonroot vertices");
		using nonroot_permutation = std::array<vertex_type, N_NONROOT>;
		using core_adjacency_rows = std::array<std::uint64_t, N_NONROOT>;

		struct canonicalization_result {
			graph_type canonical_graph;
			// This is the order of the automorphism group fixing vertex zero.
			std::uint64_t automorphism_count = 1;
			std::size_t final_attempt_count = 1;
		};

		graph_type standardize_no_sign(const graph_type& input) const {
			return standardize_with_info(input).canonical_graph;
		}

		canonicalization_result standardize_with_info(
			const graph_type& input
		) const {
			if constexpr (N_NONROOT == 0) {
				return {input, 1, 1};
			} else {
				const core_adjacency_rows adjacency = create_adjacency_rows(input);
				auto [attempts, valid_attempts]
					= create_final_attempts(input, adjacency);

				graph_type best = input.permuted_nonroot_unchecked(
					attempts[valid_attempts.front()].create_nonroot_permutation()
				);
				std::uint64_t automorphism_count = 1;

				// Refinement colors are an efficient invariant, not the graph
				// certificate.  Always materialize and compare every survivor.  This
				// both resolves rare color collisions and gives the exact rooted
				// automorphism count essentially for free.
				for (std::size_t i = 1; i < valid_attempts.size(); ++i) {
					graph_type candidate = input.permuted_nonroot_unchecked(
						attempts[valid_attempts[i]].create_nonroot_permutation()
					);
					const int comparison = best.compare(candidate);
					if (comparison < 0) {
						best = std::move(candidate);
						automorphism_count = 1;
					} else if (comparison == 0) {
						++automorphism_count;
					}
				}

				return {
					std::move(best),
					automorphism_count,
					valid_attempts.size()
				};
			}
		}

	private:
		struct vertex_bucket {
			std::array<vertex_type, N_NONROOT> data{};
			std::size_t size = 0;

			vertex_type operator[](std::size_t index) const noexcept {
				return data[index];
			}

			void push_back(vertex_type vertex) noexcept {
				data[size++] = vertex;
			}

			int compare(const vertex_bucket& other) const noexcept {
				if (size < other.size) {
					return -1;
				}
				if (size > other.size) {
					return 1;
				}
				if (size == 0) {
					return 0;
				}
				return std::memcmp(
					data.data(), other.data.data(), size * sizeof(vertex_type)
				);
			}
		};

		struct canon_builder {
			std::array<hash_type, N_NONROOT> colors{};
			std::array<vertex_type, N_NONROOT> vertex_groups{};
			vertex_bucket group_separators;

			canon_builder() {
				std::iota(
					vertex_groups.begin(), vertex_groups.end(), vertex_type{0}
				);
			}

			void initialize_colors(
				const graph_type& graph,
				const core_adjacency_rows& adjacency
			) {
				for (std::size_t local = 0; local < N_NONROOT; ++local) {
					const vertex_type vertex
						= static_cast<vertex_type>(local + 1);
					// The exact size of the 0--vertex bundle is an initial vertex
					// weight.  It is not reduced modulo two.  Core degree is packed
					// alongside it, matching standardize4's degree starter while
					// preserving the bundle multiplicity as the leading attribute.
					const hash_type multiplicity = graph.root_multiplicity(vertex);
					const hash_type core_degree = std::popcount(adjacency[local]);
					colors[local] = hash(
						(multiplicity << 8) | core_degree
					);
				}
				group_separators.size = 0;
			}

			void update_colors(const core_adjacency_rows& adjacency) {
				std::array<hash_type, N_NONROOT> next_colors{};
				for (std::size_t vertex = 0; vertex < N_NONROOT; ++vertex) {
					next_colors[vertex] = hash(colors[vertex]);
					std::uint64_t neighbors = adjacency[vertex];
					while (neighbors != 0) {
						const unsigned neighbor = std::countr_zero(neighbors);
						next_colors[vertex] += colors[neighbor];
						neighbors &= neighbors - 1;
					}
				}
				colors = next_colors;
			}

			void update_vertex_groups() {
				vertex_bucket new_group_separators;

				for (std::size_t group = 0;
					 group <= group_separators.size;
					 ++group) {
					const std::size_t begin = group == 0
						? 0
						: static_cast<std::size_t>(group_separators[group - 1]);
					const std::size_t end = group < group_separators.size
						? static_cast<std::size_t>(group_separators[group])
						: N_NONROOT;

					if (group > 0) {
						new_group_separators.push_back(
							static_cast<vertex_type>(begin)
						);
					}
					if (begin + 1 >= end) {
						continue;
					}

					std::sort(
						vertex_groups.begin() + static_cast<std::ptrdiff_t>(begin),
						vertex_groups.begin() + static_cast<std::ptrdiff_t>(end),
						[this](vertex_type first, vertex_type second) {
							return colors[first] > colors[second];
						}
					);

					for (std::size_t separator = begin + 1;
						 separator < end;
						 ++separator) {
						if (colors[vertex_groups[separator]]
							!= colors[vertex_groups[separator - 1]]) {
							new_group_separators.push_back(
								static_cast<vertex_type>(separator)
							);
						}
					}
				}
				group_separators = new_group_separators;
			}

			int compare(const canon_builder& other) const noexcept {
				const int separator_comparison
					= group_separators.compare(other.group_separators);
				if (separator_comparison != 0) {
					return separator_comparison;
				}

				std::size_t begin = 0;
				std::size_t other_begin = 0;
				for (std::size_t group = 0;
					 group <= group_separators.size;
					 ++group) {
					const std::size_t end = group < group_separators.size
						? static_cast<std::size_t>(group_separators[group])
						: N_NONROOT;
					const std::size_t other_end
						= group < other.group_separators.size
						? static_cast<std::size_t>(other.group_separators[group])
						: N_NONROOT;

					if (begin < end && other_begin < other_end) {
						const hash_type color = colors[vertex_groups[begin]];
						const hash_type other_color
							= other.colors[other.vertex_groups[other_begin]];
						if (color < other_color) {
							return -1;
						}
						if (color > other_color) {
							return 1;
						}
					}
					begin = end;
					other_begin = other_end;
				}
				return 0;
			}

			std::pair<std::size_t, std::size_t> branching_group_range() const {
				std::size_t begin = 0;
				for (std::size_t group = 0;
					 group <= group_separators.size;
					 ++group) {
					const std::size_t end = group < group_separators.size
						? static_cast<std::size_t>(group_separators[group])
						: N_NONROOT;
					if (end - begin > 1) {
						return {begin, end};
					}
					begin = end;
				}
				return {N_NONROOT, N_NONROOT};
			}

			void branch(
				std::pair<std::size_t, std::size_t> range,
				std::vector<canon_builder>& collector
			) const {
				if (range.first >= range.second || range.second > N_NONROOT) {
					return;
				}
				collector.reserve(
					collector.size() + range.second - range.first
				);
				for (std::size_t chosen = range.first;
					 chosen < range.second;
					 ++chosen) {
					collector.emplace_back(*this);
					canon_builder& child = collector.back();
					++child.colors[child.vertex_groups[chosen]];
				}
			}

			nonroot_permutation create_nonroot_permutation() const noexcept {
				nonroot_permutation permutation{};
				for (std::size_t new_vertex = 0;
					 new_vertex < N_NONROOT;
					 ++new_vertex) {
					permutation[vertex_groups[new_vertex]]
						= static_cast<vertex_type>(new_vertex);
				}
				return permutation;
			}

			static hash_type hash(hash_type number) noexcept {
				number += 0x9e3779b97f4a7c15ULL;
				number = (number ^ (number >> 30))
					* 0xbf58476d1ce4e5b9ULL;
				number = (number ^ (number >> 27))
					* 0x94d049bb133111ebULL;
				return number ^ (number >> 31);
			}
		};

		using final_attempt_set = std::pair<
			std::vector<canon_builder>,
			std::vector<std::size_t>
		>;

		static core_adjacency_rows create_adjacency_rows(
			const graph_type& graph
		) {
			core_adjacency_rows result{};
			graph.for_each_core_edge([&](vertex_type full_first,
										 vertex_type full_second) {
				const std::size_t first
					= static_cast<std::size_t>(full_first - 1);
				const std::size_t second
					= static_cast<std::size_t>(full_second - 1);
				result[first] |= std::uint64_t{1} << second;
				result[second] |= std::uint64_t{1} << first;
			});
			return result;
		}

		final_attempt_set create_final_attempts(
			const graph_type& graph,
			const core_adjacency_rows& adjacency
		) const {
			static_assert(N_NONROOT > 0);
			static constexpr std::size_t RELOAD_ITERATIONS
				= N_NONROOT > 2 ? N_NONROOT / 3 : 1;

			std::vector<canon_builder> attempts(1);
			std::vector<canon_builder> next_attempts;
			std::vector<std::size_t> valid_attempts{0};
			std::vector<std::size_t> next_valid_attempts;
			attempts.front().initialize_colors(graph, adjacency);

			while (attempts[valid_attempts.front()].group_separators.size
				< N_NONROOT - 1) {
				for (std::size_t reload = 0;
					 reload < RELOAD_ITERATIONS;
					 ++reload) {
					next_valid_attempts.clear();
					next_valid_attempts.reserve(valid_attempts.size());

					if (valid_attempts.size() == 1) {
						canon_builder& attempt
							= attempts[valid_attempts.front()];
						attempt.update_colors(adjacency);
						attempt.update_vertex_groups();
						continue;
					}

					for (const std::size_t attempt_index : valid_attempts) {
						canon_builder& attempt = attempts[attempt_index];
						attempt.update_colors(adjacency);
						attempt.update_vertex_groups();

						if (next_valid_attempts.empty()) {
							next_valid_attempts.push_back(attempt_index);
							continue;
						}
						const int comparison = attempt.compare(
							attempts[next_valid_attempts.front()]
						);
						if (comparison > 0) {
							next_valid_attempts.assign(1, attempt_index);
						} else if (comparison == 0) {
							next_valid_attempts.push_back(attempt_index);
						}
					}
					valid_attempts.swap(next_valid_attempts);
				}

				if (attempts[valid_attempts.front()].group_separators.size
					>= N_NONROOT - 1) {
					break;
				}

				const auto branch_range
					= attempts[valid_attempts.front()].branching_group_range();
				next_attempts.clear();
				next_attempts.reserve(
					(branch_range.second - branch_range.first)
					* valid_attempts.size()
				);
				for (const std::size_t attempt_index : valid_attempts) {
					attempts[attempt_index].branch(
						branch_range, next_attempts
					);
				}

				attempts.swap(next_attempts);
				valid_attempts.resize(attempts.size());
				std::iota(
					valid_attempts.begin(), valid_attempts.end(), std::size_t{0}
				);
			}

			return {std::move(attempts), std::move(valid_attempts)};
		}
};

} // namespace GraphGeneration
