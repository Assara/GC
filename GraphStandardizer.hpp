#pragma once

#include <cstring>
#include <cstdint>
#include <algorithm>
#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
#include <atomic>
#include <chrono>
#endif
#include <vector>
#include <utility>
#include "VectorSpace/BasisElement.hpp"
#include "GraphIsomorphism.hpp"

#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
namespace gc_standardizer_sort_profile {
	inline std::atomic<unsigned long long> insertion_sort_calls{0};
	inline std::atomic<unsigned long long> insertion_sort_swaps{0};
	inline std::atomic<unsigned long long> insertion_sort_nanoseconds{0};
	inline std::atomic<unsigned long long> vertex_color_update_calls{0};
	inline std::atomic<unsigned long long> vertex_color_update_nanoseconds{0};
	inline std::atomic<unsigned long long> vertex_bucket_init_calls{0};
	inline std::atomic<unsigned long long> vertex_bucket_init_nanoseconds{0};
	inline std::atomic<unsigned long long> mask_maintenance_calls{0};
	inline std::atomic<unsigned long long> mask_maintenance_nanoseconds{0};
	inline std::atomic<unsigned long long> push_next_attempts_calls{0};
	inline std::atomic<unsigned long long> push_next_attempts_nanoseconds{0};
	inline std::atomic<unsigned long long> labeling_search_calls{0};
	inline std::atomic<unsigned long long> labeling_search_nanoseconds{0};
	inline std::atomic<unsigned long long> attempts_created{0};
	inline std::atomic<unsigned long long> hash_calls{0};
	inline std::atomic<unsigned long long> hash_nanoseconds{0};

	inline void reset() {
		insertion_sort_calls.store(0, std::memory_order_relaxed);
		insertion_sort_swaps.store(0, std::memory_order_relaxed);
		insertion_sort_nanoseconds.store(0, std::memory_order_relaxed);
		vertex_color_update_calls.store(0, std::memory_order_relaxed);
		vertex_color_update_nanoseconds.store(0, std::memory_order_relaxed);
		vertex_bucket_init_calls.store(0, std::memory_order_relaxed);
		vertex_bucket_init_nanoseconds.store(0, std::memory_order_relaxed);
		mask_maintenance_calls.store(0, std::memory_order_relaxed);
		mask_maintenance_nanoseconds.store(0, std::memory_order_relaxed);
		push_next_attempts_calls.store(0, std::memory_order_relaxed);
		push_next_attempts_nanoseconds.store(0, std::memory_order_relaxed);
		labeling_search_calls.store(0, std::memory_order_relaxed);
		labeling_search_nanoseconds.store(0, std::memory_order_relaxed);
		attempts_created.store(0, std::memory_order_relaxed);
		hash_calls.store(0, std::memory_order_relaxed);
		hash_nanoseconds.store(0, std::memory_order_relaxed);
	}
}
#endif

template <
Int N_VERTICES,
	Int N_EDGES,
	Int N_OUT_HAIR,
	Int N_IN_HAIR,
	signedInt c,
	signedInt d,
	typename fieldType
	>
	class Graph;

	template <
	Int N_VERTICES,
	Int N_EDGES,
	Int N_OUT_HAIR,
	Int N_IN_HAIR,
	signedInt c,
	signedInt d,
	typename fieldType
	>
	class GraphStandardizer {
		public:
			using GraphType = Graph<N_VERTICES,N_EDGES,N_OUT_HAIR, N_IN_HAIR,c,d, fieldType>;
			using IsomorphismType = GraphIsomorphism<N_VERTICES, N_EDGES>;

			using hash_int_type = std::uint64_t;
			using assign_type = Int;


			struct VertexBucket {
				array<assign_type, N_VERTICES> data{};
				std::size_t size = 0;

				assign_type operator[](std::size_t i) const {
					return data[i];
				}

				assign_type& operator[](std::size_t i) {
					return data[i];
				}

				void push_back(assign_type v) {
					data[size++] = v;
				}

				void clean_and_push(assign_type v) {
					size = 1;
					data[0] = v;
				}

				bool empty() const {
					return size == 0;
				}

				signedInt compare(const VertexBucket& other) const {
					if (size < other.size) {
						return -1;
					}
					if (size > other.size) {
						return 1;
					}

					return std::memcmp(
						data.data(),
						other.data.data(),
						size * sizeof(assign_type)
					);
				}
			};

			struct CanonBuilder3 {
				array<hash_int_type, N_VERTICES> colors;
				Permutation<N_VERTICES> vertex_permutation;
				array<assign_type, N_VERTICES> vertex_groups;
				VertexBucket group_separators;
				assign_type next_to_assign = 0;
				assign_type active_size = N_VERTICES;

				CanonBuilder3()
					: colors{},
					  vertex_permutation{},
					  vertex_groups{},
					  group_separators{},
					  next_to_assign(0),
					  active_size(N_VERTICES) {
					vertex_permutation.p.fill(N_VERTICES);
					std::iota(vertex_groups.begin(), vertex_groups.end(), 0);
				}

				void update_vertex_groups() {
					array<assign_type, N_VERTICES> new_vertex_groups{};
					VertexBucket new_group_separators;
					std::size_t write = 0;
					std::size_t begin = 0;

					for (std::size_t sep_i = 0; sep_i <= group_separators.size; ++sep_i) {
						const std::size_t end =
							(sep_i < group_separators.size) ? group_separators[sep_i] : active_size;

						if (end - begin <= 1) {
							begin = end;
							continue;
						}

						std::sort(
							vertex_groups.begin() + begin,
							vertex_groups.begin() + end,
							[this](assign_type a, assign_type b) {
								return colors[a] > colors[b];
							}
						);

						for (std::size_t run_begin = begin; run_begin < end;) {
							std::size_t run_end = run_begin + 1;
							while (
								run_end < end &&
								colors[vertex_groups[run_begin]] == colors[vertex_groups[run_end]]
							) {
								++run_end;
							}

							if (run_end - run_begin == 1) {
								const assign_type v = vertex_groups[run_begin];
								if (vertex_permutation[v] == N_VERTICES) {
									vertex_permutation.p[v] = next_to_assign;
									++colors[v];
									++next_to_assign;
								}
							} else {
								for (std::size_t i = run_begin; i < run_end; ++i) {
									new_vertex_groups[write++] = vertex_groups[i];
								}

								if (write < N_VERTICES) {
									new_group_separators.push_back(static_cast<assign_type>(write));
								}
							}

							run_begin = run_end;
						}

						begin = end;
					}

					vertex_groups = new_vertex_groups;
					group_separators = new_group_separators;
					active_size = static_cast<assign_type>(write);
				}

				void init_starter_colors(const GraphType& G) {
					array<assign_type, (N_VERTICES > 0) ? (N_VERTICES - 1) : 0> valency_counts{};
					for (Int e = 0; e < G.half_edges.size(); ++e) {
						++colors[G.half_edges[e]];
					}

					for (assign_type v = 0; v < N_VERTICES; ++v) {
						if (colors[v] > 0 && colors[v] < N_VERTICES) {
							++valency_counts[colors[v] - 1];
						}
					}

					if constexpr (GraphType::N_HAIR > 0) {
						for (Int i = 0; i < GraphType::N_HAIR; ++i) {
							colors[G.half_edges[i]] += hash(i);
						}
					}

					next_to_assign = 0;
					for (std::size_t valency = N_VERTICES; valency-- > 1;) {
						if (valency_counts[valency - 1] != 1) {
							continue;
						}

						for (assign_type v = 0; v < N_VERTICES; ++v) {
							if (colors[v] == valency) {
								vertex_permutation.p[v] = next_to_assign;
								++colors[v];
								++next_to_assign;
								break;
							}
						}
					}
				}

				signedInt compare(const CanonBuilder3& other) const {
					if (next_to_assign < other.next_to_assign) {
						return -1;
					}
					if (next_to_assign > other.next_to_assign) {
						return 1;
					}

					const signedInt separator_comparison = group_separators.compare(other.group_separators);
					if (separator_comparison != 0) {
						return separator_comparison;
					}

					assign_type begin = 0;
					assign_type other_begin = 0;
					for (std::size_t sep_i = 0; sep_i <= group_separators.size; ++sep_i) {
						const assign_type end =
							(sep_i < group_separators.size) ? group_separators[sep_i] : active_size;
						const assign_type other_end =
							(sep_i < other.group_separators.size) ? other.group_separators[sep_i] : other.active_size;

						if (end > begin && other_end > other_begin) {
							const hash_int_type color = colors[vertex_groups[begin]];
							const hash_int_type other_color = other.colors[other.vertex_groups[other_begin]];
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

				std::pair<assign_type, assign_type> branching_group_range() const {
					assign_type begin = 0;
					for (std::size_t sep_i = 0; sep_i <= group_separators.size; ++sep_i) {
						const assign_type end =
							(sep_i < group_separators.size) ? group_separators[sep_i] : active_size;
						if (end - begin > 1) {
							return {begin, end};
						}
						begin = end;
					}
					return {N_VERTICES, N_VERTICES};
				}

				void assign_singleton_groups() {
					assign_type begin = 0;
					for (std::size_t sep_i = 0; sep_i <= group_separators.size; ++sep_i) {
						const assign_type end =
							(sep_i < group_separators.size) ? group_separators[sep_i] : active_size;

						if (end - begin == 1) {
							const assign_type v = vertex_groups[begin];
							if (vertex_permutation[v] == N_VERTICES) {
								vertex_permutation.p[v] = next_to_assign;
								++colors[v];
								++next_to_assign;
							}
						}

						begin = end;
					}
				}

				void assign_remaining_groups_in_order() {
					assign_type begin = 0;
					for (std::size_t sep_i = 0; sep_i <= group_separators.size; ++sep_i) {
						const assign_type end =
							(sep_i < group_separators.size) ? group_separators[sep_i] : active_size;

						for (assign_type i = begin; i < end; ++i) {
							const assign_type v = vertex_groups[i];
							if (vertex_permutation[v] == N_VERTICES) {
								vertex_permutation.p[v] = next_to_assign;
								++next_to_assign;
							}
						}

						begin = end;
					}
				}

				bool all_unassigned_vertices_are_singleton_groups() const {
					assign_type singleton_unassigned_count = 0;
					assign_type begin = 0;
					for (std::size_t sep_i = 0; sep_i <= group_separators.size; ++sep_i) {
						const assign_type end =
							(sep_i < group_separators.size) ? group_separators[sep_i] : active_size;
						const assign_type size = end - begin;

						if (size > 1) {
							return false;
						}
						if (size == 1) {
							const assign_type v = vertex_groups[begin];
							if (vertex_permutation[v] == N_VERTICES) {
								++singleton_unassigned_count;
							}
						}

						begin = end;
					}

					assign_type total_unassigned_count = 0;
					for (assign_type v = 0; v < N_VERTICES; ++v) {
						if (vertex_permutation[v] == N_VERTICES) {
							++total_unassigned_count;
						}
					}

					return singleton_unassigned_count == total_unassigned_count;
				}

				void branch(
					std::pair<assign_type, assign_type> branch_range,
					vector<CanonBuilder3>& collector
				) const {
					const assign_type begin = branch_range.first;
					const assign_type end = branch_range.second;
					if (begin >= end || end > active_size) {
						return;
					}

					const assign_type branch_size = end - begin;
					collector.reserve(collector.size() + branch_size);

					for (assign_type chosen = begin; chosen < end; ++chosen) {
						collector.emplace_back(*this);
						CanonBuilder3& child = collector.back();
						++child.colors[child.vertex_groups[chosen]];
					}
				}

				void update_colors(const GraphType& G) {
					array<hash_int_type, N_VERTICES> next_colors;

					for (Int v = 0; v < N_VERTICES; ++v) {
						next_colors[v] = hash(colors[v]);
					}

					for (Int e = G.N_HAIR; e < G.half_edges.size(); e += 2) {
						next_colors[G.half_edges[e]] += colors[G.half_edges[e + 1]];
						next_colors[G.half_edges[e + 1]] += colors[G.half_edges[e]];
					}

					colors = next_colors;
				}

				static hash_int_type hash(hash_int_type n) noexcept {
					n += 0x9e3779b97f4a7c15ULL;
					n = (n ^ (n >> 30)) * 0xbf58476d1ce4e5b9ULL;
					n = (n ^ (n >> 27)) * 0x94d049bb133111ebULL;
					return n ^ (n >> 31);
				}
			};

			BasisElement<GraphType, fieldType> standardize3(const BasisElement<GraphType, fieldType>& input) const {
				static constexpr Int RELOAD_ITERATIONS = (N_VERTICES > 2) ? (N_VERTICES / 3) : 1;

#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
				const auto labeling_search_start = std::chrono::steady_clock::now();
				unsigned long long local_attempts_created = 1;
#endif
#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
				unsigned long long local_update_colors_calls = 0;
				unsigned long long local_update_groups_calls = 0;
#endif
#if defined(GC_PROFILE_STANDARDIZER_SORT)
				unsigned long long local_update_colors_nanoseconds = 0;
				unsigned long long local_update_groups_nanoseconds = 0;
#endif
				vector<CanonBuilder3> attempts;
				vector<CanonBuilder3> next_attempts;
				const GraphType& G = input.getValue();

				attempts.emplace_back(CanonBuilder3());
				attempts[0].init_starter_colors(G);

				signedInt cmp;
				while (attempts.back().next_to_assign < N_VERTICES) {
					for (Int reload = 0; reload < RELOAD_ITERATIONS; ++reload) {
						next_attempts.clear();
						next_attempts.reserve(attempts.size());

						for (size_t i = 0; i < attempts.size(); ++i) {
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							const auto update_colors_start = std::chrono::steady_clock::now();
#endif
							attempts[i].update_colors(G);
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							const auto update_colors_stop = std::chrono::steady_clock::now();
							local_update_colors_nanoseconds += static_cast<unsigned long long>(
								std::chrono::duration_cast<std::chrono::nanoseconds>(
									update_colors_stop - update_colors_start
								).count()
							);
#endif
#if defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
							++local_update_colors_calls;
#endif
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							++local_update_colors_calls;
#endif
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							const auto update_groups_start = std::chrono::steady_clock::now();
#endif
							attempts[i].update_vertex_groups();
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							const auto update_groups_stop = std::chrono::steady_clock::now();
							local_update_groups_nanoseconds += static_cast<unsigned long long>(
								std::chrono::duration_cast<std::chrono::nanoseconds>(
									update_groups_stop - update_groups_start
								).count()
							);
#endif
#if defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
							++local_update_groups_calls;
#endif
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							++local_update_groups_calls;
#endif

							if (next_attempts.empty()) {
								next_attempts.push_back(attempts[i]);
							} else {
								cmp = attempts[i].compare(next_attempts.back());
								if (cmp > 0) {
									next_attempts[0] = attempts[i];
									next_attempts.resize(1);
								} else if (cmp == 0) {
									next_attempts.push_back(attempts[i]);
								}
							}
						}
						std::swap(attempts, next_attempts);
					}

					if (attempts[0].next_to_assign >= N_VERTICES) {
						break;
					}

					const auto branch_range = attempts[0].branching_group_range();
					if (branch_range.first >= branch_range.second) {
						for (auto& attempt : attempts) {
							attempt.assign_remaining_groups_in_order();
						}
						break;
					}

					next_attempts.clear();
					next_attempts.reserve(
						static_cast<std::size_t>(branch_range.second - branch_range.first) * attempts.size()
					);

					for (const auto& attempt : attempts) {
						attempt.branch(branch_range, next_attempts);
					}
#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
					local_attempts_created += next_attempts.size();
#endif

					attempts.swap(next_attempts);
				}

#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
				gc_standardizer_sort_profile::vertex_color_update_calls.fetch_add(
					local_update_colors_calls,
					std::memory_order_relaxed
				);
				gc_standardizer_sort_profile::vertex_bucket_init_calls.fetch_add(
					local_update_groups_calls,
					std::memory_order_relaxed
				);
				gc_standardizer_sort_profile::attempts_created.fetch_add(
					local_attempts_created,
					std::memory_order_relaxed
				);
#endif
#if defined(GC_PROFILE_STANDARDIZER_SORT)
				gc_standardizer_sort_profile::vertex_color_update_nanoseconds.fetch_add(
					local_update_colors_nanoseconds,
					std::memory_order_relaxed
				);
				gc_standardizer_sort_profile::vertex_bucket_init_nanoseconds.fetch_add(
					local_update_groups_nanoseconds,
					std::memory_order_relaxed
				);
#endif
#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
				const auto labeling_search_stop = std::chrono::steady_clock::now();
				gc_standardizer_sort_profile::labeling_search_calls.fetch_add(
					1,
					std::memory_order_relaxed
				);
				gc_standardizer_sort_profile::labeling_search_nanoseconds.fetch_add(
					static_cast<unsigned long long>(
						std::chrono::duration_cast<std::chrono::nanoseconds>(
							labeling_search_stop - labeling_search_start
						).count()
					),
					std::memory_order_relaxed
				);
#endif

				GraphType best_graph = input.getValue();
				bool have_best = false;
				bool containsPlus = false;
				bool containsMinus = false;
				bool checked_double_edge = false;

				for (const CanonBuilder3& attempt : attempts) {
					GraphType graph = input.getValue();
					signedInt sign = graph.permuteVertices(attempt.vertex_permutation);
					sign *= graph.directAndSortEdges();

					if constexpr (GraphType::SWAP_EDGE_SIGN == -1) {
						if (!checked_double_edge) {
							checked_double_edge = true;
							if (graph.has_double_edge()) {
								return BasisElement<GraphType, fieldType>(graph, 0);
							}
						}
					}

					const signedInt comparison = have_best ? graph.compare(best_graph) : -1;

					if (comparison < 0) {
						best_graph = graph;
						have_best = true;
						containsPlus = sign > 0;
						containsMinus = sign < 0;
					} else if (comparison == 0) {
						containsPlus = containsPlus || sign > 0;
						containsMinus = containsMinus || sign < 0;
					}
				}

				if (containsPlus && containsMinus) {
					return BasisElement<GraphType, fieldType>(best_graph, 0);
				}
				if (containsPlus) {
					return BasisElement<GraphType, fieldType>(best_graph, input.getCoefficient());
				}
				return BasisElement<GraphType, fieldType>(best_graph, -input.getCoefficient());
			}



			struct CanonBuilder2 {
				array<hash_int_type, N_VERTICES> colors;
				Permutation<N_VERTICES> vertex_permutation;

				CanonBuilder2()
					: colors{}, vertex_permutation{} {
					vertex_permutation.p.fill(N_VERTICES);
				}

				assign_type init_starter_colors(const GraphType& G) {
					array<assign_type, (N_VERTICES > 0) ? (N_VERTICES - 1) : 0> valency_counts{};

					for (Int e = 0; e < G.half_edges.size(); ++e) {
						++colors[G.half_edges[e]];
					}

					for (assign_type v = 0; v < N_VERTICES; ++v) {
						if (colors[v] > 0 && colors[v] < N_VERTICES) {
							++valency_counts[colors[v] - 1];
						}
					}

					if constexpr (GraphType::N_HAIR > 0) {
						for (Int i = 0; i < GraphType::N_HAIR; ++i) {
							colors[G.half_edges[i]] += hash(i);
						}
					}

					assign_type next_vertex_to_assign = 0;
					for (std::size_t valency = N_VERTICES; valency-- > 1;) {
						if (valency_counts[valency - 1] != 1) {
							continue;
						}

						for (assign_type v = 0; v < N_VERTICES; ++v) {
							if (colors[v] == valency) {
								vertex_permutation.p[v] = next_vertex_to_assign;
								++colors[v];
								++next_vertex_to_assign;
								break;
							}
						}
					}

					return next_vertex_to_assign;
				}

				void update_colors(const GraphType& G) {
					array<hash_int_type, N_VERTICES> next_colors;
#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
					unsigned long long hash_count = 0;
#endif
#if defined(GC_PROFILE_STANDARDIZER_SORT)
					const auto hash_start = std::chrono::steady_clock::now();
#endif
					for (Int v = 0; v < N_VERTICES; ++v) {
						next_colors[v] = hash(colors[v]);
#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
						++hash_count;
#endif
					}
					for (Int e = G.N_HAIR; e<G.half_edges.size(); e+=2) {
						next_colors[G.half_edges[e]] += colors[G.half_edges[e+1]];
						next_colors[G.half_edges[e+1]] += colors[G.half_edges[e]];
					}
#if defined(GC_PROFILE_STANDARDIZER_SORT)
					const auto hash_stop = std::chrono::steady_clock::now();
					gc_standardizer_sort_profile::hash_calls.fetch_add(
						hash_count,
						std::memory_order_relaxed
					);
					gc_standardizer_sort_profile::hash_nanoseconds.fetch_add(
						static_cast<unsigned long long>(
							std::chrono::duration_cast<std::chrono::nanoseconds>(
								hash_stop - hash_start
							).count()
						),
						std::memory_order_relaxed
					);
#endif
#if defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
					gc_standardizer_sort_profile::hash_calls.fetch_add(
						hash_count,
						std::memory_order_relaxed
					);
#endif
					colors = next_colors;
				}

				struct CanonicalBucketAnalysis {
					VertexBucket branch_bucket{};
					assign_type unique_count = 0;
					hash_int_type prune_hash = 0;
				};

				CanonicalBucketAnalysis assign_unique_colors_and_make_bucket(assign_type next_vertex_to_assign) {
					array<std::pair<hash_int_type, assign_type>, N_VERTICES> unassigned_colors{};
					std::size_t n_unassigned = 0;
					for (assign_type v = 0; v < N_VERTICES; ++v) {
						if (vertex_permutation[v] == N_VERTICES) {
							unassigned_colors[n_unassigned++] = {colors[v], v};
						}
					}

					std::sort(
						unassigned_colors.begin(),
						unassigned_colors.begin() + n_unassigned,
						[](const auto& lhs, const auto& rhs) {
							if (lhs.first != rhs.first) {
								return lhs.first > rhs.first;
							}
							return lhs.second < rhs.second;
						}
					);

					CanonicalBucketAnalysis analysis;
					for (std::size_t begin = 0; begin < n_unassigned;) {
						std::size_t end = begin + 1;
						while (end < n_unassigned && unassigned_colors[end].first == unassigned_colors[begin].first) {
							++end;
						}

						const assign_type run_size = static_cast<assign_type>(end - begin);
						if (run_size == 1) {
							const assign_type vertex = unassigned_colors[begin].second;
							vertex_permutation.p[vertex] = next_vertex_to_assign + analysis.unique_count;
							++colors[vertex];
							++analysis.unique_count;
						} else {
							for (std::size_t i = begin; i < end; ++i) {
								analysis.prune_hash = hash(analysis.prune_hash ^ unassigned_colors[i].first);
							}
							if (
								analysis.branch_bucket.empty() ||
								run_size < static_cast<assign_type>(analysis.branch_bucket.size)
							) {
								analysis.branch_bucket.size = 0;
								for (std::size_t i = begin; i < end; ++i) {
									analysis.branch_bucket.push_back(unassigned_colors[i].second);
								}
							}
						}

						begin = end;
					}

					return analysis;
				}

				void push_next_attempts(
					vector<CanonBuilder2>& collector,
					const VertexBucket& branch_bucket,
					assign_type v
				) const {
					if (branch_bucket.size == 2) {
						collector.emplace_back(*this);
						collector.back().vertex_permutation.p[branch_bucket[0]] = v;
						collector.back().vertex_permutation.p[branch_bucket[1]] = v+1;
						++collector.back().colors[branch_bucket[0]];

						collector.emplace_back(*this);
						collector.back().vertex_permutation.p[branch_bucket[1]] = v;
						collector.back().vertex_permutation.p[branch_bucket[0]] = v+1;
						++collector.back().colors[branch_bucket[1]];

						return;
					}

					for (std::size_t i = 0; i < branch_bucket.size; ++i) {
						collector.emplace_back(*this);
						collector.back().vertex_permutation.p[branch_bucket[i]] = v;
						++collector.back().colors[branch_bucket[i]];
					}
				}



				static hash_int_type hash(hash_int_type n) noexcept {
					n += 0x9e3779b97f4a7c15ULL;
					n = (n ^ (n >> 30)) * 0xbf58476d1ce4e5b9ULL;
					n = (n ^ (n >> 27)) * 0x94d049bb133111ebULL;
					return n ^ (n >> 31);
				}
			};

			BasisElement<GraphType, fieldType> standardize2(BasisElement<GraphType, fieldType>& input) const {
				static constexpr Int RELOAD_ITERATIONS = (N_VERTICES > 2) ? (N_VERTICES / 3) : 1;
#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
				const auto labeling_search_start = std::chrono::steady_clock::now();
#endif

				vector<CanonBuilder2> attempts;
				vector<CanonBuilder2> next_attempts;
				vector<typename CanonBuilder2::CanonicalBucketAnalysis> analyses;
				vector<size_t> next_attempt_mask;
				GraphType& G = input.getValue();

				attempts.emplace_back();
				analyses.resize(1);
#if defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
				gc_standardizer_sort_profile::attempts_created.fetch_add(
					1,
					std::memory_order_relaxed
				);
#endif

				assign_type next_vertex_to_assign = attempts[0].init_starter_colors(G);
  
				assign_type unique_count;
				assign_type min_bucket_size;
				hash_int_type max_color;

				while (next_vertex_to_assign < N_VERTICES) {
					unique_count = 0;
					for (Int reload = 0; reload < RELOAD_ITERATIONS; ++reload) {
						analyses.resize(attempts.size());
						unique_count = 0;
						min_bucket_size = N_VERTICES + 1;
						max_color = 0;
						next_attempt_mask.clear();

						for (size_t i = 0; i < attempts.size(); ++i) {
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							const auto color_update_start = std::chrono::steady_clock::now();
#endif
							attempts[i].update_colors(G);
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							const auto color_update_stop = std::chrono::steady_clock::now();
							gc_standardizer_sort_profile::vertex_color_update_calls.fetch_add(
								1,
								std::memory_order_relaxed
							);
							gc_standardizer_sort_profile::vertex_color_update_nanoseconds.fetch_add(
								static_cast<unsigned long long>(
									std::chrono::duration_cast<std::chrono::nanoseconds>(
										color_update_stop - color_update_start
									).count()
								),
								std::memory_order_relaxed
							);

							const auto bucket_init_start = std::chrono::steady_clock::now();
#endif
							analyses[i] = attempts[i].assign_unique_colors_and_make_bucket(next_vertex_to_assign);
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							const auto bucket_init_stop = std::chrono::steady_clock::now();
							gc_standardizer_sort_profile::vertex_bucket_init_calls.fetch_add(
								1,
								std::memory_order_relaxed
							);
							gc_standardizer_sort_profile::vertex_bucket_init_nanoseconds.fetch_add(
								static_cast<unsigned long long>(
									std::chrono::duration_cast<std::chrono::nanoseconds>(
										bucket_init_stop - bucket_init_start
									).count()
								),
								std::memory_order_relaxed
							);
#endif

							if (analyses[i].unique_count > unique_count) {
								unique_count = analyses[i].unique_count;
								min_bucket_size = analyses[i].branch_bucket.size;
								max_color = analyses[i].prune_hash;
								next_attempt_mask.clear();
								next_attempt_mask.push_back(i);
							} else if (analyses[i].unique_count == unique_count) {
								if (analyses[i].branch_bucket.size < min_bucket_size) {
									min_bucket_size = analyses[i].branch_bucket.size;
									max_color = analyses[i].prune_hash;
									next_attempt_mask.clear();
									next_attempt_mask.push_back(i);
								} else if (analyses[i].branch_bucket.size == min_bucket_size) {
									if (max_color < analyses[i].prune_hash) {
										max_color = analyses[i].prune_hash;
										next_attempt_mask.clear();
										next_attempt_mask.push_back(i);
									} else if (max_color == analyses[i].prune_hash) {
										next_attempt_mask.push_back(i);
									}
								}
							}
						}

						
						if (next_attempt_mask.size() < attempts.size()) {
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							const auto mask_maintenance_start = std::chrono::steady_clock::now();
#endif
							size_t write = 0;
							for (size_t read : next_attempt_mask) {
								if (write != read) {
									attempts[write] = attempts[read];
									analyses[write] = analyses[read];
								}
								++write;
							}
							attempts.resize(next_attempt_mask.size());
							analyses.resize(next_attempt_mask.size());
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							const auto mask_maintenance_stop = std::chrono::steady_clock::now();
							gc_standardizer_sort_profile::mask_maintenance_calls.fetch_add(
								1,
								std::memory_order_relaxed
							);
							gc_standardizer_sort_profile::mask_maintenance_nanoseconds.fetch_add(
								static_cast<unsigned long long>(
									std::chrono::duration_cast<std::chrono::nanoseconds>(
										mask_maintenance_stop - mask_maintenance_start
									).count()
								),
								std::memory_order_relaxed
							);
#endif
						}

						if (unique_count > 0) {
							break;
						}
					}

						if (unique_count > 0) {
							next_vertex_to_assign += unique_count;
							continue;
						}

					if (min_bucket_size == 0) {
						break;
					}

					next_attempts.clear();
					next_attempts.reserve(next_attempt_mask.size()*min_bucket_size);

#if defined(GC_PROFILE_STANDARDIZER_SORT)
					const auto push_next_attempts_start = std::chrono::steady_clock::now();
#endif
					for (size_t i = 0; i < attempts.size(); ++i) {
						attempts[i].push_next_attempts(
							next_attempts,
							analyses[i].branch_bucket,
							next_vertex_to_assign
						);
					}
#if defined(GC_PROFILE_STANDARDIZER_SORT)
					const auto push_next_attempts_stop = std::chrono::steady_clock::now();
					gc_standardizer_sort_profile::push_next_attempts_calls.fetch_add(
						1,
						std::memory_order_relaxed
					);
					gc_standardizer_sort_profile::push_next_attempts_nanoseconds.fetch_add(
						static_cast<unsigned long long>(
							std::chrono::duration_cast<std::chrono::nanoseconds>(
								push_next_attempts_stop - push_next_attempts_start
							).count()
						),
						std::memory_order_relaxed
					);
#endif
#if defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
					gc_standardizer_sort_profile::attempts_created.fetch_add(
						next_attempts.size(),
						std::memory_order_relaxed
					);
#endif

					if (min_bucket_size == 2) {
						next_vertex_to_assign += 2;
					} else {
						++next_vertex_to_assign;
					}

					attempts.swap(next_attempts);
				}

#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
				const auto labeling_search_stop = std::chrono::steady_clock::now();
				gc_standardizer_sort_profile::labeling_search_calls.fetch_add(
					1,
					std::memory_order_relaxed
				);
				gc_standardizer_sort_profile::labeling_search_nanoseconds.fetch_add(
					static_cast<unsigned long long>(
						std::chrono::duration_cast<std::chrono::nanoseconds>(
							labeling_search_stop - labeling_search_start
						).count()
					),
					std::memory_order_relaxed
				);
#endif

				GraphType best_graph = input.getValue();
				bool have_best = false;
				bool containsPlus = false;
				bool containsMinus = false;
				bool checked_double_edge = false;

				for (const CanonBuilder2& attempt : attempts) {
					GraphType graph = input.getValue();
					signedInt sign = graph.permuteVertices(attempt.vertex_permutation);
					sign *= graph.directAndSortEdges();

					if constexpr (GraphType::SWAP_EDGE_SIGN == -1) {
						if (!checked_double_edge) {
							checked_double_edge = true;
							if (graph.has_double_edge()) {
								return BasisElement<GraphType, fieldType>(graph, 0);
							}
						}
					}

					const signedInt comparison = have_best ? graph.compare(best_graph) : -1;
					if (comparison < 0) {
						best_graph = graph;
						have_best = true;
						containsPlus = sign > 0;
						containsMinus = sign < 0;
					} else if (comparison == 0) {
						containsPlus = containsPlus || sign > 0;
						containsMinus = containsMinus || sign < 0;
					}

					if (containsPlus && containsMinus) {
						return BasisElement<GraphType, fieldType>(best_graph, 0);
					}
				}

				return BasisElement<GraphType, fieldType>(
					best_graph,
					containsPlus ? input.getCoefficient() : -input.getCoefficient()
				);
			}


			class CanonBuilder {
				public:
					GraphType G;
					Int n_assignedVertices;
					signedInt sign;

					CanonBuilder(const GraphType& initialGraph, Int n, signedInt s)
						: G(initialGraph), n_assignedVertices(n) {
							n_assignedVertices = n;

							sign = s*G.directAndSortEdges();
						}

					signedInt sort_edges() {
						sign *= G.directAndSortEdges();
						return sign;
					}

					Int vertex_value(Int v) const {
						if(v < n_assignedVertices) {
							return v;
						}
						return N_VERTICES;
					}

					array<array<Int, N_VERTICES+1>,N_VERTICES> score_vertices() {
						array<array<Int, N_VERTICES+1>,N_VERTICES> result{};
						for (Int i = G.N_HAIR; i < G.SIZE; i+=2) {
							Int a = G.half_edges[i];
							Int b = G.half_edges[i+1];
							++result[a][vertex_value(b)];
							++result[b][vertex_value(a)];
						}

						return result;
					}

					GraphType vertex_values_graph() const {
						GraphType fakeGraph = G;

						for (Int i = 0; i < GraphType::SIZE; ++i) {
							fakeGraph.half_edges[i] = vertex_value(G.half_edge(i));
						}

						return fakeGraph;
					}

					int compare(const CanonBuilder& other) {
						return combutils::compareHalfEdges(vertex_values_graph().half_edges, other.vertex_values_graph().half_edges);
					} 

					//we are assigning j <- n_assignedVertices
					CanonBuilder with_assigned_next(Int j) {
						GraphType copied_graph = G;
						signedInt s = copied_graph.swapVertices(j, n_assignedVertices);
						return CanonBuilder(copied_graph, n_assignedVertices +1, sign*s);
					}
			};

			class IsoCanonBuilder {
				public:
					GraphType G;
					Int n_assignedVertices;
					signedInt sign;
					IsomorphismType iso;

					IsoCanonBuilder(const GraphType& initialGraph, Int n, signedInt s, IsomorphismType input_iso = IsomorphismType{})
						: G(initialGraph), n_assignedVertices(n), sign(s), iso(input_iso) {
							n_assignedVertices = n;
							sign = s * direct_and_sort_edges_with_isomorphism(G, iso);
						}

					Int vertex_value(Int v) const {
						if(v < n_assignedVertices) {
							return v;
						}
						return N_VERTICES;
					}

					GraphType vertex_values_graph() const {
						GraphType fakeGraph = G;

						for (Int i = 0; i < GraphType::SIZE; ++i) {
							fakeGraph.half_edges[i] = vertex_value(G.half_edge(i));
						}

						return fakeGraph;
					}

					int compare(const IsoCanonBuilder& other) {
						return combutils::compareHalfEdges(vertex_values_graph().half_edges, other.vertex_values_graph().half_edges);
					}

					IsoCanonBuilder with_assigned_next(Int j) {
						GraphType copied_graph = G;
						IsomorphismType copied_iso = iso;
						signedInt s = copied_graph.swapVertices(j, n_assignedVertices);
						copied_iso = copied_iso.compose(vertex_swap_isomorphism(j, n_assignedVertices));
						return IsoCanonBuilder(copied_graph, n_assignedVertices + 1, sign * s, copied_iso);
					}

					static IsomorphismType vertex_swap_isomorphism(Int v, Int w) {
						IsomorphismType iso;
						std::swap(iso.vertex_permutation_data()[v], iso.vertex_permutation_data()[w]);
						iso.template compute_signs<N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>();
						return iso;
					}

					static IsomorphismType edge_flip_isomorphism(Int e) {
						IsomorphismType iso;
						iso.edge_flip_data()[e] = true;
						iso.template compute_signs<N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>();
						return iso;
					}

					static IsomorphismType edge_swap_isomorphism(Int i, Int j) {
						IsomorphismType iso;
						std::swap(iso.edge_permutation_data()[i], iso.edge_permutation_data()[j]);
						iso.template compute_signs<N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>();
						return iso;
					}

					static signedInt direct_and_sort_edges_with_isomorphism(GraphType& graph, IsomorphismType& iso) {
						signedInt overall_sign = 1;

						for (Int edge_index = 0; edge_index < N_EDGES; ++edge_index) {
							const Int base = GraphType::N_HAIR + 2 * edge_index;
							if (graph.half_edges[base] > graph.half_edges[base + 1]) {
								overall_sign *= graph.flipEdge(edge_index);
								iso = iso.compose(edge_flip_isomorphism(edge_index));
							}
						}

						for (Int i = 1; i < N_EDGES; ++i) {
							Int j = i;
							while (j > 0 && graph.compareEdge(j - 1, j) > 0) {
								overall_sign *= graph.swapEdges(j - 1, j);
								iso = iso.compose(edge_swap_isomorphism(j - 1, j));
								--j;
							}
						}

						return overall_sign;
					}
			};

			BasisElement<GraphType, fieldType> standardize(BasisElement<GraphType, fieldType>& input) {
				return standardize(input.getValue(), input.getCoefficient());
			}    

			bigInt automorphism_group_size(const GraphType& input_graph) const {
				return minimizing_isomorphisms(input_graph).size();
			}

			std::vector<IsomorphismType> minimizing_isomorphisms(const GraphType& input_graph) const {
				GraphType graph = input_graph;
				auto attempts = canonical_attempts_with_isomorphisms(graph);
				const auto& minimizers = attempts[N_VERTICES % 2];

				std::vector<IsomorphismType> result;
				result.reserve(minimizers.size());
				for (const auto& attempt : minimizers) {
					result.push_back(attempt.iso);
				}
				return result;
			}


			BasisElement<GraphType, fieldType> standardize(GraphType& graph, fieldType k) const {
				if (GraphType::SWAP_EDGE_SIGN == -1 ) {
					k *= graph.directAndSortEdges();
					if (graph.has_double_edge()) {
						return BasisElement<GraphType, fieldType>(graph, 0); 
					}
				}

				CanonBuilder G = assignHair(graph);
				vector<CanonBuilder> attempts[2];

				attempts[G.n_assignedVertices%2].push_back(G);

				for (Int n = G.n_assignedVertices; n < N_VERTICES; n++) {                    

					attempts[(n+1)%2].clear();
					attempts[(n+1)%2].reserve(attempts[n%2].size() * (N_VERTICES - n));
					for (CanonBuilder attempt : attempts[n%2]) {
						for (Int l = attempt.n_assignedVertices; l < N_VERTICES; ++l) {
							CanonBuilder next = attempt.with_assigned_next(l);
							if (attempts[(n+1)%2].empty()) {
								attempts[(n+1)%2].push_back(next);
								continue;
							}
							int comparison = next.compare(attempts[(n+1)%2].back());
							if (comparison  > 0 ) {
								continue;
							}
							if (comparison  < 0 ) {
								attempts[(n+1)%2].clear();
							}
							attempts[(n+1)%2].push_back(next);
						}

					}
				}

				//cout << "Aut size = " << attempts[N_VERTICES%2].size() << endl;

				bool containsPlus = false;
				bool containsMinus = false;

				for (bigInt i = 0; i < attempts[N_VERTICES%2].size(); i++) {
					if (attempts[N_VERTICES%2][i].sign > 0) {
						containsPlus = true;
					} else {
						containsMinus = true;
					};

					if (containsPlus && containsMinus) {

						return BasisElement<GraphType, fieldType>(attempts[N_VERTICES%2][0].G, 0);
					}
				}


				BasisElement<GraphType, fieldType> final_form;

				if (containsPlus) {
					final_form = BasisElement<GraphType, fieldType>(attempts[N_VERTICES%2][0].G, k);
				} else {
					final_form= BasisElement<GraphType, fieldType>(attempts[N_VERTICES%2][0].G, -k);
				}

				return final_form;
			}

			CanonBuilder assignHair(GraphType G) const {
				Int n_assigned = 0;
				signedInt sign = 1;
				for (Int i = 0; i < G.N_HAIR ; ++i) {
					if (G.half_edges[i] > n_assigned) {
						sign *= G.swapVertices(G.half_edges[i], n_assigned);
						n_assigned ++;
					} 
				}
				return CanonBuilder(G, n_assigned, sign);
			}

			IsoCanonBuilder assignHair_with_isomorphisms(GraphType G) const {
				Int n_assigned = 0;
				signedInt sign = 1;
				IsomorphismType iso;
				for (Int i = 0; i < G.N_HAIR ; ++i) {
					const Int hair_vertex = G.half_edges[i];
					if (hair_vertex > n_assigned) {
						sign *= G.swapVertices(hair_vertex, n_assigned);
						iso = iso.compose(IsoCanonBuilder::vertex_swap_isomorphism(hair_vertex, n_assigned));
						n_assigned ++;
					}
				}
				return IsoCanonBuilder(G, n_assigned, sign, iso);
			}

			array<vector<IsoCanonBuilder>, 2> canonical_attempts_with_isomorphisms(GraphType& graph) const {
				IsoCanonBuilder G = assignHair_with_isomorphisms(graph);
				array<vector<IsoCanonBuilder>, 2> attempts;

				attempts[G.n_assignedVertices % 2].push_back(G);

				for (Int n = G.n_assignedVertices; n < N_VERTICES; n++) {
					attempts[(n + 1) % 2].clear();
					for (IsoCanonBuilder attempt : attempts[n % 2]) {
						for (Int l = attempt.n_assignedVertices; l < N_VERTICES; ++l) {
							IsoCanonBuilder next = attempt.with_assigned_next(l);
							if (attempts[(n + 1) % 2].empty()) {
								attempts[(n + 1) % 2].push_back(next);
								continue;
							}
							int comparison = next.compare(attempts[(n + 1) % 2].back());
							if (comparison > 0) {
								continue;
							}
							if (comparison < 0) {
								attempts[(n + 1) % 2].clear();
							}
							attempts[(n + 1) % 2].push_back(next);
						}
					}
				}

				return attempts;
			}
	};
