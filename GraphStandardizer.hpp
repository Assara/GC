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
	inline std::atomic<unsigned long long> create_final_attempts_nanoseconds{0};
	inline std::atomic<unsigned long long> sort_and_filter_nanoseconds{0};

	inline void reset() {
		create_final_attempts_nanoseconds.store(0, std::memory_order_relaxed);
		sort_and_filter_nanoseconds.store(0, std::memory_order_relaxed);
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
					for (Int e = 0; e < G.half_edges.size(); ++e) {
						++colors[G.half_edges[e]];
					}
					if constexpr (GraphType::N_HAIR > 0) {
						for (Int i = 0; i < GraphType::N_HAIR; ++i) {
							colors[G.half_edges[i]] += hash(i);
						}
					}
					for (auto& color : colors) {
						color = hash(color);
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

			using FinalAttemptSet = std::pair<vector<CanonBuilder3>, vector<std::size_t>>;

			FinalAttemptSet create_final_attempts(const GraphType& G) const {
				static constexpr Int RELOAD_ITERATIONS = (N_VERTICES > 2) ? (N_VERTICES / 3) : 1;

				vector<CanonBuilder3> attempts;
				vector<CanonBuilder3> next_attempts;
				vector<std::size_t> valid_attempts;
				vector<std::size_t> next_valid_attempts;
				attempts.emplace_back(CanonBuilder3());
				attempts[0].init_starter_colors(G);
				valid_attempts.push_back(0);

				signedInt cmp;

				while (attempts[valid_attempts[0]].next_to_assign < N_VERTICES) {
					for (Int reload = 0; reload < RELOAD_ITERATIONS; ++reload) {
						next_valid_attempts.clear();
						next_valid_attempts.reserve(valid_attempts.size());

						if (valid_attempts.size() == 1) {
							CanonBuilder3& attempt = attempts[valid_attempts[0]];
							attempt.update_colors(G);
							attempt.update_vertex_groups();
							continue;
						}

						for (const std::size_t attempt_index : valid_attempts) {
							attempts[attempt_index].update_colors(G);
							attempts[attempt_index].update_vertex_groups();

							if (next_valid_attempts.empty()) {
								next_valid_attempts.push_back(attempt_index);
							} else {
								cmp = attempts[attempt_index].compare(attempts[next_valid_attempts.back()]);
								if (cmp > 0) {
									next_valid_attempts[0] = attempt_index;
									next_valid_attempts.resize(1);
								} else if (cmp == 0) {
									next_valid_attempts.push_back(attempt_index);
								}
							}
						}
						valid_attempts.swap(next_valid_attempts);
					}

					if (attempts[valid_attempts[0]].next_to_assign >= N_VERTICES) {
						break;
					}

					const auto branch_range = attempts[valid_attempts[0]].branching_group_range();
					if (branch_range.first >= branch_range.second) {
						for (const std::size_t attempt_index : valid_attempts) {
							attempts[attempt_index].assign_remaining_groups_in_order();
						}
						break;
					}

					next_attempts.clear();
					next_attempts.reserve(
						static_cast<std::size_t>(branch_range.second - branch_range.first) * valid_attempts.size()
					);

					for (const std::size_t attempt_index : valid_attempts) {
						attempts[attempt_index].branch(branch_range, next_attempts);
					}

					attempts.swap(next_attempts);
					valid_attempts.resize(attempts.size());
					for (std::size_t i = 0; i < attempts.size(); ++i) {
						valid_attempts[i] = i;
					}
				}

				return {std::move(attempts), std::move(valid_attempts)};
			}

			BasisElement<GraphType, fieldType> standardize3(const BasisElement<GraphType, fieldType>& input) const {
#if defined(GC_PROFILE_STANDARDIZER_SORT)
				const auto create_final_attempts_start = std::chrono::steady_clock::now();
#endif
				const GraphType& G = input.getValue();
				auto [attempts, valid_attempts] = create_final_attempts(G);

#if defined(GC_PROFILE_STANDARDIZER_SORT)
				const auto create_final_attempts_stop = std::chrono::steady_clock::now();
				gc_standardizer_sort_profile::create_final_attempts_nanoseconds.fetch_add(
					static_cast<unsigned long long>(
						std::chrono::duration_cast<std::chrono::nanoseconds>(
							create_final_attempts_stop - create_final_attempts_start
						).count()
					),
					std::memory_order_relaxed
				);
				const auto sort_and_filter_start = std::chrono::steady_clock::now();
				const auto record_sort_and_filter_stop = [&]() {
					const auto sort_and_filter_stop = std::chrono::steady_clock::now();
					gc_standardizer_sort_profile::sort_and_filter_nanoseconds.fetch_add(
						static_cast<unsigned long long>(
							std::chrono::duration_cast<std::chrono::nanoseconds>(
								sort_and_filter_stop - sort_and_filter_start
							).count()
						),
						std::memory_order_relaxed
					);
				};
#endif

				typename GraphType::Basis best_basis = GraphType::assignPermutedDirectedSortedEdgesBasis(
						input,
						attempts[valid_attempts[0]].vertex_permutation
				);

				if (valid_attempts.size() == 1 || best_basis.getCoefficient() == fieldType{0}) {
#if defined(GC_PROFILE_STANDARDIZER_SORT)
					record_sort_and_filter_stop();
#endif
					return best_basis;
				}

				typename GraphType::Basis attempt_basis;
				signedInt comparison;
			
				for (std::size_t i = 1; i < valid_attempts.size(); ++i) {
					const std::size_t attempt_index = valid_attempts[i];
					attempt_basis = GraphType::assignPermutedDirectedSortedEdgesBasis(
						input,
						attempts[attempt_index].vertex_permutation);
					
					comparison = best_basis.compare(attempt_basis);	
					if  (comparison == 0) {
						if (attempt_basis.getCoefficient() == -best_basis.getCoefficient()) {
							best_basis.set_coefficient(fieldType(0));
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							record_sort_and_filter_stop();
#endif
							return best_basis;
						}

					} else if (comparison < 0) {
						best_basis = attempt_basis;
					} 					
				}

#if defined(GC_PROFILE_STANDARDIZER_SORT)
				record_sort_and_filter_stop();
#endif
				return best_basis;
			}

			class LexicographicalBuilder {
				public:
					GraphType G;
					Int n_assignedVertices;
					signedInt sign;

					LexicographicalBuilder(const GraphType& initialGraph, Int n, signedInt s)
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

					int compare(const LexicographicalBuilder& other) {
						return combutils::compareHalfEdges(vertex_values_graph().half_edges, other.vertex_values_graph().half_edges);
					} 

					//we are assigning j <- n_assignedVertices
					LexicographicalBuilder with_assigned_next(Int j) {
						GraphType copied_graph = G;
						signedInt s = copied_graph.swapVertices(j, n_assignedVertices);
						return LexicographicalBuilder(copied_graph, n_assignedVertices +1, sign*s);
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

			BasisElement<GraphType, fieldType> lexicographical_standardize(BasisElement<GraphType, fieldType>& input) {
				return lexicographical_standardize(input.getValue(), input.getCoefficient());
			}    

			bigInt automorphism_group_size(const GraphType& input_graph) const {
				return minimizing_isomorphisms(input_graph).size();
			}

			std::vector<IsomorphismType> minimizing_isomorphisms(const GraphType& input_graph) const {
				const typename GraphType::Basis input(input_graph, fieldType{1});
				auto [attempts, valid_attempts] = create_final_attempts(input_graph);

				std::vector<IsomorphismType> minimizers;
				minimizers.reserve(valid_attempts.size());
				typename GraphType::Basis best_basis;
				bool have_best = false;

				for (const std::size_t attempt_index : valid_attempts) {
					typename GraphType::ThisGraph graph;
					IsomorphismType iso;
					const signedInt sign = graph.assignPermutedDirectedSortedEdgesWithIsomorphism(
						input_graph,
						attempts[attempt_index].vertex_permutation,
						iso
					);
					if (sign == 0) {
						continue;
					}

					typename GraphType::Basis attempt_basis(std::move(graph), fieldType(sign) * input.getCoefficient());
					if (!have_best) {
						best_basis = attempt_basis;
						minimizers.clear();
						minimizers.push_back(iso);
						have_best = true;
						continue;
					}

					const signedInt comparison = best_basis.compare(attempt_basis);
					if (comparison < 0) {
						best_basis = attempt_basis;
						minimizers.clear();
						minimizers.push_back(iso);
					} else if (comparison == 0) {
						minimizers.push_back(iso);
					}
				}

				return minimizers;
			}


			BasisElement<GraphType, fieldType> lexicographical_standardize(GraphType& graph, fieldType k) const {
				if (GraphType::SWAP_EDGE_SIGN == -1 ) {
					k *= graph.directAndSortEdges();
					if (graph.has_double_edge()) {
						return BasisElement<GraphType, fieldType>(graph, 0); 
					}
				}

				LexicographicalBuilder G = assignHair(graph);
				vector<LexicographicalBuilder> attempts[2];

				attempts[G.n_assignedVertices%2].push_back(G);

				for (Int n = G.n_assignedVertices; n < N_VERTICES; n++) {                    

					attempts[(n+1)%2].clear();
					attempts[(n+1)%2].reserve(attempts[n%2].size() * (N_VERTICES - n));
					for (LexicographicalBuilder attempt : attempts[n%2]) {
						for (Int l = attempt.n_assignedVertices; l < N_VERTICES; ++l) {
							LexicographicalBuilder next = attempt.with_assigned_next(l);
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

			LexicographicalBuilder assignHair(GraphType G) const {
				Int n_assigned = 0;
				signedInt sign = 1;
				for (Int i = 0; i < G.N_HAIR ; ++i) {
					if (G.half_edges[i] > n_assigned) {
						sign *= G.swapVertices(G.half_edges[i], n_assigned);
						n_assigned ++;
					} 
				}
				return LexicographicalBuilder(G, n_assigned, sign);
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
