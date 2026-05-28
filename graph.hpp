#pragma once

#include <array>
#include <vector>
#include <memory>
#include <utility>
#include <cstring>
#include <unordered_set>
#include <iostream>
#include <algorithm>
#if defined(GC_PROFILE_STANDARDIZER_SORT)
#include <atomic>
#include <chrono>
#endif
#include "types.hpp"
#include "permutation.hpp"
#include "CombinatorialUtils.hpp"
#include "VectorSpace/LinComb.hpp"
#include "graph_hash.hpp"
#include "GraphStandardizer.hpp"
#include <ranges>

using namespace std;

template<typename T, typename k>
class BasisElement;

template <
Int N_VERTICES,
	Int N_EDGES,
	Int N_OUT_HAIR,
	Int N_IN_HAIR,
	signedInt c,
	signedInt d,
	typename fieldType
	>
	class Graph {
		public:
			static constexpr Int N_VERTICES_ = N_VERTICES;
			static constexpr Int SIZE = N_IN_HAIR + N_OUT_HAIR + (Int)2 * N_EDGES;
			static constexpr Int N_HAIR = N_IN_HAIR + N_OUT_HAIR;
			static constexpr Int N_OUT_HAIR_ = N_OUT_HAIR;
			static constexpr Int N_IN_HAIR_ = N_IN_HAIR;
			static constexpr Int N_EDGES_ = N_EDGES;
			static constexpr signedInt C_ = c;
			static constexpr signedInt D_ = d;
			static constexpr signedInt FLIP_EDGE_SIGN = ((c % 2) != 0 && (d % 2) != 0) ? -1 : 1;
			static constexpr signedInt SWAP_EDGE_SIGN = (((c + d) % 2) != 0) ? -1 : 1;
			static constexpr signedInt SWAP_VERTICES_SIGN = -1 * SWAP_EDGE_SIGN;
			using Field = fieldType;

			using SplitGraph = Graph<N_VERTICES + 1, N_EDGES + 1, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>;

			using ContGraph = Graph<N_VERTICES - 1, N_EDGES - 1, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>;

			using ThisGraph = Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>;

			using ExtraEdgeGraph = Graph<N_VERTICES, N_EDGES + 1, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>;

			using Basis = BasisElement<ThisGraph, fieldType>;
			
			template <typename NewFieldType>
			using RebindField = Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, NewFieldType>;

			template <signedInt NewC, signedInt NewD>
			using RebindDegree = Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, NewC, NewD, fieldType>;


			Graph() = default;

			explicit Graph(const array<Int, SIZE>& arr) : half_edges(arr) {}

			inline signedInt custom_filter() const {
				return n_odd_pairs();
			}

			inline signedInt n_odd_pairs() const {
				return std::ranges::count_if(valence_array(), [](Int v){ return v % 2 != 0; }) / 2;
			}

			inline signedInt n_4valent_vertices() const {
				return std::ranges::count_if(valence_array(), [](Int v){ return v == 4; });
			}

			inline Int find_edge_index(Int u, Int v) const {
				for (Int i = 0; i < N_EDGES_; ++i) {
					auto [a, b] = getEdge(i);
					if ((a == u && b == v) || (a == v && b == u)) {
						return i;
					}
				}
				return N_EDGES_;
			}



			inline signedInt total_order_comp(const ThisGraph& other) const {
				return combutils::compareHalfEdges(half_edges, other.half_edges);
			}

			signedInt compare(const ThisGraph& other) const {

				// for homological filtration
				/*
				   signedInt r = other.custom_filter() - custom_filter();
				   if (r != 0) {
				   return r;
				   }
				 */

				// For absolute ordering
				return total_order_comp(other);
			}

			bool operator<(const Graph& o) const { return compare(o) < 0; }

			static void std(Basis& b) {
				GraphStandardizer<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType> s;
				b = s.standardize3(b);
			}

			static Basis canonized(Basis& b) {
				GraphStandardizer<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType> s;
				return s.standardize3(b);
			}

			ThisGraph canonical_represesentation() {
				Basis be = Basis(*this);
				return ThisGraph::canonized(be).getValue();
			}

			template <typename NewFieldType>
			RebindField<NewFieldType> cast_field() const {
				RebindField<NewFieldType> result;
				result.half_edges = half_edges;
				return result;
			}


			pair<Int, Int> getEdge(Int i) const {
				return { half_edges[N_HAIR + 2 * i], half_edges[N_HAIR + 2 * i + 1] };
			}

			void setEdge(Int i, Int v, Int w) {
				half_edges[N_HAIR + 2 * i] = v;
				half_edges[N_HAIR + 2 * i + 1] = w;
			}

			vector<Int> adjacent(Int v) const {
				vector<Int> result;
				result.reserve(SIZE);
				for (Int i = 0; i < SIZE; ++i) {
					if (half_edges[i] == v) {
						result.push_back(i);
					}
				}
				return result;
			}

			array<Int, N_VERTICES> valence_array() const {
				array<Int, N_VERTICES> valence{};
				for (Int v : half_edges) {
					++valence[v];
				}
				return valence;
			}

			bool is_connected_after_removing_vertex(Int removed_vertex) const {
				if constexpr (N_VERTICES <= 2) {
					return true;
				}

				if (removed_vertex >= N_VERTICES) {
					return false;
				}

				Int start = 0;
				while (start == removed_vertex) {
					++start;
				}

				array<bool, N_VERTICES> visited{};
				vector<Int> stack;
				stack.reserve(N_VERTICES);
				visited[removed_vertex] = true;
				visited[start] = true;
				stack.push_back(start);

				while (!stack.empty()) {
					const Int current = stack.back();
					stack.pop_back();

					for (Int e = 0; e < N_EDGES; ++e) {
						auto [u, v] = getEdge(e);
						if (u == removed_vertex || v == removed_vertex) {
							continue;
						}

						Int next = N_VERTICES;
						if (u == current && !visited[v]) {
							next = v;
						} else if (v == current && !visited[u]) {
							next = u;
						}

						if (next != N_VERTICES) {
							visited[next] = true;
							stack.push_back(next);
						}
					}
				}

				for (bool seen : visited) {
					if (!seen) {
						return false;
					}
				}
				return true;
			}

			bigInt count_triangles() const {
				array<array<bool, N_VERTICES>, N_VERTICES> adjacency{};

				for (Int e = 0; e < N_EDGES; ++e) {
					auto [u, v] = getEdge(e);
					adjacency[u][v] = true;
					adjacency[v][u] = true;
				}

				bigInt count = 0;
				for (Int u = 0; u < N_VERTICES; ++u) {
					for (Int v = u + 1; v < N_VERTICES; ++v) {
						if (!adjacency[u][v]) {
							continue;
						}
						for (Int w = v + 1; w < N_VERTICES; ++w) {
							if (adjacency[u][w] && adjacency[v][w]) {
								++count;
							}
						}
					}
				}

				return count;
			}

			vector<vector<Int>> triangle_banks_at_vertex(Int center) const {
				array<array<bool, N_VERTICES>, N_VERTICES> adjacency{};
				for (Int e = 0; e < N_EDGES; ++e) {
					auto [u, v] = getEdge(e);
					adjacency[u][v] = true;
					adjacency[v][u] = true;
				}

				vector<Int> neighbors;
				neighbors.reserve(N_VERTICES);
				for (Int v = 0; v < N_VERTICES; ++v) {
					if (v != center && adjacency[center][v]) {
						neighbors.push_back(v);
					}
				}

				vector<vector<Int>> banks;
				vector<bool> visited(neighbors.size(), false);
				vector<bool> has_triangle(neighbors.size(), false);

				for (size_t i = 0; i < neighbors.size(); ++i) {
					for (size_t j = 0; j < neighbors.size(); ++j) {
						if (i != j && adjacency[neighbors[i]][neighbors[j]]) {
							has_triangle[i] = true;
							break;
						}
					}
				}

				for (size_t start = 0; start < neighbors.size(); ++start) {
					if (visited[start] || !has_triangle[start]) {
						continue;
					}

					vector<Int> bank;
					vector<size_t> stack{start};
					visited[start] = true;

					while (!stack.empty()) {
						const size_t i = stack.back();
						stack.pop_back();
						const Int u = neighbors[i];
						bank.push_back(u);

						for (size_t j = 0; j < neighbors.size(); ++j) {
							if (visited[j] || !has_triangle[j]) {
								continue;
							}
							if (adjacency[u][neighbors[j]]) {
								visited[j] = true;
								stack.push_back(j);
							}
						}
					}

					banks.push_back(std::move(bank));
				}

				return banks;
			}

			vector<vector<Int>> triangle_banks_at_zero() const {
				return triangle_banks_at_vertex(0);
			}

			VectorSpace::LinComb<SplitGraph, fieldType> split_vertex_differential(fieldType coef) const {
				VectorSpace::LinComb<SplitGraph, fieldType> splits = unsorted_splits(coef);
				splits.standardize_and_sort();
				return splits;
			}

			VectorSpace::LinComb<SplitGraph, fieldType> split_vertex_differential_4valent(fieldType coef) const {
				VectorSpace::LinComb<SplitGraph, fieldType> splits = unsorted_4valent_splits(coef);
				splits.standardize_and_sort();
				return splits;
			}

			VectorSpace::LinComb<SplitGraph, fieldType> split_vertex_differential_4valent_preserving_bank_at_zero(fieldType coef) const {
				VectorSpace::LinComb<SplitGraph, fieldType> splits = unsorted_4valent_splits_preserving_bank_at_zero(coef);
				splits.standardize_and_sort();
				return splits;
			}

			bool staying_neighbors_form_single_bank_at_zero(const vector<Int>& staying_neighbors) const {
				if (staying_neighbors.empty()) {
					return false;
				}

				array<array<bool, N_VERTICES>, N_VERTICES> adjacency{};
				for (Int e = 0; e < N_EDGES; ++e) {
					auto [u, v] = getEdge(e);
					adjacency[u][v] = true;
					adjacency[v][u] = true;
				}

				vector<bool> in_bank(staying_neighbors.size(), false);
				for (size_t i = 0; i < staying_neighbors.size(); ++i) {
					for (size_t j = 0; j < staying_neighbors.size(); ++j) {
						if (i != j && adjacency[staying_neighbors[i]][staying_neighbors[j]]) {
							in_bank[i] = true;
							break;
						}
					}
				}

				size_t bank_vertex_count = 0;
				size_t extra_count = 0;
				for (bool flag : in_bank) {
					if (flag) {
						++bank_vertex_count;
					} else {
						++extra_count;
					}
				}

				if (bank_vertex_count < 2 || extra_count > 1) {
					return false;
				}

				vector<bool> visited(staying_neighbors.size(), false);
				size_t n_components = 0;
				for (size_t start = 0; start < staying_neighbors.size(); ++start) {
					if (visited[start] || !in_bank[start]) {
						continue;
					}
					++n_components;
					vector<size_t> stack{start};
					visited[start] = true;
					while (!stack.empty()) {
						const size_t i = stack.back();
						stack.pop_back();
						for (size_t j = 0; j < staying_neighbors.size(); ++j) {
							if (visited[j] || !in_bank[j]) {
								continue;
							}
							if (adjacency[staying_neighbors[i]][staying_neighbors[j]]) {
								visited[j] = true;
								stack.push_back(j);
							}
						}
					}
				}

				return n_components == 1;
			}



			VectorSpace::LinComb<SplitGraph, fieldType> unsorted_splits(fieldType coef) const {
				vector<vector<Int>> adjRepresentation;
				adjRepresentation.reserve(N_VERTICES);

				bigInt resultSize = 0;
				for (Int v = 0; v < N_VERTICES; v++) {
					adjRepresentation.push_back(adjacent(v));
					resultSize += combutils::n_splits(adjRepresentation.back().size());
				}

				VectorSpace::LinComb<SplitGraph, fieldType> result;
				result.reserve(resultSize);

				for (Int v = 0; v < N_VERTICES; v++) {
					split_vertex(v, adjRepresentation[v], result, coef);
				}
				return result;
			}

			VectorSpace::LinComb<SplitGraph, fieldType> unsorted_full_splits(fieldType coef) const {
				vector<vector<Int>> adjRepresentation;
				adjRepresentation.reserve(N_VERTICES);

				VectorSpace::LinComb<SplitGraph, fieldType> result;

				for (Int v = 0; v < N_VERTICES; v++) {
					adjRepresentation.push_back(adjacent(v));
				}

				for (Int v = 0; v < N_VERTICES; v++) {
					split_vertex_full(v, adjRepresentation[v], result, coef);
				}
				return result;
			}

			VectorSpace::LinComb<SplitGraph, fieldType> unsorted_4valent_splits(fieldType coef) const {
				vector<vector<Int>> adjRepresentation;
				adjRepresentation.reserve(N_VERTICES);

				bigInt resultSize = 0;
				for (Int v = 0; v < N_VERTICES; v++) {
					adjRepresentation.push_back(adjacent(v));
					if (adjRepresentation.back().size() > 4) {
						const bigInt valence = adjRepresentation.back().size();
						resultSize += (valence * (valence - 1) * (valence - 2)) / 6;
					}
				}

				VectorSpace::LinComb<SplitGraph, fieldType> result;
				result.reserve(resultSize);

				for (Int v = 0; v < N_VERTICES; v++) {
					if (adjRepresentation[v].size() <= 4) {
						continue;
					}

					Int max_index = static_cast<Int>(adjRepresentation[v].size()) - 1;
					vector<Int> moved_half_edges = combutils::firstSubset(0, 3);
					do {
						result.append_in_basis_order(splitGraph(v, adjRepresentation[v], moved_half_edges), coef);
					} while (combutils::nextSubset(moved_half_edges, max_index));
				}

				return result;
			}

			VectorSpace::LinComb<SplitGraph, fieldType> unsorted_4valent_splits_preserving_bank_at_zero(fieldType coef) const {
				VectorSpace::LinComb<SplitGraph, fieldType> result;

				const Int split_vertex = 0;
				const vector<Int> adjacent_zero = adjacent(split_vertex);
				if (adjacent_zero.size() <= 4) {
					return result;
				}

				const auto banks = triangle_banks_at_zero();
				if (banks.empty()) {
					return result;
				}

				const auto opposite_vertex = [&](Int half_edge_index) -> Int {
					const Int mate =
						((half_edge_index - N_HAIR) % 2 == 0) ? half_edge_index + 1 : half_edge_index - 1;
					return half_edges[mate];
				};

				std::unordered_map<Int, Int> neighbor_to_adjacent_index;
				neighbor_to_adjacent_index.reserve(adjacent_zero.size());
				for (Int i = 0; i < static_cast<Int>(adjacent_zero.size()); ++i) {
					neighbor_to_adjacent_index.emplace(opposite_vertex(adjacent_zero[i]), i);
				}

				if (banks.size() == 1) {
					const auto& bank = banks.front();
					vector<Int> outsiders;
					for (const auto& [neighbor, idx] : neighbor_to_adjacent_index) {
						(void)idx;
						if (std::find(bank.begin(), bank.end(), neighbor) == bank.end()) {
							outsiders.push_back(neighbor);
						}
					}

					if (outsiders.size() == 1) {
						array<array<bool, N_VERTICES>, N_VERTICES> adjacency{};
						for (Int e = 0; e < N_EDGES; ++e) {
							auto [u, v] = getEdge(e);
							adjacency[u][v] = true;
							adjacency[v][u] = true;
						}

						std::unordered_map<Int, vector<Int>> bank_neighbors;
						vector<Int> endpoints;
						for (Int u : bank) {
							auto& nbrs = bank_neighbors[u];
							for (Int v : bank) {
								if (u != v && adjacency[u][v]) {
									nbrs.push_back(v);
								}
							}
							if (nbrs.size() == 1) {
								endpoints.push_back(u);
							}
						}

						if (endpoints.size() == 2) {
							const Int outsider = outsiders.front();
							for (Int endpoint : endpoints) {
								const Int next_in_bank = bank_neighbors[endpoint].front();
								vector<Int> moved_half_edges{
									neighbor_to_adjacent_index[outsider],
									neighbor_to_adjacent_index[endpoint],
									neighbor_to_adjacent_index[next_in_bank]
								};
								std::sort(moved_half_edges.begin(), moved_half_edges.end());
								result.append_in_basis_order(splitGraph(split_vertex, adjacent_zero, moved_half_edges), coef);
							}
							return result;
						}
					}

					if (outsiders.empty() && static_cast<Int>(bank.size()) == static_cast<Int>(adjacent_zero.size())) {
						vector<Int> moved_left{0, 1, 2};
						vector<Int> moved_right;
						moved_right.push_back(0);
						moved_right.push_back(static_cast<Int>(adjacent_zero.size() - 1));
						moved_right.push_back(static_cast<Int>(adjacent_zero.size() - 2));
						std::sort(moved_left.begin(), moved_left.end());
						std::sort(moved_right.begin(), moved_right.end());
						result.append_in_basis_order(splitGraph(split_vertex, adjacent_zero, moved_left), coef);
						if (moved_right != moved_left) {
							result.append_in_basis_order(splitGraph(split_vertex, adjacent_zero, moved_right), coef);
						}
						return result;
					}
				}

				return result;
			}


			void split_vertex(Int split_vertex, const vector<Int>& adjacent, VectorSpace::LinComb<SplitGraph, fieldType>& result, fieldType coef) const {
				if (adjacent.size() < 4) {

					if (adjacent.size() == 2 ) {

						result.append_in_basis_order(splitGraph(split_vertex, adjacent, vector<Int>(adjacent.back())), -coef);
					}
					else if (adjacent.size() < 2) {
						result.append_in_basis_order(splitGraph(split_vertex, adjacent, vector<Int>()), coef);
					}
					//do nothing for adjacent.size() == 3
					return; 
				}

				Int max_index = adjacent.size() - 1;
				for (Int i = 2; i < max_index; i++) {
					vector<Int> S = combutils::firstSubset(1, i);

					do {
						result.append_in_basis_order(splitGraph(split_vertex, adjacent, S), coef);
					} while (combutils::nextSubset(S, max_index));
				}
			}

			void split_vertex_full(Int split_vertex, const vector<Int>& adjacent, VectorSpace::LinComb<SplitGraph, fieldType>& result, fieldType coef) const {
				if (adjacent.empty()) {
					result.append_in_basis_order(splitGraph(split_vertex, adjacent, vector<Int>()), coef);
					return;
				}

				if (adjacent.size() == 1) {
					result.append_in_basis_order(splitGraph(split_vertex, adjacent, vector<Int>()), coef);
					result.append_in_basis_order(splitGraph(split_vertex, adjacent, vector<Int>{0}), coef);
					return;
				}

				if (adjacent.size() == 2) {
					result.append_in_basis_order(splitGraph(split_vertex, adjacent, vector<Int>{0}), -coef);
					result.append_in_basis_order(splitGraph(split_vertex, adjacent, vector<Int>{1}), -coef);
					return;
				}

				if (adjacent.size() == 3) {
					return;
				}

				const Int max_index = adjacent.size() - 1;
				for (Int i = 2; i < max_index; i++) {
					vector<Int> S = combutils::firstSubset(0, i);

					do {
						result.append_in_basis_order(splitGraph(split_vertex, adjacent, S), coef);
					} while (combutils::nextSubset(S, max_index));
				}
			}

			SplitGraph splitGraph(Int split_vertex, const vector<Int>& adjacent, const vector<Int>& S) const {
				auto sg = SplitGraph();

				std::copy_n(half_edges.begin(), SIZE, sg.half_edges.begin());
				for (auto s : S) {
					sg.half_edges[adjacent[s]] = N_VERTICES;
				}

				sg.half_edges[SplitGraph::SIZE - 2] = split_vertex;
				sg.half_edges[SplitGraph::SIZE - 1] = N_VERTICES;

				return sg;
			}

			void add_splits_to_set(unordered_set<SplitGraph>& result) {
				VectorSpace::LinComb<SplitGraph, fieldType> splits = unsorted_splits();
				splits.standardize_all();


				constexpr fieldType zero{}; 
				for (auto const& be : splits) {
					if (be.getCoefficient() != zero) {
						result.emplace(be.getValue());
					}
				}
			}




			VectorSpace::LinComb<SplitGraph, fieldType> unsorted_even_splits(fieldType coef) const {
				vector<vector<Int>> adjRepresentation;
				adjRepresentation.reserve(N_VERTICES);

				bigInt resultSize = 0;
				for (Int v = 0; v < N_VERTICES; v++) {
					adjRepresentation.push_back(adjacent(v));
					resultSize += combutils::n_splits(adjRepresentation.back().size());
				}

				VectorSpace::LinComb<SplitGraph, fieldType> result;
				result.reserve(resultSize);

				for (Int v = 0; v < N_VERTICES; v++) {
					split_vertex_even(v, adjRepresentation[v], result, coef);
				}
				return result;
			}


			VectorSpace::LinComb<SplitGraph, fieldType> split_vertex_differential_even(fieldType coef) const {
				VectorSpace::LinComb<SplitGraph, fieldType> splits = unsorted_even_splits(coef);
				splits.standardize_and_sort();
				return splits;
			}




			void split_vertex_even(Int vertex_to_split, const vector<Int>& adjacent, VectorSpace::LinComb<SplitGraph,fieldType>& result, fieldType coef) const {
				//odd vertices and bivalent vertices are handled the same way

				if (adjacent.size()%2 == 1 || adjacent.size() == 2) {
					split_vertex(vertex_to_split, adjacent, result, coef);
					return;
				}

				if (adjacent.size() < 5) {
					// do nothing for 4
					return; 
				}

				Int max_index = adjacent.size() - 1;
				for (Int i = 2; i < max_index; i++) {
					vector<Int> S = combutils::firstSubset(1, i);

					do {
						result.append_in_basis_order(splitGraph(vertex_to_split, adjacent, S), coef);
					} while (combutils::nextSubset(S, max_index));
				}
			}

			void add_even_splits_to_set(unordered_set<SplitGraph>& result) const {
				VectorSpace::LinComb<SplitGraph, fieldType> splits = unsorted_even_splits(fieldType{1});
				splits.standardize_all();

				constexpr fieldType zero{}; 
				for (auto const& be : splits) {
					if (be.getCoefficient() != zero) {
						result.emplace(be.getValue());
					}
				}
			}


			vector<SplitGraph> even_splits_vector() const {
				VectorSpace::LinComb<SplitGraph, fieldType> splits = unsorted_even_splits(fieldType{1});
				splits.standardize_all();

				splits.sort_without_deduplicate();
				vector<SplitGraph> result;
				result.reserve(splits.size());

				for (const auto& be : splits) {
					if (be.getCoefficient() == fieldType{0}) continue;
					if (!result.empty() && be.getValue() == result.back()) continue;
					result.push_back(be.getValue);

				}

				return result;
			}


			VectorSpace::LinComb<ExtraEdgeGraph, fieldType> add_edge_differential(fieldType coef) const {
				VectorSpace::LinComb<ExtraEdgeGraph, fieldType> result;
				result.reserve(N_VERTICES * (N_VERTICES - 1) - N_EDGES);


				ExtraEdgeGraph base_graph;
				std::copy_n(half_edges.begin(), SIZE, base_graph.half_edges.begin() + 2);

				for (Int u = 0; u < N_VERTICES - 1; u++) {
					for (Int v = u+1; v < N_VERTICES; v++) {
						base_graph.half_edges[0] = u;
						base_graph.half_edges[1] = v;
						result.append_in_basis_order(base_graph, coef);
					}
				}
				result.standardize_and_sort();

				return result;
			}


			//only suitable after sorting and directing edges
			bool has_double_edge() {
				for (Int i = 0; i < ThisGraph::N_EDGES_ - 1; i++) {
					if (getEdge(i) == getEdge(i+1)) {
						return true;
					}
				}
				return false;

			}

			VectorSpace::LinComb<ContGraph, fieldType> contraction_differential(fieldType k) const {

				VectorSpace::LinComb<ContGraph, fieldType> result;

				for(Int i = 0; i < N_EDGES_; i++) {
					result.append_in_basis_order(contract_edge(i, k));
				}
				result.standardize_and_sort();
				return result;
			}

			unordered_set<ContGraph> contraction_set() const {

				unordered_set<ContGraph> result;
				result.reserve(N_EDGES_);

				for(Int i = 0; i < N_EDGES_; i++) {
					auto be = contract_edge(i, fieldType{1});
					ContGraph::std(be);

					if (be.getCoefficient() != fieldType{0}) {
						result.emplace(be.getValue());
					}		 
				}
				return result;
			}



			BasisElement<ContGraph, fieldType> contract_edge(Int i, fieldType k) const {
				const Int edge_index = N_HAIR + 2 * i;
				const Int contraction_vertex = min(half_edges[edge_index], half_edges[edge_index + 1]);
				const Int deletion_vertex = max(half_edges[edge_index], half_edges[edge_index + 1]);

				// skip for tadpoles
				if (contraction_vertex == deletion_vertex) {
					return BasisElement<ContGraph, fieldType>(ContGraph{}, static_cast<fieldType>(0));
				}

				BasisElement<ContGraph, fieldType> contracted(ContGraph(), k);

				for (Int j = 0; j < edge_index; ++j) {
					contracted.getValue().half_edges[j] =  
						contraction_value(half_edges[j], contraction_vertex, deletion_vertex);
				} 

				for (Int j = edge_index + 2 ; j < ThisGraph::SIZE; ++j) {
					contracted.getValue().half_edges[j - 2] =  // shift down by 2 to account for removed edge
						contraction_value(half_edges[j], contraction_vertex, deletion_vertex);
				}


				if (!contracted.getValue().is_connected_after_removing_vertex(contraction_vertex)) {
					contracted.multiplyCoefficient(fieldType{0});
					return contracted;
				}

				if ( (N_EDGES - i) % 2 == 0) {
					contracted.multiplyCoefficient(SWAP_EDGE_SIGN);
				} 

				if (deletion_vertex != N_VERTICES - 1) {
					contracted.multiplyCoefficient(SWAP_VERTICES_SIGN);
				} 

				if (half_edges[edge_index] > half_edges[edge_index + 1]) {
					contracted.multiplyCoefficient(FLIP_EDGE_SIGN);
				}


				return contracted;
			}

			static BasisElement<ContGraph, fieldType> contract_edge(const BasisElement<ThisGraph, fieldType>& be, Int i) {
				return be.getValue().contract_edge(i, be.getCoefficient());
			}

			static Int contraction_value(Int v, Int contraction_vertex, Int deletion_vertex) {
				if (v == deletion_vertex) {
					return contraction_vertex;
				}
				if (v == N_VERTICES - 1) {
					return deletion_vertex;
				}
				return v;

			}

			signedInt flipEdge(Int i) {
				const Int base = N_HAIR + 2 * i;
				std::swap(half_edges[base], half_edges[base + 1]);
				return FLIP_EDGE_SIGN;
			}

			signedInt swapEdges(Int i, Int j) {
				const Int base_i = N_HAIR + 2 * i;
				const Int base_j = N_HAIR + 2 * j;
				std::swap(half_edges[base_i], half_edges[base_j]);
				std::swap(half_edges[base_i + 1], half_edges[base_j + 1]);
				return SWAP_EDGE_SIGN;
			}

			signedInt directEdges() {
				signedInt sign = 1;
				for (Int edgeIndex = 0; edgeIndex < N_EDGES; ++edgeIndex) {
					Int base = N_HAIR + 2 * edgeIndex;
					if (half_edges[base] > half_edges[base + 1]) {
						sign *= flipEdge(edgeIndex);
					}
				}
				return sign;
			}

			signedInt compareEdge(Int e1, Int e2) const {
				Int base1 = N_HAIR + 2 * e1;
				Int base2 = N_HAIR + 2 * e2;
				return std::memcmp(&half_edges[base1], &half_edges[base2], 2 * sizeof(Int));
			}

			signedInt compareEdgeToPair(int edge, Int pivot0, Int pivot1) const {
				const Int base = N_HAIR + 2 * static_cast<Int>(edge);
				if (half_edges[base] != pivot0) {
					return static_cast<signedInt>(half_edges[base]) - static_cast<signedInt>(pivot0);
				}
				return static_cast<signedInt>(half_edges[base + 1]) - static_cast<signedInt>(pivot1);
			}

			Int medianEdgeIndex(Int a, Int b, Int c_edge) const {
				const bool ab = compareEdge(a, b) < 0;
				const bool ac = compareEdge(a, c_edge) < 0;
				const bool bc = compareEdge(b, c_edge) < 0;
				if (ab) {
					return bc ? b : (ac ? c_edge : a);
				}
				return bc ? (ac ? a : c_edge) : b;
			}

			signedInt sortEdgesQuick() {
#if defined(GC_PROFILE_STANDARDIZER_SORT)
				const auto sort_start = std::chrono::steady_clock::now();
				unsigned long long swap_count = 0;
#endif
				signedInt overallSign = 1;

				auto insertion_sort_range = [&](int left, int right) {
					for (int i = left + 1; i <= right; ++i) {
						int j = i;
						while (j > left && compareEdge(static_cast<Int>(j - 1), static_cast<Int>(j)) > 0) {
							overallSign *= swapEdges(static_cast<Int>(j - 1), static_cast<Int>(j));
#if defined(GC_PROFILE_STANDARDIZER_SORT)
							++swap_count;
#endif
							--j;
						}
					}
				};

				auto quick_sort = [&](auto&& self, int left, int right) -> void {
					while (right - left > 8) {
						const Int pivot_index = medianEdgeIndex(
							static_cast<Int>(left),
							static_cast<Int>(left + (right - left) / 2),
							static_cast<Int>(right)
						);
						const Int pivot_base = N_HAIR + 2 * pivot_index;
						const Int pivot0 = half_edges[pivot_base];
						const Int pivot1 = half_edges[pivot_base + 1];

						int i = left;
						int j = right;
						while (i <= j) {
							while (compareEdgeToPair(i, pivot0, pivot1) < 0) {
								++i;
							}
							while (compareEdgeToPair(j, pivot0, pivot1) > 0) {
								--j;
							}
							if (i <= j) {
								if (i < j && compareEdge(static_cast<Int>(i), static_cast<Int>(j)) != 0) {
									overallSign *= swapEdges(static_cast<Int>(i), static_cast<Int>(j));
#if defined(GC_PROFILE_STANDARDIZER_SORT)
									++swap_count;
#endif
								}
								++i;
								--j;
							}
						}

						if (j - left < right - i) {
							if (left < j) {
								self(self, left, j);
							}
							left = i;
						} else {
							if (i < right) {
								self(self, i, right);
							}
							right = j;
						}
					}
					if (left < right) {
						insertion_sort_range(left, right);
					}
				};

				if constexpr (N_EDGES > 1) {
					quick_sort(quick_sort, 0, N_EDGES - 1);
				}

#if defined(GC_PROFILE_STANDARDIZER_SORT)
				const auto sort_stop = std::chrono::steady_clock::now();
				const auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(sort_stop - sort_start).count();
				gc_standardizer_sort_profile::insertion_sort_calls.fetch_add(1, std::memory_order_relaxed);
				gc_standardizer_sort_profile::insertion_sort_swaps.fetch_add(swap_count, std::memory_order_relaxed);
				gc_standardizer_sort_profile::insertion_sort_nanoseconds.fetch_add(static_cast<unsigned long long>(elapsed), std::memory_order_relaxed);
#endif
				return overallSign;
			}

			struct EdgeRadixEntry {
				Int u = 0;
				Int v = 0;
				Int original_index = 0;
			};

			template <typename KeyFn>
			static void stable_bucket_pass(
				const std::array<EdgeRadixEntry, N_EDGES>& input,
				std::array<EdgeRadixEntry, N_EDGES>& output,
				KeyFn key_fn
			) {
				std::array<Int, N_VERTICES> counts{};
				for (const auto& entry : input) {
					++counts[key_fn(entry)];
				}

				std::array<Int, N_VERTICES> offsets{};
				Int next_offset = 0;
				for (Int key = 0; key < N_VERTICES; ++key) {
					offsets[key] = next_offset;
					next_offset += counts[key];
				}

				std::array<Int, N_VERTICES> next_position = offsets;
				for (const auto& entry : input) {
					output[next_position[key_fn(entry)]++] = entry;
				}
			}

			static signedInt permutation_sign_from_destinations(
				const std::array<Int, N_EDGES>& destinations
			) {
				std::array<bool, N_EDGES> visited{};
				signedInt sign = 1;

				for (Int edge_index = 0; edge_index < N_EDGES; ++edge_index) {
					if (visited[edge_index]) {
						continue;
					}

					Int cycle_length = 0;
					Int current = edge_index;
					while (!visited[current]) {
						visited[current] = true;
						current = destinations[current];
						++cycle_length;
					}

					if ((cycle_length & 1) == 0) {
						sign *= -1;
					}
				}

				return sign;
			}

			signedInt sortEdgesRadix() {
				if constexpr (N_EDGES <= 1) {
					return 1;
				}

				std::array<EdgeRadixEntry, N_EDGES> original_entries{};
				for (Int edge_index = 0; edge_index < N_EDGES; ++edge_index) {
					const Int base = N_HAIR + 2 * edge_index;
					original_entries[edge_index] = EdgeRadixEntry{
						.u = half_edges[base],
						.v = half_edges[base + 1],
						.original_index = edge_index
					};
				}

				std::array<EdgeRadixEntry, N_EDGES> by_second{};
				std::array<EdgeRadixEntry, N_EDGES> sorted_entries{};
				stable_bucket_pass(original_entries, by_second, [](const EdgeRadixEntry& entry) {
					return entry.v;
				});
				stable_bucket_pass(by_second, sorted_entries, [](const EdgeRadixEntry& entry) {
					return entry.u;
				});

				if constexpr (SWAP_EDGE_SIGN == 1) {
					for (Int sorted_index = 0; sorted_index < N_EDGES; ++sorted_index) {
						const Int base = N_HAIR + 2 * sorted_index;
						half_edges[base] = sorted_entries[sorted_index].u;
						half_edges[base + 1] = sorted_entries[sorted_index].v;
					}
					return 1;
				}

				std::array<Int, N_EDGES> destinations{};
				for (Int sorted_index = 0; sorted_index < N_EDGES; ++sorted_index) {
					destinations[sorted_entries[sorted_index].original_index] = sorted_index;
					const Int base = N_HAIR + 2 * sorted_index;
					half_edges[base] = sorted_entries[sorted_index].u;
					half_edges[base + 1] = sorted_entries[sorted_index].v;
				}

				return permutation_sign_from_destinations(destinations);
			}

			signedInt directAndSortEdges() {
				signedInt sign = directEdges();
				return sign * sortEdgesRadix();
			}

			signedInt sortEdgesActive() {
				return sortEdgesRadix();
			}

			signedInt assignPermutedDirectedSortedEdges(
				const ThisGraph& source,
				const Permutation<N_VERTICES>& perm
			) {
				signedInt sign = 1;
				if constexpr (SWAP_VERTICES_SIGN == -1) {
					sign = perm.sign();
				}

				for (Int i = 0; i < N_HAIR; ++i) {
					half_edges[i] = perm[source.half_edges[i]];
				}

				std::array<EdgeRadixEntry, N_EDGES> original_entries{};
				for (Int edge_index = 0; edge_index < N_EDGES; ++edge_index) {
					const Int base = N_HAIR + 2 * edge_index;
					Int u = perm[source.half_edges[base]];
					Int v = perm[source.half_edges[base + 1]];
					if (u <= v) {
						original_entries[edge_index] = EdgeRadixEntry{u, v, edge_index};
					} else {
						original_entries[edge_index] = EdgeRadixEntry{v, u, edge_index};
						if constexpr (FLIP_EDGE_SIGN == -1) {
							sign *= -1;
						}
					}
				}

				if constexpr (N_EDGES <= 1) {
					if constexpr (N_EDGES == 1) {
						half_edges[N_HAIR] = original_entries[0].u;
						half_edges[N_HAIR + 1] = original_entries[0].v;
					}
					return sign;
				}

				std::array<EdgeRadixEntry, N_EDGES> by_second{};
				std::array<EdgeRadixEntry, N_EDGES> sorted_entries{};
				stable_bucket_pass(original_entries, by_second, [](const EdgeRadixEntry& entry) {
					return entry.v;
				});
				stable_bucket_pass(by_second, sorted_entries, [](const EdgeRadixEntry& entry) {
					return entry.u;
				});

				if constexpr (SWAP_EDGE_SIGN == 1) {
					for (Int sorted_index = 0; sorted_index < N_EDGES; ++sorted_index) {
						const Int base = N_HAIR + 2 * sorted_index;
						half_edges[base] = sorted_entries[sorted_index].u;
						half_edges[base + 1] = sorted_entries[sorted_index].v;
					}
					return sign;
				}

				std::array<Int, N_EDGES> destinations{};
				for (Int sorted_index = 0; sorted_index < N_EDGES; ++sorted_index) {
					destinations[sorted_entries[sorted_index].original_index] = sorted_index;
					const Int base = N_HAIR + 2 * sorted_index;
					half_edges[base] = sorted_entries[sorted_index].u;
					half_edges[base + 1] = sorted_entries[sorted_index].v;
				}
				sign *= permutation_sign_from_destinations(destinations);
				return sign;
			}

			signedInt permuteVertices(Permutation<N_VERTICES> perm) {
				for (auto& h : half_edges) {
					h = perm[h];
				}
				return (SWAP_VERTICES_SIGN == -1) ? perm.sign() : 1;
			}

			signedInt swapVertices(Int v, Int w) {
				if (v == w) return 1;
				for (Int& h : half_edges) {
					if (h == v) h = w;
					else if (h == w) h = v;
				}
				return SWAP_VERTICES_SIGN;
			}

			void print(std::ostream& out = std::cout) const {
				out << "printing graph:\n";

				if (N_OUT_HAIR > 0) {
					out << "out_hair: ";
					for (Int i = 0; i < N_OUT_HAIR; ++i) {
						out << +half_edges[i];
						if (i < N_OUT_HAIR - 1) out << ", ";
					}
					out << "\n";
				}

				if (N_IN_HAIR > 0) {
					out << "in_hair: ";
					for (Int i = 0; i < N_IN_HAIR; ++i) {
						out << +half_edges[N_OUT_HAIR + i];
						if (i < N_IN_HAIR - 1) out << ", ";
					}
					out << "\n";
				}

				out << "edges: ";
				for (Int e = 0; e < N_EDGES; ++e) {
					Int base = N_HAIR + 2 * e;
					out << "(" << +half_edges[base] << "," << +half_edges[base + 1] << ")";
					if (e < N_EDGES - 1) out << ", ";
				}
				out << "\n";

				out << "grade = " << custom_filter() << "\n";
			}


			bool operator==(Graph const& other) const {
				return half_edges == other.half_edges;
			}

			bool operator!=(Graph const& other) const {
				return !(*this == other);
			}

			std::array<Int, SIZE> half_edges{};


			Int half_edge(Int i) const {
				return half_edges[i];
			}

	};
