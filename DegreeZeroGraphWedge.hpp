#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "coalgebra_utils.hpp"
#include "graph.hpp"
#include "VectorSpace/BasisElement.hpp"
#include "VectorSpace/LinComb.hpp"

template <Int TotalHalfEdges, Int NumberOfWedges>
class DegreeZeroGraphWedge {
	static_assert(NumberOfWedges > 0, "A wedge must contain at least one factor.");
	static_assert(TotalHalfEdges > 0, "A wedge must contain half-edges.");

	public:
		using ThisType = DegreeZeroGraphWedge<TotalHalfEdges, NumberOfWedges>;
		using Basis = BasisElement<ThisType, fieldType>;
		using LinComb = VectorSpace::LinComb<ThisType, fieldType>;

		static constexpr Int TOTAL_HALF_EDGES = TotalHalfEdges;
		static constexpr Int N_WEDGES = NumberOfWedges;
		static constexpr std::size_t TOTAL_HALF_EDGES_SIZE_T = static_cast<std::size_t>(TotalHalfEdges);
		static constexpr std::size_t N_WEDGES_SIZE_T = static_cast<std::size_t>(NumberOfWedges);

		DegreeZeroGraphWedge() {
			split_offsets.fill(0);
			split_offsets[N_WEDGES] = TOTAL_HALF_EDGES;
		}

		explicit DegreeZeroGraphWedge(
			const std::array<Int, N_WEDGES + 1>& split_offsets_,
			const std::array<Int, TOTAL_HALF_EDGES_SIZE_T>& half_edges_
		) : split_offsets(split_offsets_), half_edges(half_edges_) {}

		template <Int N_VERTICES, Int N_EDGES>
		static ThisType from_graph(const Graph<N_VERTICES, N_EDGES, 0, 0, 0, 1, fieldType>& graph) {
			static_assert(NumberOfWedges == 1, "from_graph is only available for one-factor wedges.");
			static_assert(static_cast<Int>(2 * N_EDGES) == TotalHalfEdges,
				"graph size does not match the wedge total half-edge count");
			static_assert(N_EDGES - 2 * N_VERTICES + 2 == 0,
				"from_graph expects a degree-zero GC_2 graph");

			ThisType wedge;
			wedge.split_offsets[0] = 0;
			wedge.split_offsets[1] = TotalHalfEdges;
			for (std::size_t i = 0; i < TOTAL_HALF_EDGES_SIZE_T; ++i) {
				wedge.half_edges[i] = graph.half_edges[i];
			}
			return wedge;
		}

		signedInt compare(const ThisType& other) const noexcept {
			for (std::size_t i = 0; i < N_WEDGES_SIZE_T + 1; ++i) {
				if (split_offsets[i] < other.split_offsets[i]) {
					return -1;
				}
				if (split_offsets[i] > other.split_offsets[i]) {
					return 1;
				}
			}

			for (std::size_t i = 0; i < TOTAL_HALF_EDGES_SIZE_T; ++i) {
				if (half_edges[i] < other.half_edges[i]) {
					return -1;
				}
				if (half_edges[i] > other.half_edges[i]) {
					return 1;
				}
			}

			return 0;
		}

		bool operator<(const ThisType& other) const noexcept {
			return compare(other) < 0;
		}

		bool operator==(const ThisType& other) const noexcept {
			return split_offsets == other.split_offsets
				&& half_edges == other.half_edges;
		}

		static void std(Basis& b) {
			b = canonized(b);
		}

		static Basis canonized(Basis& b) {
			Basis result = b;
			if (result.getCoefficient() == fieldType{0}) {
				return result;
			}
			auto& wedge = result.getValue();

			assert(wedge.has_valid_splits());

			for (std::size_t i = 0; i < N_WEDGES_SIZE_T; ++i) {
				result.multiplyCoefficient(wedge.standardize_factor(i));
			}

			if (wedge.sort_factors_and_count_sign() % 2 != 0) {
				result.multiplyCoefficient(fieldType{-1});
			}

			return result;
		}

		bool has_valid_splits() const noexcept {
			if (split_offsets[0] != 0 || split_offsets[N_WEDGES] != TOTAL_HALF_EDGES) {
				return false;
			}

			for (std::size_t i = 0; i < N_WEDGES_SIZE_T; ++i) {
				if (split_offsets[i] >= split_offsets[i + 1]) {
					return false;
				}
				if (!is_degree_zero_factor_half_edge_count(factor_half_edge_count(i))) {
					return false;
				}
			}

			return true;
		}

		Int factor_half_edge_count(std::size_t factor_index) const noexcept {
			return split_offsets[factor_index + 1] - split_offsets[factor_index];
		}

		std::size_t factor_offset(std::size_t factor_index) const noexcept {
			return split_offsets[factor_index];
		}

		std::size_t factor_edge_count(std::size_t factor_index) const noexcept {
			return static_cast<std::size_t>(factor_half_edge_count(factor_index)) / 2;
		}

		std::size_t factor_vertex_count(std::size_t factor_index) const noexcept {
			return (static_cast<std::size_t>(factor_half_edge_count(factor_index)) + 4) / 4;
		}

		const std::array<Int, N_WEDGES + 1>& splits() const noexcept {
			return split_offsets;
		}

		std::array<Int, N_WEDGES + 1>& splits() noexcept {
			return split_offsets;
		}

		const std::array<Int, TOTAL_HALF_EDGES_SIZE_T>& data() const noexcept {
			return half_edges;
		}

		std::array<Int, TOTAL_HALF_EDGES_SIZE_T>& data() noexcept {
			return half_edges;
		}

		using NextWedge = DegreeZeroGraphWedge<TotalHalfEdges, NumberOfWedges + 1>;
		using NextBasis = BasisElement<NextWedge, fieldType>;
		using NextLinComb = VectorSpace::LinComb<NextWedge, fieldType>;

		NextBasis contract_subgraph(std::size_t graph_index, const std::vector<Int>& selected_edges) const {
			if (graph_index >= N_WEDGES_SIZE_T) {
				return NextBasis(NextWedge{}, fieldType{0});
			}

			const Int source_half_edge_count = factor_half_edge_count(graph_index);
			const std::size_t source_edge_count = factor_edge_count(graph_index);
			const std::size_t source_vertex_count = factor_vertex_count(graph_index);
			const std::size_t source_offset = factor_offset(graph_index);

			if (selected_edges.empty()) {
				return NextBasis(NextWedge{}, fieldType{0});
			}

			std::vector<bool> is_selected_edge(source_edge_count, false);
			for (const Int edge : selected_edges) {
				const std::size_t edge_index = static_cast<std::size_t>(edge);
				if (edge_index >= source_edge_count || is_selected_edge[edge_index]) {
					return NextBasis(NextWedge{}, fieldType{0});
				}
				is_selected_edge[edge_index] = true;
			}

			const std::size_t subgraph_edge_count = selected_edges.size();
			const Int subgraph_half_edge_count = static_cast<Int>(2 * subgraph_edge_count);
			const Int quotient_half_edge_count = source_half_edge_count - subgraph_half_edge_count;

			if (!is_degree_zero_factor_half_edge_count(subgraph_half_edge_count)
				|| !is_degree_zero_factor_half_edge_count(quotient_half_edge_count)) {
				return NextBasis(NextWedge{}, fieldType{0});
			}

			std::vector<std::vector<Int>> selected_adjacency(source_vertex_count);
			std::vector<Int> selected_valence(source_vertex_count, 0);
			std::vector<bool> is_selected_vertex(source_vertex_count, false);

			for (std::size_t edge_index = 0; edge_index < source_edge_count; ++edge_index) {
				if (!is_selected_edge[edge_index]) {
					continue;
				}

				const Int u = half_edges[source_offset + 2 * edge_index];
				const Int v = half_edges[source_offset + 2 * edge_index + 1];

				is_selected_vertex[u] = true;
				is_selected_vertex[v] = true;
				++selected_valence[u];
				++selected_valence[v];
				selected_adjacency[u].push_back(v);
				selected_adjacency[v].push_back(u);
			}

			std::vector<Int> selected_vertices;
			selected_vertices.reserve(source_vertex_count);
			for (std::size_t v = 0; v < source_vertex_count; ++v) {
				if (is_selected_vertex[v]) {
					selected_vertices.push_back(static_cast<Int>(v));
				}
			}

			if (selected_vertices.empty()) {
				return NextBasis(NextWedge{}, fieldType{0});
			}

			for (const Int v : selected_vertices) {
				if (selected_valence[v] < 3) {
					return NextBasis(NextWedge{}, fieldType{0});
				}
			}

			if (static_cast<signedInt>(subgraph_edge_count)
				- 2 * static_cast<signedInt>(selected_vertices.size()) + 2 != 0) {
				return NextBasis(NextWedge{}, fieldType{0});
			}

			std::vector<bool> visited(source_vertex_count, false);
			std::vector<Int> stack{selected_vertices.front()};
			visited[selected_vertices.front()] = true;
			std::size_t visited_count = 0;

			while (!stack.empty()) {
				const Int v = stack.back();
				stack.pop_back();
				++visited_count;

				for (const Int w : selected_adjacency[v]) {
					if (!visited[w]) {
						visited[w] = true;
						stack.push_back(w);
					}
				}
			}

			if (visited_count != selected_vertices.size()) {
				return NextBasis(NextWedge{}, fieldType{0});
			}

			std::vector<Int> subgraph_vertex_map(source_vertex_count, static_cast<Int>(255));
			for (std::size_t i = 0; i < selected_vertices.size(); ++i) {
				subgraph_vertex_map[selected_vertices[i]] = static_cast<Int>(i);
			}

			std::array<Int, TOTAL_HALF_EDGES_SIZE_T> next_half_edges{};
			std::array<Int, NextWedge::N_WEDGES + 1> next_splits{};

			std::vector<Int> quotient_vertex_map(source_vertex_count, static_cast<Int>(255));
			Int next_outside_label = 0;
			for (std::size_t v = 0; v < source_vertex_count; ++v) {
				if (!is_selected_vertex[v]) {
					quotient_vertex_map[v] = next_outside_label++;
				}
			}
			const Int contracted_vertex_label = next_outside_label;

			fieldType coefficient{1};

			signedInt edge_inversions = 0;
			signedInt unselected_seen_to_right = 0;
			for (std::size_t idx = source_edge_count; idx-- > 0;) {
				if (is_selected_edge[idx]) {
					edge_inversions += unselected_seen_to_right;
				} else {
					++unselected_seen_to_right;
				}
			}
			if (edge_inversions % 2 != 0) {
				coefficient *= fieldType{-1};
			}

			const signedInt wedge_swaps = static_cast<signedInt>(N_WEDGES) - static_cast<signedInt>(graph_index);
			if (wedge_swaps % 2 != 0) {
				coefficient *= fieldType{-1};
			}

			std::size_t write_factor = 0;
			std::size_t write_offset = 0;
			next_splits[0] = 0;

			for (std::size_t original_factor = 0; original_factor < N_WEDGES_SIZE_T; ++original_factor) {
				if (original_factor != graph_index) {
					const Int factor_count = factor_half_edge_count(original_factor);
					const std::size_t read_offset = factor_offset(original_factor);
					for (std::size_t i = 0; i < static_cast<std::size_t>(factor_count); ++i) {
						next_half_edges[write_offset + i] = half_edges[read_offset + i];
					}
					write_offset += factor_count;
					next_splits[++write_factor] = static_cast<Int>(write_offset);
					continue;
				}

				for (std::size_t edge_index = 0; edge_index < source_edge_count; ++edge_index) {
					if (is_selected_edge[edge_index]) {
						continue;
					}

					const Int u = half_edges[source_offset + 2 * edge_index];
					const Int v = half_edges[source_offset + 2 * edge_index + 1];
					const Int new_u = is_selected_vertex[u] ? contracted_vertex_label : quotient_vertex_map[u];
					const Int new_v = is_selected_vertex[v] ? contracted_vertex_label : quotient_vertex_map[v];

					next_half_edges[write_offset++] = new_u;
					next_half_edges[write_offset++] = new_v;
				}
				next_splits[++write_factor] = static_cast<Int>(write_offset);
			}

			for (std::size_t edge_index = 0; edge_index < source_edge_count; ++edge_index) {
				if (!is_selected_edge[edge_index]) {
					continue;
				}

				const Int u = half_edges[source_offset + 2 * edge_index];
				const Int v = half_edges[source_offset + 2 * edge_index + 1];

				next_half_edges[write_offset++] = subgraph_vertex_map[u];
				next_half_edges[write_offset++] = subgraph_vertex_map[v];
			}
			next_splits[++write_factor] = static_cast<Int>(write_offset);

			assert(write_factor == NextWedge::N_WEDGES);
			assert(write_offset == TOTAL_HALF_EDGES_SIZE_T);

			return NextBasis(NextWedge(next_splits, next_half_edges), coefficient);
		}

		NextLinComb cobracket() const {
			NextLinComb result;

			for (std::size_t factor_index = 0; factor_index < N_WEDGES_SIZE_T; ++factor_index) {
				append_factor_cobracket(result, factor_index);
			}

			result.standardize_and_sort();
			return result;
		}

	private:
		static constexpr bool is_degree_zero_factor_half_edge_count(Int count) noexcept {
			return count > 0 && (count % 4) == 0;
		}

		template <Int FactorHalfEdges>
		static constexpr Int vertex_count_for_factor() noexcept {
			return static_cast<Int>((FactorHalfEdges + 4) / 4);
		}

		template <Int FactorHalfEdges>
		static constexpr Int edge_count_for_factor() noexcept {
			return static_cast<Int>(FactorHalfEdges / 2);
		}

		template <Int CurrentHalfEdges, typename Fn>
		static decltype(auto) dispatch_degree_zero_graph_type(Int factor_half_edge_count, Fn&& fn) {
			if (factor_half_edge_count == CurrentHalfEdges) {
				using GraphType = Graph<
					vertex_count_for_factor<CurrentHalfEdges>(),
					edge_count_for_factor<CurrentHalfEdges>(),
					0,
					0,
					0,
					1,
					fieldType
				>;
				return fn.template operator()<GraphType>();
			}

			if constexpr (CurrentHalfEdges + 4 <= TotalHalfEdges) {
				return dispatch_degree_zero_graph_type<CurrentHalfEdges + 4>(
					factor_half_edge_count,
					std::forward<Fn>(fn)
				);
			} else {
				assert(false && "Factor size is not a supported degree-zero GC_2 size.");
				using GraphType = Graph<
					vertex_count_for_factor<CurrentHalfEdges>(),
					edge_count_for_factor<CurrentHalfEdges>(),
					0,
					0,
					0,
					1,
					fieldType
				>;
				return fn.template operator()<GraphType>();
			}
		}

		fieldType standardize_factor(std::size_t factor_index) {
			const Int factor_count = factor_half_edge_count(factor_index);
			const std::size_t offset = factor_offset(factor_index);

			auto standardize = [&]<typename GraphType>() -> fieldType {
				GraphType graph;
				for (std::size_t i = 0; i < static_cast<std::size_t>(factor_count); ++i) {
					graph.half_edges[i] = half_edges[offset + i];
				}

				BasisElement<GraphType, fieldType> graph_basis(graph, fieldType{1});
				auto canon = GraphType::canonized(graph_basis);

				for (std::size_t i = 0; i < static_cast<std::size_t>(factor_count); ++i) {
					half_edges[offset + i] = canon.getValue().half_edges[i];
				}

				return canon.getCoefficient();
			};

			return dispatch_degree_zero_graph_type<4>(factor_count, standardize);
		}

		void append_factor_cobracket(NextLinComb& result, std::size_t factor_index) const {
			const Int factor_count = factor_half_edge_count(factor_index);
			dispatch_degree_zero_graph_type<4>(factor_count, [&]<typename GraphType>() {
				GraphType factor_graph;
				const std::size_t offset = factor_offset(factor_index);
				for (std::size_t i = 0; i < static_cast<std::size_t>(factor_count); ++i) {
					factor_graph.half_edges[i] = half_edges[offset + i];
				}

				const auto banks = coalgebra_utils::connected_triangle_banks(factor_graph);
				std::vector<Int> selected_edges;

				const std::function<void(std::size_t)> recurse = [&](std::size_t bank_index) {
					if (bank_index == banks.size()) {
						auto term = contract_subgraph(factor_index, selected_edges);
						if (term.getCoefficient() != fieldType{0}) {
							result.append_in_basis_order(term);
						}
						return;
					}

					recurse(bank_index + 1);

					const std::size_t old_size = selected_edges.size();
					selected_edges.insert(
						selected_edges.end(),
						banks[bank_index].begin(),
						banks[bank_index].end()
					);
					recurse(bank_index + 1);
					selected_edges.resize(old_size);
				};

				recurse(0);
				return 0;
			});
		}

		signedInt compare_factors(std::size_t lhs, std::size_t rhs) const noexcept {
			const Int lhs_count = factor_half_edge_count(lhs);
			const Int rhs_count = factor_half_edge_count(rhs);

			if (lhs_count < rhs_count) {
				return -1;
			}
			if (lhs_count > rhs_count) {
				return 1;
			}

			const std::size_t lhs_offset = factor_offset(lhs);
			const std::size_t rhs_offset = factor_offset(rhs);
			for (std::size_t i = 0; i < static_cast<std::size_t>(lhs_count); ++i) {
				const Int a = half_edges[lhs_offset + i];
				const Int b = half_edges[rhs_offset + i];
				if (a < b) {
					return -1;
				}
				if (a > b) {
					return 1;
				}
			}

			return 0;
		}

		signedInt sort_factors_and_count_sign() {
			std::array<std::size_t, N_WEDGES_SIZE_T> permutation{};
			for (std::size_t i = 0; i < N_WEDGES_SIZE_T; ++i) {
				permutation[i] = i;
			}

			signedInt inversion_count = 0;
			for (std::size_t i = 1; i < N_WEDGES_SIZE_T; ++i) {
				std::size_t j = i;
				while (j > 0 && compare_factors(permutation[j - 1], permutation[j]) > 0) {
					std::swap(permutation[j - 1], permutation[j]);
					++inversion_count;
					--j;
				}
			}

			const auto original_offsets = split_offsets;
			const auto original_half_edges = half_edges;

			std::size_t write_offset = 0;
			split_offsets[0] = 0;
			for (std::size_t i = 0; i < N_WEDGES_SIZE_T; ++i) {
				const std::size_t original_index = permutation[i];
				const Int factor_count =
					original_offsets[original_index + 1] - original_offsets[original_index];
				const std::size_t read_offset = original_offsets[original_index];

				for (std::size_t j = 0; j < static_cast<std::size_t>(factor_count); ++j) {
					half_edges[write_offset + j] = original_half_edges[read_offset + j];
				}

				write_offset += factor_count;
				split_offsets[i + 1] = static_cast<Int>(write_offset);
			}

			return inversion_count;
		}

		std::array<Int, N_WEDGES + 1> split_offsets{};
		std::array<Int, TOTAL_HALF_EDGES_SIZE_T> half_edges{};
};
