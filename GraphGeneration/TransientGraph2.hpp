#pragma once

#include <algorithm>
#include <array>
#include <limits>
#include <utility>

#include "graph.hpp"

namespace GraphGeneration {

template <
	Int N_VERTICES,
	Int N_EDGES
>
class transient_graph2 {
	public:
		static_assert(N_VERTICES > 0);

		using graph_type = Graph<
			N_VERTICES,
			N_EDGES,
			0,
			0,
			0,
			0,
			fieldType
		>;
		// Inputs are simple graphs, so at most V-1 incident edges can be shared.
		static constexpr std::size_t MAX_SHARED_EDGES =
			std::min<std::size_t>(N_EDGES, N_VERTICES - 1);
		static_assert(N_VERTICES < std::numeric_limits<Int>::max());
		static_assert(2 * (N_EDGES + 1 + MAX_SHARED_EDGES)
			<= std::numeric_limits<Int>::max());

		template <std::size_t SharedEdges = 0>
		using split_graph_type = Graph<
			N_VERTICES + 1, N_EDGES + 1 + SharedEdges, 0, 0, 0, 0, fieldType
		>;

		static constexpr Int N_VERTICES_ = N_VERTICES;
		static constexpr Int N_EDGES_ = N_EDGES;

		transient_graph2() = default;

		explicit transient_graph2(graph_type graph)
			: graph_(std::move(graph)) {}

		const graph_type& graph() const noexcept {
			return graph_;
		}

		graph_type& graph() noexcept {
			return graph_;
		}

		// Split one vertex.  The caller supplies a valid vertex label; whether
		// that vertex is maximal is deliberately not part of this operation.
		// Collector provides add(child) for every split_graph_type<K>. Children
		// are passed as rvalues; the collector owns storage and allocation policy.
		// Both child vertices must reach min_valence, including the connecting edge.
		// At least one must reach preserve_valence; neither may exceed max_valence.
		// Child loop numbers lie in the inclusive [min_loop_number, max_loop_number] range.
		template <typename Collector>
		void split(
			Int vertex,
			Int preserve_valence,
			Int min_valence,
			Int max_valence,
			Int min_loop_number,
			Int max_loop_number,
			Collector& collector
		) const {
			split_impl<false>(vertex, preserve_valence, min_valence, max_valence,
				min_loop_number, max_loop_number, collector);
		}

		// The caller selects this entry point when vertex has bivalent neighbours.
		// Only emit assignments retaining a common non-bivalent neighbour for all bivalents.
		// Create a bivalent split vertex, or leave no bivalent vertices in the child.
		template <typename Collector>
		void split_with_bivalent_neighbours(
			Int vertex,
			Int preserve_valence,
			Int min_valence,
			Int max_valence,
			Int min_loop_number,
			Int max_loop_number,
			Collector& collector
		) const {
			split_impl<true>(vertex, preserve_valence, min_valence, max_valence,
				min_loop_number, max_loop_number, collector);
		}

		// For each incident edge vertex--v, retain it and add vertex--b--v.
		// Emit one graph per edge; each adds one vertex, two edges, and one loop.
		// Only the chosen vertex is checked against the inclusive valence bounds.
		// The caller owns standardization and deduplication, as for split().
		template <typename Collector>
		void create_bivalent_vertices(
			Int vertex,
			Int preserve_valence,
			Int max_valence,
			Int min_loop_number,
			Int max_loop_number,
			Collector& collector
		) const {
			constexpr int child_loops = N_EDGES - N_VERTICES + 2;
			if (child_loops < min_loop_number || child_loops > max_loop_number) return;
			const auto incidences = graph_.adjacent(vertex);
			const auto resulting_valence = incidences.size() + 1;
			if (resulting_valence < preserve_valence || resulting_valence > max_valence) return;
			for (const auto position : incidences) {
				const auto other = position % 2 == 0 ? position + 1 : position - 1;
				split_graph_type<1> child;
				std::copy(graph_.half_edges.begin(), graph_.half_edges.end(), child.half_edges.begin());
				child.setEdge(N_EDGES, vertex, N_VERTICES);
				child.setEdge(N_EDGES + 1, N_VERTICES, graph_.half_edges[other]);
				collector.add(std::move(child));
			}
		}

	private:
		template <bool RequireRoot, typename Collector>
		void split_impl(Int vertex, Int preserve_valence, Int min_valence, Int max_valence,
			Int min_loop_number, Int max_loop_number, Collector& collector) const {
			if (min_valence > max_valence || preserve_valence > max_valence) return;
			// For connected graphs L = E-V+1; each shared edge adds one loop.
			constexpr int parent_loops = N_EDGES - N_VERTICES + 1;
			if (min_loop_number > max_loop_number || parent_loops > max_loop_number) return;
			const auto shared_minimum = static_cast<std::size_t>(
				std::max(0, static_cast<int>(min_loop_number) - parent_loops));
			const auto shared_limit = static_cast<std::size_t>(max_loop_number - parent_loops);
			const auto incidences = graph_.adjacent(vertex);
			std::array<Int, 2 * (N_EDGES + 1 + MAX_SHARED_EDGES)> edges{};
			std::copy(graph_.half_edges.begin(), graph_.half_edges.end(), edges.begin());

			// Keep the first exclusively assigned incidence on the old vertex
			// to identify assignments related by exchanging the two new vertices.
			// Shared edges before that incidence must still be enumerated.
			auto enumerate = [&](auto&& self, std::size_t index,
				std::size_t shared, std::size_t left, std::size_t right,
				bool distinguished) -> void {
				// Assigned incidences cannot be removed later in this branch.
				if (left + 1 > max_valence || right + 1 > max_valence) return;
				// Stop as soon as even sharing every remaining incidence cannot reach the minimum.
				if (shared + incidences.size() - index < shared_minimum) return;
				if (index == incidences.size()) {
					// The connecting edge adds one to each new vertex valence.
					if (left + 1 < min_valence || right + 1 < min_valence) return;
					if (std::max(left, right) + 1 < preserve_valence) return;
					const auto connection = graph_type::SIZE + 2 * shared;
					edges[connection] = vertex;
					edges[connection + 1] = N_VERTICES;
					if constexpr (RequireRoot) {
						if (!bivalent_assignment_allowed(edges, connection + 2,
							left == 1 || right == 1)) return;
					}
					append_split(shared, edges, collector,
						std::make_index_sequence<MAX_SHARED_EDGES + 1>{});
					return;
				}

				const auto position = incidences[index];
				edges[position] = vertex;
				self(self, index + 1, shared, left + 1, right, true);
				if (distinguished) {
					edges[position] = N_VERTICES;
					self(self, index + 1, shared, left, right + 1, true);
				}

				if (shared == shared_limit) return;
				edges[position] = vertex;
				const auto other = position % 2 == 0 ? position + 1 : position - 1;
				const auto extra = graph_type::SIZE + 2 * shared;
				edges[extra] = N_VERTICES;
				edges[extra + 1] = graph_.half_edges[other];
				self(self, index + 1, shared + 1, left + 1, right + 1,
					distinguished);
			};
			enumerate(enumerate, 0, 0, 0, 0, false);
		}


		// Inspect the completed incidence assignment before constructing any child graph.
		template <std::size_t Size>
		static bool bivalent_assignment_allowed(const std::array<Int, Size>& edges,
			std::size_t used, bool creates_bivalent) {
			std::array<Int, N_VERTICES + 1> valences{};
			for (std::size_t i = 0; i < used; ++i) ++valences[edges[i]];
			const auto bivalents = std::count(valences.begin(), valences.end(), Int{2});
			if (bivalents == 0) return true;
			if (!creates_bivalent) return false;
			std::array<Int, N_VERTICES + 1> bivalent_neighbours{};
			for (std::size_t i = 0; i < used; i += 2) {
				const auto a = edges[i], b = edges[i + 1];
				bivalent_neighbours[a] += valences[b] == 2;
				bivalent_neighbours[b] += valences[a] == 2;
			}
			for (std::size_t v = 0; v < valences.size(); ++v)
				if (valences[v] != 2 && bivalent_neighbours[v] == bivalents) return true;
			return false;
		}

		template <typename Collector, std::size_t Size, std::size_t... SharedEdges>
		static void append_split(std::size_t shared, const std::array<Int, Size>& edges,
			Collector& collector, std::index_sequence<SharedEdges...>) {
			auto append = [&]<std::size_t K>() {
				if (shared != K) return;
				split_graph_type<K> child;
				std::copy_n(edges.begin(), child.SIZE, child.half_edges.begin());
				collector.add(std::move(child));
			};
			(append.template operator()<SharedEdges>(), ...);
		}

		graph_type graph_{};
};

} // namespace GraphGeneration
