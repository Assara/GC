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
		// Both child vertices must reach min_valence, and at least one must
		// reach preserve_valence, including the new connecting edge.
		template <typename Collector>
		void split(
			Int vertex,
			Int preserve_valence,
			Int min_valence,
			Int max_loop_number,
			Collector& collector
		) const {
			// For connected graphs L = E-V+1; each shared edge adds one loop.
			constexpr int parent_loops = N_EDGES - N_VERTICES + 1;
			if (parent_loops > max_loop_number) return;
			const auto shared_limit = static_cast<std::size_t>(max_loop_number - parent_loops);
			const auto incidences = graph_.adjacent(vertex);
			const bool is_universal = incidences.size() == N_VERTICES - 1;
			std::array<Int, 2 * (N_EDGES + 1 + MAX_SHARED_EDGES)> edges{};
			std::copy(graph_.half_edges.begin(), graph_.half_edges.end(), edges.begin());

			// Keep the first exclusively assigned incidence on the old vertex
			// to identify assignments related by exchanging the two new vertices.
			// Shared edges before that incidence must still be enumerated.
			auto enumerate = [&](auto&& self, std::size_t index,
				std::size_t shared, std::size_t left, std::size_t right,
				bool distinguished) -> void {
				if (index == incidences.size()) {
					// The connecting edge adds one to each new vertex valence.
					if (left + 1 < min_valence || right + 1 < min_valence) return;
					const auto child_max_valence = std::max(left, right) + 1;
					if (child_max_valence < preserve_valence) return;
					if (!is_universal && child_max_valence > incidences.size()) return;
					const auto connection = graph_type::SIZE + 2 * shared;
					edges[connection] = vertex;
					edges[connection + 1] = N_VERTICES;
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

	private:
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
