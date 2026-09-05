#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>

#include "graph.hpp"

namespace GraphGeneration {

// Legacy survival labels. The odd-vertex case includes edge reversals when
// the canonicalizer's OddEdgeFlips option is enabled (ordinary odd GC).
enum class final_graph_survival : std::uint8_t {
	none = 0,
	odd_edges = 1U << 0,
	even_edges_odd_vertices = 1U << 1,
	both = (1U << 0) | (1U << 1)
};

constexpr bool survives(
	final_graph_survival flags,
	final_graph_survival sign_case
) noexcept {
	return (static_cast<std::uint8_t>(flags)
		& static_cast<std::uint8_t>(sign_case)) != 0;
}

constexpr final_graph_survival without_survival(
	final_graph_survival flags,
	final_graph_survival sign_case
) noexcept {
	return static_cast<final_graph_survival>(
		static_cast<std::uint8_t>(flags)
		& ~static_cast<std::uint8_t>(sign_case)
	);
}

template <typename GraphType>
struct final_canonicalization_result {
	GraphType canonical_graph{};
	std::uint64_t automorphism_order = 0;
	final_graph_survival survival = final_graph_survival::none;

	constexpr bool survives_odd_edges() const noexcept {
		return survives(survival, final_graph_survival::odd_edges);
	}

	constexpr bool survives_even_edges_odd_vertices() const noexcept {
		return survives(
			survival,
			final_graph_survival::even_edges_odd_vertices
		);
	}
};

// Canonicalize a final, simple, hairless graph once and collect the metadata
// needed by both sign conventions. Automorphisms here mean vertex
// automorphisms.  This intentionally does not add formal tadpole flips or
// parallel-edge swaps; such graphs are rejected before the search.
// Preserve the legacy convention by default; the new GC pipeline opts into
// the half-edge sign required for odd GC.
template <typename GraphType, bool OddEdgeFlips = false>
class final_graph_canonicalizer {
	public:
		static_assert(GraphType::N_HAIR == 0,
			"final graph canonicalization requires a hairless graph");
		static_assert(GraphType::N_VERTICES_ > 0,
			"final graph canonicalization requires at least one vertex");

		using result_type = final_canonicalization_result<GraphType>;
		using standardizer_type = GraphStandardizer<
			GraphType::N_VERTICES_,
			GraphType::N_EDGES_,
			GraphType::N_OUT_HAIR_,
			GraphType::N_IN_HAIR_,
			GraphType::C_,
			GraphType::D_,
			typename GraphType::Field
		>;
		using isomorphism_type = typename standardizer_type::IsomorphismType;

		result_type operator()(const GraphType& input) const {
			if (!is_simple(input)) {
				throw std::invalid_argument(
					"final graph canonicalization requires a simple graph"
				);
			}

			standardizer_type standardizer;
			auto [attempts, valid_attempts]
				= standardizer.create_final_attempts4(input);
			if (valid_attempts.empty()) {
				throw std::logic_error("canonical search produced no labeling");
			}

			result_type result;
			bool have_best = false;
			signedInt reference_edge_parity = 1;
			signedInt reference_odd_parity = 1;

			for (const std::size_t attempt_index : valid_attempts) {
				GraphType candidate;
				isomorphism_type isomorphism;
				(void)candidate.assignPermutedDirectedSortedEdgesWithIsomorphism(
					input,
					attempts[attempt_index].create_vertex_permutation(),
					isomorphism
				);

				signedInt odd_parity = isomorphism.vertex_permutation_sign();
				if constexpr (OddEdgeFlips) {
					for (bool flipped : isomorphism.edge_flip_data())
						if (flipped) odd_parity = -odd_parity;
				}

				const bool new_best = !have_best
					|| result.canonical_graph.compare(candidate) < 0;
				if (new_best) {
					result.canonical_graph = std::move(candidate);
					result.automorphism_order = 1;
					result.survival = final_graph_survival::both;
					reference_edge_parity
						= isomorphism.edge_permutation_sign();
					reference_odd_parity = odd_parity;
					have_best = true;
					continue;
				}

				if (result.canonical_graph.compare(candidate) != 0) {
					continue;
				}

				++result.automorphism_order;
				if (isomorphism.edge_permutation_sign()
					!= reference_edge_parity) {
					result.survival = without_survival(
						result.survival,
						final_graph_survival::odd_edges
					);
				}
				if (odd_parity != reference_odd_parity) {
					result.survival = without_survival(
						result.survival,
						final_graph_survival::even_edges_odd_vertices
					);
				}
			}

			return result;
		}

	private:
		static bool is_simple(const GraphType& graph) noexcept {
			constexpr std::size_t vertex_count = GraphType::N_VERTICES_;
			constexpr std::size_t edge_slots
				= vertex_count * (vertex_count - 1) / 2;
			std::array<bool, edge_slots> seen{};

			for (Int edge = 0; edge < GraphType::N_EDGES_; ++edge) {
				auto [first, second] = graph.getEdge(edge);
				if (first == second) {
					return false;
				}
				if (second < first) {
					std::swap(first, second);
				}
				const std::size_t first_index = first;
				const std::size_t second_index = second;
				const std::size_t slot
					= first_index * (2 * vertex_count - first_index - 1) / 2
					+ second_index - first_index - 1;
				if (seen[slot]) {
					return false;
				}
				seen[slot] = true;
			}
			return true;
		}
};

template <typename GraphType>
final_canonicalization_result<GraphType> canonicalize_final_graph(
	const GraphType& graph
) {
	return final_graph_canonicalizer<GraphType>{}(graph);
}

} // namespace GraphGeneration
