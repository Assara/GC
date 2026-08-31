#pragma once

#include <algorithm>
#include <array>
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <optional>
#include <span>
#include <stdexcept>
#include <utility>
#include <vector>

#include "graph.hpp"

namespace GraphGeneration {

// Classification after applying all inexpensive generation filters.  Only the
// two active cases belong in the next-stage frontier.  Both admissible cases
// belong in the permanent (hairless) graph set.
enum class transient_classification : std::uint8_t {
	discard,
	active_transient,
	active_admissible,
	terminal_admissible
};

constexpr bool is_active(transient_classification classification) noexcept {
	return classification == transient_classification::active_transient
		|| classification == transient_classification::active_admissible;
}

constexpr bool is_admissible(transient_classification classification) noexcept {
	return classification == transient_classification::active_admissible
		|| classification == transient_classification::terminal_admissible;
}

// A one-hair graph used only during generation.  The hair is a distinguished
// root marker: it is deliberately excluded from all graph-theoretic valences,
// adjacency lists, stability checks, and split subsets.
template <
	Int N_VERTICES,
	Int N_EDGES,
	typename FieldType = fieldType
>
class transient_graph {
	public:
		static_assert(N_VERTICES > 0);

		using graph_type = Graph<N_VERTICES, N_EDGES, 1, 0, 0, 0, FieldType>;
		using hairless_graph_type = Graph<N_VERTICES, N_EDGES, 0, 0, 0, 0, FieldType>;
		using split_type = transient_graph<
			static_cast<Int>(N_VERTICES + 1),
			static_cast<Int>(N_EDGES + 1),
			FieldType
		>;

		static constexpr Int N_VERTICES_ = N_VERTICES;
		static constexpr Int N_EDGES_ = N_EDGES;

		struct analysis_type {
			std::array<Int, N_VERTICES> valences{};
			std::array<Int, N_VERTICES> reduced_valences{};
			Int maximum_valence = 0;
			Int maximum_reduced_valence = 0;
			Int maximum_valence_count = 0;
			bool connected = true;
			bool at_least_trivalent = true;
			bool simple = true;
			bool defects_at_root = true;
			bool root_is_maximum = false;
			bool root_is_unique_maximum = false;
			bool root_is_reduced_maximum = false;
		};

		transient_graph() {
			graph_.half_edges.fill(0);
			graph_.half_edges[0] = 0;
		}

		explicit transient_graph(graph_type graph)
			: graph_(std::move(graph)) {
			validate_root();
		}

		transient_graph(const hairless_graph_type& graph, Int root_vertex) {
			set_root(root_vertex);
			std::copy(
				graph.half_edges.begin(),
				graph.half_edges.end(),
				graph_.half_edges.begin() + 1
			);
		}

		const graph_type& graph() const noexcept {
			return graph_;
		}

		graph_type& graph() noexcept {
			return graph_;
		}

		Int root() const noexcept {
			return graph_.half_edges[0];
		}

		void set_root(Int root_vertex) {
			if (root_vertex >= N_VERTICES) {
				throw std::out_of_range("transient graph root");
			}
			graph_.half_edges[0] = root_vertex;
		}

		std::array<Int, N_VERTICES> edge_valence_array() const noexcept {
			std::array<Int, N_VERTICES> valences{};
			for (Int edge = 0; edge < N_EDGES; ++edge) {
				const auto [first, second] = graph_.getEdge(edge);
				++valences[static_cast<std::size_t>(first)];
				++valences[static_cast<std::size_t>(second)];
			}
			return valences;
		}

		std::array<Int, N_VERTICES> ordinary_valences() const noexcept {
			return edge_valence_array();
		}

		std::vector<Int> edge_incidence_positions(Int vertex) const {
			if (vertex >= N_VERTICES) {
				throw std::out_of_range("transient graph vertex");
			}
			std::vector<Int> result;
			result.reserve(static_cast<std::size_t>(2 * N_EDGES));
			// Position zero is the artificial hair and must never enter a split.
			for (Int position = 1; position < graph_type::SIZE; ++position) {
				if (graph_.half_edges[position] == vertex) {
					result.push_back(position);
				}
			}
			return result;
		}

		analysis_type analyze() const {
			validate_root();
			analysis_type result;
			std::array<std::bitset<N_VERTICES>, N_VERTICES> neighbors{};

			for (Int edge = 0; edge < N_EDGES; ++edge) {
				auto [first, second] = graph_.getEdge(edge);
				++result.valences[static_cast<std::size_t>(first)];
				++result.valences[static_cast<std::size_t>(second)];

				if (first == second) {
					result.simple = false;
					if (first != root()) {
						result.defects_at_root = false;
					}
					continue;
				}

				const bool duplicate = neighbors[static_cast<std::size_t>(first)].test(
					static_cast<std::size_t>(second)
				);
				if (duplicate) {
					result.simple = false;
					if (first != root() && second != root()) {
						result.defects_at_root = false;
					}
				}
				neighbors[static_cast<std::size_t>(first)].set(
					static_cast<std::size_t>(second)
				);
				neighbors[static_cast<std::size_t>(second)].set(
					static_cast<std::size_t>(first)
				);
			}

			for (Int vertex = 0; vertex < N_VERTICES; ++vertex) {
				const std::size_t index = static_cast<std::size_t>(vertex);
				result.reduced_valences[index] = static_cast<Int>(neighbors[index].count());
				result.maximum_valence = std::max(
					result.maximum_valence, result.valences[index]
				);
				result.maximum_reduced_valence = std::max(
					result.maximum_reduced_valence,
					result.reduced_valences[index]
				);
				result.at_least_trivalent &= result.valences[index] >= 3;
			}

			for (const Int valence : result.valences) {
				result.maximum_valence_count += valence == result.maximum_valence;
			}

			std::array<Int, N_VERTICES> stack{};
			std::array<bool, N_VERTICES> visited{};
			std::size_t stack_size = 1;
			std::size_t visited_count = 1;
			stack[0] = 0;
			visited[0] = true;
			while (stack_size != 0) {
				const Int vertex = stack[--stack_size];
				for (Int neighbor = 0; neighbor < N_VERTICES; ++neighbor) {
					if (neighbors[static_cast<std::size_t>(vertex)].test(
							static_cast<std::size_t>(neighbor))
						&& !visited[static_cast<std::size_t>(neighbor)]) {
						visited[static_cast<std::size_t>(neighbor)] = true;
						stack[stack_size++] = neighbor;
						++visited_count;
					}
				}
			}
			result.connected = visited_count == N_VERTICES;

			const std::size_t root_index = static_cast<std::size_t>(root());
			result.root_is_maximum
				= result.valences[root_index] == result.maximum_valence;
			result.root_is_unique_maximum
				= result.root_is_maximum && result.maximum_valence_count == 1;
			// Ties are intentionally allowed for the reduced valence.
			result.root_is_reduced_maximum
				= result.reduced_valences[root_index]
				== result.maximum_reduced_valence;
			return result;
		}

		transient_classification classification() const {
			const analysis_type state = analyze();
			if (!state.connected || !state.at_least_trivalent
				|| !state.root_is_maximum) {
				return transient_classification::discard;
			}

			if (!state.root_is_unique_maximum) {
				return state.simple
					? transient_classification::terminal_admissible
					: transient_classification::discard;
			}

			if (!state.root_is_reduced_maximum) {
				return transient_classification::discard;
			}

			if (state.simple) {
				return transient_classification::active_admissible;
			}
			if (!state.defects_at_root) {
				return transient_classification::discard;
			}
			return transient_classification::active_transient;
		}

		bool can_reduce_to_admissible() const {
			return classification() != transient_classification::discard;
		}

		bool is_admissible_after_removing_hair() const {
			const analysis_type state = analyze();
			return state.simple && state.connected && state.at_least_trivalent;
		}

		// This is a label-preserving conversion, not canonicalization.  A caller
		// inserting the result into the permanent set must run the hairless
		// standardizer; the rooted canonical labeling generally stops being
		// canonical after the marker is removed.
		hairless_graph_type without_hair() const noexcept {
			hairless_graph_type result;
			std::copy(
				graph_.half_edges.begin() + 1,
				graph_.half_edges.end(),
				result.half_edges.begin()
			);
			return result;
		}

		hairless_graph_type forget_root() const noexcept {
			return without_hair();
		}

		std::optional<hairless_graph_type> remove_hair_if_admissible() const {
			if (!is_admissible_after_removing_hair()) {
				return std::nullopt;
			}
			return without_hair();
		}

		// Split using indices into edge_incidence_positions(root()).  Index zero
		// is fixed on the old side, which chooses one orientation of an
		// unoriented split.  The marker follows the strictly more-valent child
		// and stays on the old child in a tie.
		split_type split_root(std::span<const Int> moved_incidence_indices) const {
			const std::vector<Int> adjacent = edge_incidence_positions(root());
			if (moved_incidence_indices.size() < 2
				|| moved_incidence_indices.size() + 2 > adjacent.size()) {
				throw std::invalid_argument("unstable transient root split");
			}

			std::vector<bool> selected(adjacent.size());
			for (const Int incidence_index : moved_incidence_indices) {
				if (incidence_index == 0
					|| static_cast<std::size_t>(incidence_index) >= adjacent.size()
					|| selected[static_cast<std::size_t>(incidence_index)]) {
					throw std::invalid_argument(
						"transient split indices must be distinct and exclude zero"
					);
				}
				selected[static_cast<std::size_t>(incidence_index)] = true;
			}

			const std::vector<Int> moved(
				moved_incidence_indices.begin(), moved_incidence_indices.end()
			);
			return split_root_unchecked(adjacent, moved);
		}

		// Enumerate every stable unoriented split, without applying generation
		// filters.  This is useful to differential code and to test the split
		// convention independently of the frontier policy.
		template <typename Emit>
		std::uint64_t for_each_root_split(Emit&& emit) const {
			const Int old_root = root();
			const std::vector<Int> adjacent = edge_incidence_positions(old_root);
			const std::size_t valence = adjacent.size();
			if (valence < 4) {
				return 0;
			}

			std::uint64_t emitted = 0;
			const Int maximum_incidence_index = static_cast<Int>(valence - 1);
			for (Int moved_count = 2;
				 moved_count < maximum_incidence_index;
				 ++moved_count) {
				std::vector<Int> moved = combutils::firstSubset(1, moved_count);
				do {
					std::invoke(emit, split_root_unchecked(adjacent, moved));
					++emitted;
				} while (combutils::nextSubset(moved, maximum_incidence_index));
			}
			return emitted;
		}

		// Enumerate only children that may contribute either to the permanent
		// admissible set or to a later generation stage.  Unlike the general
		// differential-style iterator above, this generator works with groups of
		// indistinguishable root incidences:
		//
		//   (n tadpoles), (k_1 edges to v_1), ..., (k_s edges to v_s).
		//
		// A relevant oriented split is determined by the number of tadpoles that
		// are separated and by one bit per non-tadpole group.  Moving a whole
		// tadpole to the frozen child would leave a tadpole away from the marked
		// root; moving two edges from one group would leave a parallel bundle
		// there.  Such children can be rejected without constructing them.  The
		// lowest-indexed members provide one deterministic representative of each
		// group allocation.
		template <typename Emit>
		std::uint64_t for_each_relevant_root_split(Emit&& emit) const {
			const Int old_root = root();
			const std::vector<Int> adjacent = edge_incidence_positions(old_root);
			const std::size_t valence = adjacent.size();
			if (valence < 4) {
				return 0;
			}

			std::array<Int, graph_type::SIZE> adjacent_index_by_position{};
			for (std::size_t index = 0; index < adjacent.size(); ++index) {
				adjacent_index_by_position[static_cast<std::size_t>(adjacent[index])]
					= static_cast<Int>(index);
			}

			struct tadpole_incidence_indices {
				Int first = 0;
				Int second = 0;
			};
			std::vector<tadpole_incidence_indices> tadpoles;
			std::array<std::vector<Int>, N_VERTICES> incidences_by_neighbor;
			for (Int edge = 0; edge < N_EDGES; ++edge) {
				const Int base = static_cast<Int>(graph_type::N_HAIR + 2 * edge);
				const Int first = graph_.half_edges[base];
				const Int second = graph_.half_edges[static_cast<Int>(base + 1)];
				if (first == old_root && second == old_root) {
					tadpoles.push_back({
						adjacent_index_by_position[base],
						adjacent_index_by_position[static_cast<Int>(base + 1)]
					});
				} else if (first == old_root) {
					incidences_by_neighbor[static_cast<std::size_t>(second)].push_back(
						adjacent_index_by_position[base]
					);
				} else if (second == old_root) {
					incidences_by_neighbor[static_cast<std::size_t>(first)].push_back(
						adjacent_index_by_position[static_cast<Int>(base + 1)]
					);
				}
			}

			std::vector<Int> bundle_vertices;
			std::vector<Int> bundle_sizes;
			for (Int vertex = 0; vertex < N_VERTICES; ++vertex) {
				if (!incidences_by_neighbor[static_cast<std::size_t>(vertex)].empty()) {
					bundle_vertices.push_back(vertex);
					bundle_sizes.push_back(static_cast<Int>(
						incidences_by_neighbor[static_cast<std::size_t>(vertex)].size()
					));
				}
			}
			const std::vector<Int> moved_count_bounds(bundle_sizes.size(), Int{1});

			const auto parent_valences = edge_valence_array();
			Int maximum_untouched_valence = 0;
			for (Int vertex = 0; vertex < N_VERTICES; ++vertex) {
				if (vertex != old_root) {
					maximum_untouched_valence = std::max(
						maximum_untouched_valence,
						parent_valences[static_cast<std::size_t>(vertex)]
					);
				}
			}

			std::uint64_t emitted = 0;
			for (std::size_t separated_tadpoles = 0;
				 separated_tadpoles <= tadpoles.size();
				 ++separated_tadpoles) {
				combutils::for_each_bounded_count_vector(
					moved_count_bounds,
					[&](std::span<const Int> moved_from_bundle) {
						std::size_t moved_count = separated_tadpoles;
						for (const Int count : moved_from_bundle) {
							moved_count += count;
						}

						// Both split vertices must be at least trivalent.  Orient an
						// unequal split so the larger child remains the old/root side.
						const std::size_t retained_count = valence - moved_count;
						if (moved_count < 2 || retained_count < 2
							|| moved_count > retained_count) {
							return;
						}

						const Int old_child_valence
							= static_cast<Int>(retained_count + 1);
						const Int new_child_valence
							= static_cast<Int>(moved_count + 1);
						if (old_child_valence < maximum_untouched_valence
							&& new_child_valence < maximum_untouched_valence) {
							return;
						}

						if (moved_count == retained_count) {
							// Equal split children can only survive as a simple terminal
							// graph.  Tadpoles make that impossible, and every parallel
							// group must have at most one edge on either side.
							if (!tadpoles.empty()) {
								return;
							}
							for (std::size_t i = 0; i < bundle_sizes.size(); ++i) {
								if (bundle_sizes[i] - moved_from_bundle[i] > 1) {
									return;
								}
							}
							if (!combutils::is_first_grouped_bipartition(
									moved_from_bundle, bundle_sizes
								)) {
								return;
							}
						}

						std::vector<Int> moved;
						moved.reserve(moved_count);
						for (std::size_t i = 0; i < separated_tadpoles; ++i) {
							// One fixed endpoint represents both directions of splitting
							// an undirected tadpole.
							moved.push_back(tadpoles[i].second);
						}
						for (std::size_t i = 0; i < bundle_vertices.size(); ++i) {
							if (moved_from_bundle[i] != 0) {
								moved.push_back(
									incidences_by_neighbor[static_cast<std::size_t>(
										bundle_vertices[i]
									)].front()
								);
							}
						}

						split_type child = split_root_unchecked(adjacent, moved);
						const transient_classification child_classification
							= child.classification();
						if (child_classification != transient_classification::discard) {
							std::invoke(
								emit,
								std::move(child),
								child_classification
							);
							++emitted;
						}
					}
				);
			}
			return emitted;
		}

	private:
		split_type split_root_unchecked(
			const std::vector<Int>& adjacent,
			const std::vector<Int>& moved
		) const {
			const Int old_root = root();
			typename graph_type::SplitGraph child_graph
				= graph_.splitGraph(old_root, adjacent, moved);
			split_type child(std::move(child_graph));

			const std::size_t moved_count = moved.size();
			const std::size_t old_child_valence
				= adjacent.size() - moved_count + 1;
			const std::size_t new_child_valence = moved_count + 1;
			child.set_root(
				new_child_valence > old_child_valence
					? N_VERTICES
					: old_root
			);
			return child;
		}

		void validate_root() const {
			if (root() >= N_VERTICES) {
				throw std::out_of_range("transient graph root");
			}
		}

		graph_type graph_{};
};

} // namespace GraphGeneration
