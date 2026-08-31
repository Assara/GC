#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <stdexcept>
#include <utility>

#include "GraphGeneration/SupportTransientGraph.hpp"

namespace GraphGeneration {

template <Int N_VERTICES, Int N_EDGES>
class unrooted_support_transient_graph {
	public:
		using support_type = support_transient_graph<
			N_VERTICES, N_EDGES, fieldType
		>;
		using split_type = unrooted_support_transient_graph<
			static_cast<Int>(N_VERTICES + 1),
			static_cast<Int>(N_EDGES + 1)
		>;
		using symbolic_split_type = support_symbolic_split<split_type>;
		using support_adjacency_array = typename support_type::support_adjacency_array;

		static constexpr Int N_VERTICES_ = N_VERTICES;
		static constexpr Int N_EDGES_ = N_EDGES;

		unrooted_support_transient_graph() = default;
		explicit unrooted_support_transient_graph(support_type support)
			: support_(std::move(support)) {}

		static unrooted_support_transient_graph rose() {
			static_assert(N_VERTICES == 1);
			return unrooted_support_transient_graph(support_type::rose());
		}

		const support_type& support() const noexcept { return support_; }
		support_type& support() noexcept { return support_; }

		std::array<Int, N_VERTICES> support_valences() const noexcept {
			std::array<Int, N_VERTICES> result{};
			const auto adjacency = support_.support_adjacency();
			for (Int vertex = 0; vertex < N_VERTICES; ++vertex) {
				result[static_cast<std::size_t>(vertex)] = static_cast<Int>(
					std::popcount(adjacency[static_cast<std::size_t>(vertex)])
				);
			}
			return result;
		}

		Int maximum_support_valence() const noexcept {
			const auto degree = support_valences();
			return *std::max_element(degree.begin(), degree.end());
		}

		Int maximum_support_valence_count() const noexcept {
			const auto degree = support_valences();
			const Int maximum = *std::max_element(degree.begin(), degree.end());
			return static_cast<Int>(std::count(degree.begin(), degree.end(), maximum));
		}

		std::pair<Int, Int> storage_partition() const noexcept {
			return {maximum_support_valence(), maximum_support_valence_count()};
		}

		unrooted_support_transient_graph permuted(
			const std::array<Int, N_VERTICES>& old_to_new
		) const {
			std::array<bool, N_VERTICES> seen{};
			for (const Int image : old_to_new) {
				if (image < 0 || image >= N_VERTICES
					|| seen[static_cast<std::size_t>(image)]) {
					throw std::invalid_argument("unrooted support permutation");
				}
				seen[static_cast<std::size_t>(image)] = true;
			}
			support_type result;
			support_.for_each_support_edge([&](Int first, Int second) {
				result.set_support_edge(
					old_to_new[static_cast<std::size_t>(first)],
					old_to_new[static_cast<std::size_t>(second)]
				);
			});
			result.validate();
			return unrooted_support_transient_graph(std::move(result));
		}

		bool has_active_root_choice() const {
			for (Int root = 0; root < N_VERTICES; ++root) {
				if (rooted_at(root).has_active_allocation()) {
					return true;
				}
			}
			return false;
		}

		support_type rooted_at(Int root) const {
			if (root < 0 || root >= N_VERTICES) {
				throw std::out_of_range("unrooted support root choice");
			}
			std::array<Int, N_VERTICES> old_to_new{};
			old_to_new[static_cast<std::size_t>(root)] = 0;
			Int next = 1;
			for (Int vertex = 0; vertex < N_VERTICES; ++vertex) {
				if (vertex != root) {
					old_to_new[static_cast<std::size_t>(vertex)] = next++;
				}
			}
			support_type result;
			support_.for_each_support_edge([&](Int first, Int second) {
				result.set_support_edge(
					old_to_new[static_cast<std::size_t>(first)],
					old_to_new[static_cast<std::size_t>(second)]
				);
			});
			result.validate();
			return result;
		}

		// Try every vertex which can be the maximum-valence root of an exact
		// surplus allocation.  This is deliberately broader than maximum support
		// degree: hidden copies can make a lower support-degree vertex maximal.
		template <typename Emit>
		std::uint64_t for_each_allowed_root_split(Emit&& emit) const {
			const Int parent_maximum = maximum_support_valence();
			const bool parent_has_universal_vertex
				= parent_maximum == static_cast<Int>(N_VERTICES - 1);
			std::uint64_t emitted = 0;
			for (Int root = 0; root < N_VERTICES; ++root) {
				const support_type rooted = rooted_at(root);
				if (!rooted.has_active_allocation()) {
					continue;
				}
				rooted.for_each_symbolic_relevant_root_split_impl(
					[&](typename support_type::symbolic_split_type split) {
						symbolic_split_type result{
							split_type(std::move(split.child)),
							split.parent_tadpoles,
							split.separated_tadpoles,
							split.selected_groups,
							split.exhausted_groups,
							split.exact_image_classifications,
							split.active_root_mask,
							split.has_simple_image,
							split.witness_count
						};
						// Order the augmentation before canonicalization.  Creating a
						// support leaf, or increasing the maximum support valence, is
						// needed only after a graph has reached a universal vertex.
						if (result.selected_groups == 0
							&& !parent_has_universal_vertex) {
							return;
						}
						if (result.child.maximum_support_valence() > parent_maximum
							&& !parent_has_universal_vertex) {
							return;
						}
						std::invoke(emit, std::move(result), root);
						++emitted;
					},
					false,
					parent_has_universal_vertex
				);
			}
			return emitted;
		}

		bool operator==(const unrooted_support_transient_graph&) const noexcept
			= default;
		bool operator<(const unrooted_support_transient_graph& other) const noexcept {
			return support_ < other.support_;
		}
		std::size_t hash_value() const noexcept { return support_.hash_value(); }
		std::size_t hash() const noexcept { return hash_value(); }
		bool empty() const noexcept { return support_.empty(); }

	private:
		support_type support_{};
};

} // namespace GraphGeneration

namespace std {
	template <Int V, Int E>
	struct hash<GraphGeneration::unrooted_support_transient_graph<V, E>> {
		std::size_t operator()(
			const GraphGeneration::unrooted_support_transient_graph<V, E>& graph
		) const noexcept {
			return graph.hash_value();
		}
	};
}
