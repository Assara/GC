#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "GraphGeneration/RootedTransientGraph.hpp"

namespace GraphGeneration {

template <typename ChildGraph>
struct support_symbolic_split {
	ChildGraph child;
	std::size_t parent_tadpoles = 0;
	std::size_t separated_tadpoles = 0;
	std::uint64_t selected_groups = 0;
	std::uint64_t exhausted_groups = 0;
	std::uint8_t exact_image_classifications = 0;
	// Bits use the raw labels of child.  A set bit v means that a concrete
	// one-step image exists whose unique maximum is v and whose remaining
	// defects can all be represented with v as the distinguished root.
	std::uint64_t active_root_mask = 0;
	// A simple proper image is admissible independently of which raw child
	// label happened to be the parent root.
	bool has_simple_image = false;
	std::uint64_t witness_count = 1;
};

// Experimental quotient of rooted_transient_graph.
//
// Only the rooted simple support is retained.  If the support contains S edges,
// TOTAL_SURPLUS = N_EDGES-S represents every allocation
//
//     t + sum_{0--v in support} x_v = TOTAL_SURPLUS,
//
// where t is the number of root tadpoles and the exact multiplicity of 0--v is
// 1+x_v.  Thus one value denotes a family of exact transient graphs.  In
// particular, transitions below have existential/quotient semantics: after a
// child support is formed, constraints carried by the witnessing allocation
// are intentionally forgotten.
template <Int N_VERTICES, Int N_EDGES, typename FieldType = fieldType>
class support_transient_graph {
	public:
		static_assert(N_VERTICES > 0);
		static_assert(N_VERTICES <= 64,
			"support transient masks support at most 64 vertices");

		using field_type = FieldType;
		using exact_type
			= rooted_transient_graph<N_VERTICES, N_EDGES, FieldType>;
		using split_type = support_transient_graph<
			static_cast<Int>(N_VERTICES + 1),
			static_cast<Int>(N_EDGES + 1),
			FieldType
		>;
		using vertex_mask = std::uint64_t;
		using classification_mask_type = std::uint8_t;
		using symbolic_split_type = support_symbolic_split<split_type>;

		static constexpr Int N_VERTICES_ = N_VERTICES;
		static constexpr Int N_EDGES_ = N_EDGES;
		static constexpr std::size_t N_NONROOT_VERTICES
			= static_cast<std::size_t>(N_VERTICES - 1);
		static constexpr std::size_t SUPPORT_BIT_COUNT
			= static_cast<std::size_t>(N_VERTICES)
				* static_cast<std::size_t>(N_VERTICES - 1) / 2;
		static constexpr std::size_t SUPPORT_WORD_COUNT
			= (SUPPORT_BIT_COUNT + 63) / 64;

		using support_word_array
			= std::array<std::uint64_t, SUPPORT_WORD_COUNT>;
		using support_adjacency_array
			= std::array<std::uint64_t, static_cast<std::size_t>(N_VERTICES)>;
		// old_to_new_local[i] is the new local label of old full vertex i+1.
		using nonroot_permutation
			= std::array<Int, N_NONROOT_VERTICES>;

		constexpr support_transient_graph() noexcept = default;

		explicit support_transient_graph(support_word_array support_words)
			: support_words_(std::move(support_words)) {
			validate();
		}

		explicit support_transient_graph(const exact_type& source)
			: support_transient_graph(from_exact(source)) {}

		static support_transient_graph rose() {
			static_assert(N_VERTICES == 1);
			return {};
		}

		// Seed for the vertex-irreducible generator.  At loop order L this is
		// instantiated with E=L+2, so the rooted triangle carries L-1 units of
		// allocation-free surplus.
		static support_transient_graph triangle() {
			static_assert(N_VERTICES == 3);
			static_assert(N_EDGES >= 3);
			support_transient_graph result;
			result.set_support_edge_unchecked(0, 1, true);
			result.set_support_edge_unchecked(0, 2, true);
			result.set_support_edge_unchecked(1, 2, true);
			result.validate();
			return result;
		}

		// Convenience constructor for the four-vertex wheel W3=K4.  It is not
		// the production seed of the immutable-root generator.
		static support_transient_graph wheel3() {
			static_assert(N_VERTICES == 4);
			static_assert(N_EDGES >= 6);
			support_transient_graph result;
			for (Int first = 0; first < N_VERTICES; ++first) {
				for (Int second = static_cast<Int>(first + 1);
					 second < N_VERTICES;
					 ++second) {
					result.set_support_edge_unchecked(first, second, true);
				}
			}
			result.validate();
			return result;
		}

		static support_transient_graph from_exact(const exact_type& source) {
			source.validate();
			support_transient_graph result;
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				if (source.root_multiplicity(vertex) != 0) {
					result.set_support_edge_unchecked(0, vertex, true);
				}
			}
			source.for_each_core_edge([&](Int first, Int second) {
				result.set_support_edge_unchecked(first, second, true);
			});
			result.validate();
			return result;
		}

		static constexpr Int root() noexcept { return 0; }

		constexpr const support_word_array& support_words() const noexcept {
			return support_words_;
		}

		bool has_support_edge(Int first, Int second) const {
			const std::size_t bit = support_bit_index(first, second);
			return (support_words_[bit / 64]
				& (std::uint64_t{1} << (bit % 64))) != 0;
		}

		bool has_root_support(Int vertex) const {
			if (vertex == 0 || vertex >= N_VERTICES) {
				return false;
			}
			return has_support_edge(0, vertex);
		}

		bool has_core_support_edge(Int first, Int second) const {
			if (first == 0 || second == 0) {
				return false;
			}
			return has_support_edge(first, second);
		}

		void set_support_edge(Int first, Int second, bool present = true) {
			set_support_edge_unchecked(first, second, present);
		}

		void set_root_support(Int vertex, bool present = true) {
			if (vertex == 0 || vertex >= N_VERTICES) {
				throw std::out_of_range("support transient root neighbor");
			}
			set_support_edge_unchecked(0, vertex, present);
		}

		void set_core_support_edge(
			Int first, Int second, bool present = true
		) {
			if (first == 0 || second == 0) {
				throw std::invalid_argument("support core edge touches root");
			}
			set_support_edge_unchecked(first, second, present);
		}

		std::size_t support_edge_count() const noexcept {
			std::size_t result = 0;
			for (const std::uint64_t word : support_words_) {
				result += std::popcount(word);
			}
			return result;
		}

		std::size_t total_surplus() const {
			const std::size_t support_count = support_edge_count();
			if (support_count > static_cast<std::size_t>(N_EDGES)) {
				throw std::logic_error("support has more edges than expanded graph");
			}
			return static_cast<std::size_t>(N_EDGES) - support_count;
		}

		std::size_t surplus_edge_count() const { return total_surplus(); }

		support_adjacency_array support_adjacency() const noexcept {
			support_adjacency_array result{};
			for_each_support_edge([&](Int first, Int second) {
				result[static_cast<std::size_t>(first)]
					|= std::uint64_t{1} << second;
				result[static_cast<std::size_t>(second)]
					|= std::uint64_t{1} << first;
			});
			return result;
		}

		vertex_mask root_support_mask() const noexcept {
			return support_adjacency()[0];
		}

		Int support_degree(Int vertex) const {
			validate_vertex(vertex);
			return static_cast<Int>(
				std::popcount(support_adjacency()[static_cast<std::size_t>(vertex)])
			);
		}

		bool has_minimum_support_degree_two() const noexcept {
			const support_adjacency_array adjacency = support_adjacency();
			for (const vertex_mask neighbors : adjacency) {
				if (std::popcount(neighbors) < 2) {
					return false;
				}
			}
			return true;
		}

		bool has_minimum_support_degree_three() const noexcept {
			const support_adjacency_array adjacency = support_adjacency();
			for (const vertex_mask neighbors : adjacency) {
				if (std::popcount(neighbors) < 3) {
					return false;
				}
			}
			return true;
		}

		// A low-valence support vertex is repairable only when it is adjacent to
		// the temporary root.  Hidden copies of that bundle supply its deficit
		// 3-degree.  In particular a support leaf needs two surplus copies.
		bool has_repairable_support_valences() const noexcept {
			const support_adjacency_array adjacency = support_adjacency();
			if constexpr (N_VERTICES > 1) {
				if (std::popcount(adjacency[0]) < 1) {
					return false;
				}
			}
			std::size_t required_surplus = 0;
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				const std::size_t degree = std::popcount(
					adjacency[static_cast<std::size_t>(vertex)]
				);
				if (degree == 0) {
					return false;
				}
				if (degree < 3) {
					if ((adjacency[0]
							& (std::uint64_t{1} << vertex)) == 0) {
						return false;
					}
					required_surplus += 3 - degree;
				}
			}
			return required_surplus <= total_surplus();
		}

		bool has_vertex_irreducible_support() const noexcept {
			return support_is_vertex_irreducible(support_adjacency());
		}

		// Test the canonical-parent prefilter on raw labels.  Contracting 0--v
		// changes the root support degree to the degree of the corresponding
		// contracted parent.  Only contractions which produce an active,
		// support-repairable generator parent participate in the comparison.
		// Tied maximizers are all retained; canonization/map deduplication still
		// resolves those cases.
		bool root_contraction_is_maximal(Int distinguished_vertex) const {
			validate_vertex(distinguished_vertex);
			if (distinguished_vertex == 0
				|| !has_root_support(distinguished_vertex)) {
				return false;
			}
			if constexpr (N_VERTICES <= 1 || N_EDGES <= 0) {
				return false;
			} else {
				const auto distinguished_parent
					= contract_root_neighbor(distinguished_vertex);
				if (!contraction_parent_is_eligible(
						distinguished_parent, distinguished_vertex
					)) {
					return false;
				}
				const std::size_t distinguished_degree
					= distinguished_parent.support_degree(0);
				vertex_mask neighbors = root_support_mask();
				while (neighbors != 0) {
					const Int vertex = static_cast<Int>(
						std::countr_zero(neighbors)
					);
					neighbors &= neighbors - 1;
					if (vertex == distinguished_vertex) {
						continue;
					}
					const auto parent = contract_root_neighbor(vertex);
					if (contraction_parent_is_eligible(parent, vertex)
						&& parent.support_degree(0) > distinguished_degree) {
						return false;
					}
				}
				return true;
			}
		}

		template <typename Emit>
		void for_each_support_edge(Emit&& emit) const {
			for (Int first = 0; first < N_VERTICES; ++first) {
				for (Int second = static_cast<Int>(first + 1);
					 second < N_VERTICES;
					 ++second) {
					const std::size_t bit = support_bit_index_unchecked(first, second);
					if ((support_words_[bit / 64]
							& (std::uint64_t{1} << (bit % 64))) != 0) {
						std::invoke(emit, first, second);
					}
				}
			}
		}

		template <typename Emit>
		void for_each_root_neighbor(Emit&& emit) const {
			vertex_mask neighbors = root_support_mask();
			while (neighbors != 0) {
				const Int vertex = static_cast<Int>(std::countr_zero(neighbors));
				std::invoke(emit, vertex);
				neighbors &= neighbors - 1;
			}
		}

		void validate() const {
			if constexpr (SUPPORT_WORD_COUNT != 0 && SUPPORT_BIT_COUNT % 64 != 0) {
				constexpr std::uint64_t used_mask
					= (std::uint64_t{1} << (SUPPORT_BIT_COUNT % 64)) - 1;
				if ((support_words_.back() & ~used_mask) != 0) {
					throw std::invalid_argument("bits outside support triangle");
				}
			}
			if (support_edge_count() > static_cast<std::size_t>(N_EDGES)) {
				throw std::invalid_argument(
					"support has more edges than expanded graph"
				);
			}
		}

		bool structurally_valid() const noexcept {
			if constexpr (SUPPORT_WORD_COUNT != 0 && SUPPORT_BIT_COUNT % 64 != 0) {
				constexpr std::uint64_t used_mask
					= (std::uint64_t{1} << (SUPPORT_BIT_COUNT % 64)) - 1;
				if ((support_words_.back() & ~used_mask) != 0) {
					return false;
				}
			}
			return support_edge_count() <= static_cast<std::size_t>(N_EDGES);
		}

		bool has_simple_expansion() const noexcept {
			return support_edge_count() == static_cast<std::size_t>(N_EDGES);
		}

		exact_type expand_simple() const {
			validate();
			if (!has_simple_expansion()) {
				throw std::logic_error(
					"support family is not a singleton simple graph"
				);
			}
			exact_type result;
			for_each_support_edge([&](Int first, Int second) {
				if (first == 0) {
					result.set_root_multiplicity(second, 1);
				} else {
					result.set_core_edge(first, second);
				}
			});
			result.validate();
			return result;
		}

		exact_type to_simple_rooted_transient() const { return expand_simple(); }

		support_transient_graph permuted_nonroot(
			const nonroot_permutation& old_to_new_local
		) const {
			validate_nonroot_permutation(old_to_new_local);
			return permuted_nonroot_unchecked(old_to_new_local);
		}

		signedInt compare(const support_transient_graph& other) const noexcept {
			for (std::size_t bit = 0; bit < SUPPORT_BIT_COUNT; ++bit) {
				const bool ours = (support_words_[bit / 64]
					& (std::uint64_t{1} << (bit % 64))) != 0;
				const bool theirs = (other.support_words_[bit / 64]
					& (std::uint64_t{1} << (bit % 64))) != 0;
				if (ours != theirs) {
					return ours ? 1 : -1;
				}
			}
			return 0;
		}

		bool operator==(const support_transient_graph&) const noexcept = default;
		bool operator<(const support_transient_graph& other) const noexcept {
			return compare(other) < 0;
		}

		std::size_t hash_value() const noexcept {
			std::uint64_t hash = 0x9e3779b97f4a7c15ULL;
			for (std::uint64_t value : support_words_) {
				value += 0x9e3779b97f4a7c15ULL;
				value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
				value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
				value ^= value >> 31;
				hash ^= value + 0x9e3779b97f4a7c15ULL
					+ (hash << 6) + (hash >> 2);
			}
			return static_cast<std::size_t>(hash);
		}

		std::size_t hash() const noexcept { return hash_value(); }

		bool empty() const noexcept {
			if constexpr (N_VERTICES == 1) {
				// The zero-bit value is the one-vertex rose, not an empty slot.
				return false;
			} else {
				return std::ranges::all_of(
					support_words_,
					[](std::uint64_t word) { return word == 0; }
				);
			}
		}

		static constexpr classification_mask_type classification_bit(
			transient_classification classification
		) noexcept {
			return static_cast<classification_mask_type>(
				classification_mask_type{1}
				<< static_cast<unsigned>(classification)
			);
		}

		// Classifications are existential over the represented allocations.  More
		// than one bit can be set: the same support family may contain both active
		// and discarded exact graphs.
		classification_mask_type classification_mask() const {
			validate();
			classification_mask_type result = 0;
			const std::size_t surplus = total_surplus();
			for (std::size_t tadpoles = 0; tadpoles <= surplus; ++tadpoles) {
				result = static_cast<classification_mask_type>(
					result | classification_mask_for_domain(
						tadpoles, false, 0, 0
					)
				);
			}
			return result;
		}

		classification_mask_type classification_mask_for_tadpoles(
			std::size_t tadpoles
		) const {
			validate();
			return classification_mask_for_domain(tadpoles, false, 0, 0);
		}

		// Restrict one root bundle to exact excess `excess` while leaving the
		// other bundle excesses existential.  Symbolic splitting uses this for
		// the new 0--new_vertex bundle, whose excess is known to be c.
		classification_mask_type classification_mask_with_fixed_root_excess(
			std::size_t tadpoles,
			Int vertex,
			std::size_t excess
		) const {
			validate();
			if (!has_root_support(vertex)) {
				return 0;
			}
			return classification_mask_for_domain(
				tadpoles, true, vertex, excess
			);
		}

		bool has_classification(
			transient_classification classification
		) const {
			return (classification_mask() & classification_bit(classification)) != 0;
		}

		bool has_active_allocation() const {
			const classification_mask_type mask = classification_mask();
			return (mask & classification_bit(
				transient_classification::active_transient
			)) != 0 || (mask & classification_bit(
				transient_classification::active_admissible
			)) != 0;
		}

		bool has_admissible_allocation() const {
			const classification_mask_type mask = classification_mask();
			return (mask & classification_bit(
				transient_classification::active_admissible
			)) != 0 || (mask & classification_bit(
				transient_classification::terminal_admissible
			)) != 0;
		}

		bool can_reduce_to_admissible() const {
			return (classification_mask()
				& ~classification_bit(transient_classification::discard)) != 0;
		}

		// Enumerate quotient children.  A is the set of root bundles from which
		// one edge is moved; D_exhaust is the subset whose exact multiplicity was
		// one, so its old root-support edge disappears.  For fixed (A,D_exhaust)
		// the child support is independent of the witnessing tadpole counts t,c.
		// We therefore scan t,c internally and emit that labeled child support at
		// most once.  The first witness is retained for diagnostics and the exact
		// classification bits are ORed over every feasible witness.
		//
		// The emitted child is deliberately the whole support family.  In
		// particular, the exact constraint "new bundle excess == c" is not part of
		// its key.  This path-switching over-approximation is the experiment being
		// tested; exact_image_classifications describes only the actual one-step
		// images, not every allocation represented by child.
		template <typename Emit>
		std::uint64_t for_each_symbolic_relevant_root_split_impl(
			Emit&& emit,
			bool apply_max_contraction_prefilter,
			bool allow_special_split_for_nonuniversal_root = false
		) const {
			validate();
			// Bivalent root-neighbours are permitted only when the hidden surplus
			// can supply the missing exact incidence in each corresponding bundle.
			if (!has_repairable_support_valences()) {
				return 0;
			}
			const std::size_t surplus = total_surplus();
			const support_adjacency_array parent_adjacency
				= support_adjacency();
			const vertex_mask root_neighbors = parent_adjacency[0];
			const std::size_t bundle_count = std::popcount(root_neighbors);
			const bool root_is_universal
				= bundle_count == static_cast<std::size_t>(N_VERTICES - 1);
			const Int new_vertex = N_VERTICES;
			std::uint64_t emitted = 0;

			vertex_mask selected = root_neighbors;
			while (true) {
				vertex_mask exhausted = selected;
				while (true) {
					const std::size_t selected_count = std::popcount(selected);
					const std::size_t exhausted_count = std::popcount(exhausted);
					// The split edge supplies one support incidence to each child.
					// The new side may therefore be support-bivalent when it receives
					// one old group; exact feasibility below then requires a hidden
					// extra copy in its root bundle.  The retained marked-root side
					// also remains at least support-bivalent.
					if ((!root_is_universal
							&& !allow_special_split_for_nonuniversal_root
							&& selected_count == 0)
						|| bundle_count - exhausted_count < 1) {
						if (exhausted == 0) {
							break;
						}
						exhausted = (exhausted - 1) & selected;
						continue;
					}
					const std::size_t positive_selected_count
						= selected_count - exhausted_count;
					const std::size_t retained_support_group_count
						= bundle_count - exhausted_count;
					// Root orientation is determined structurally as well as by exact
					// incidence count.  The diverging side may not have more support
					// groups.  A support tie can continue only when some selected
					// bundle remains at the root too: that is the virtual-double case.
					if (selected_count > retained_support_group_count) {
						if (exhausted == 0) {
							break;
						}
						exhausted = (exhausted - 1) & selected;
						continue;
					}
					const bool support_tie
						= selected_count == retained_support_group_count;
					const bool support_tie_has_virtual_double
						= positive_selected_count != 0;
					// Each selected non-exhausted group needs at least one unit of
					// surplus.  Besides being infeasible, violating this inequality
					// would make the provisional child support larger than E+1.
					if (positive_selected_count > surplus) {
						if (exhausted == 0) {
							break;
						}
						exhausted = (exhausted - 1) & selected;
						continue;
					}

					split_type child;
					for_each_support_edge([&](Int first, Int second) {
						child.set_support_edge(first, second);
					});
					vertex_mask removed = exhausted;
					while (removed != 0) {
						const Int vertex
							= static_cast<Int>(std::countr_zero(removed));
						removed &= removed - 1;
						child.set_root_support(vertex, false);
					}
					child.set_root_support(new_vertex);
					vertex_mask moved_groups = selected;
					while (moved_groups != 0) {
						const Int vertex = static_cast<Int>(
							std::countr_zero(moved_groups)
						);
						moved_groups &= moved_groups - 1;
						child.set_core_support_edge(vertex, new_vertex);
					}
					child.validate();

					classification_mask_type image_mask = 0;
					vertex_mask active_root_mask = 0;
					bool has_simple_image = false;
					std::uint64_t witness_count = 0;
					std::size_t first_tadpoles = 0;
					std::size_t first_separated = 0;
					bool found_witness = false;

					for (std::size_t tadpoles = 0;
						 tadpoles <= surplus;
						 ++tadpoles) {
						const std::size_t root_excess = surplus - tadpoles;
						if (!active_allocation_exists_with_constraints(
								parent_adjacency,
								tadpoles,
								exhausted,
								selected & ~exhausted
							)) {
							continue;
						}

						const std::size_t root_valence
							= bundle_count + surplus + tadpoles;
						for (std::size_t separated = 0;
							 separated <= tadpoles;
							 ++separated) {
							const std::size_t moved_count
								= separated + selected_count;
							if (moved_count > root_valence) {
								continue;
							}
							const std::size_t retained_count
								= root_valence - moved_count;
							if (moved_count < 2 || retained_count < 2
								|| moved_count > retained_count) {
								continue;
							}
							const bool equal_split
								= moved_count == retained_count;

							if (equal_split
								&& !equal_split_has_witness(
									tadpoles,
									root_excess,
									selected,
									exhausted,
									root_neighbors,
									positive_selected_count
								)) {
								continue;
							}

							const classification_mask_type exact_mask
								= exact_split_image_classification_mask(
									parent_adjacency,
									child,
									tadpoles,
									separated,
									selected,
									exhausted
								);
							const bool witness_is_simple
								= child.total_surplus() == 0;
							// Vertex zero is immutable.  A strict split may continue only
							// with that same root; an equal split is terminal even when its
							// simple image is admissible.
							constexpr classification_mask_type active_bits
								= classification_bit(
									transient_classification::active_transient
								) | classification_bit(
									transient_classification::active_admissible
								);
							const bool support_orientation_can_continue
								= !support_tie || support_tie_has_virtual_double;
							const vertex_mask witness_active_roots
								= !equal_split
									&& support_orientation_can_continue
									&& (exact_mask & active_bits) != 0
									? vertex_mask{1}
									: vertex_mask{0};
							if (equal_split && !witness_is_simple) {
								continue;
							}
							if (witness_active_roots == 0 && !witness_is_simple) {
								continue;
							}

							if (!found_witness) {
								found_witness = true;
								first_tadpoles = tadpoles;
								first_separated = separated;
							}
							image_mask = static_cast<classification_mask_type>(
								image_mask | exact_mask
							);
							active_root_mask |= witness_active_roots;
							has_simple_image |= witness_is_simple;
							++witness_count;
						}
					}

					if (found_witness
#if !defined(GC_DISABLE_MAX_CONTRACTION_PREFILTER)
						&& (!apply_max_contraction_prefilter
							|| child.root_contraction_is_maximal(new_vertex))
#endif
					) {
						std::invoke(
							emit,
							symbolic_split_type{
								std::move(child),
								first_tadpoles,
								first_separated,
								selected,
								exhausted,
								image_mask,
								active_root_mask,
								has_simple_image,
								witness_count
							}
						);
						++emitted;
					}

					if (exhausted == 0) {
						break;
					}
					exhausted = (exhausted - 1) & selected;
				}
				if (selected == 0) {
					break;
				}
				selected = (selected - 1) & root_neighbors;
			}
			return emitted;
		}

		template <typename Emit>
		std::uint64_t for_each_symbolic_relevant_root_split(Emit&& emit) const {
			return for_each_symbolic_relevant_root_split_impl(
				std::forward<Emit>(emit), false
			);
		}

		template <typename Emit>
		std::uint64_t for_each_relevant_root_split(Emit&& emit) const {
			return for_each_symbolic_relevant_root_split(
				std::forward<Emit>(emit)
			);
		}

	private:
		template <Int, Int, typename>
		friend class support_transient_graph;

		auto contract_root_neighbor(Int removed_vertex) const {
			static_assert(N_VERTICES > 1);
			static_assert(N_EDGES > 0);
			using parent_type = support_transient_graph<
				static_cast<Int>(N_VERTICES - 1),
				static_cast<Int>(N_EDGES - 1),
				FieldType
			>;
			parent_type parent;
			const auto parent_label = [removed_vertex](Int vertex) {
				return vertex < removed_vertex
					? vertex
					: static_cast<Int>(vertex - 1);
			};
			for_each_support_edge([&](Int first, Int second) {
				if (first == 0 && second == removed_vertex) {
					return;
				}
				if (first == removed_vertex || second == removed_vertex) {
					const Int other = first == removed_vertex ? second : first;
					if (other != 0) {
						parent.set_root_support(parent_label(other));
					}
					return;
				}
				parent.set_support_edge(
					parent_label(first), parent_label(second)
				);
			});
			parent.validate();
			return parent;
		}

		template <typename CandidateChild>
		bool same_support_after_reinserting(
			const CandidateChild& candidate,
			Int removed_vertex
		) const {
			if (candidate.support_edge_count() != support_edge_count()) {
				return false;
			}
			const Int appended_vertex = static_cast<Int>(N_VERTICES - 1);
			const auto original_label = [removed_vertex](
				Int vertex
			) {
				if (vertex == appended_vertex) {
					return removed_vertex;
				}
				return vertex < removed_vertex
					? vertex
					: static_cast<Int>(vertex + 1);
			};
			bool equal = true;
			candidate.for_each_support_edge([&](Int first, Int second) {
				if (!has_support_edge(
						original_label(first), original_label(second)
					)) {
					equal = false;
				}
			});
			return equal;
		}

		template <typename Parent>
		bool contraction_parent_is_eligible(
			const Parent& parent,
			Int removed_vertex
		) const {
			if (!parent.has_repairable_support_valences()
				|| !parent.has_vertex_irreducible_support()
				|| !parent.has_active_allocation()) {
				return false;
			}
			bool regenerates = false;
			parent.for_each_symbolic_relevant_root_split_impl(
				[&](typename Parent::symbolic_split_type split) {
					if (!regenerates
						&& (split.active_root_mask != 0 || split.has_simple_image)
						&& same_support_after_reinserting(
							split.child, removed_vertex
						)) {
						regenerates = true;
					}
				},
				false
			);
			return regenerates;
		}

		classification_mask_type exact_split_image_classification_mask(
			const support_adjacency_array& parent_adjacency,
			const split_type& child,
			std::size_t parent_tadpoles,
			std::size_t separated_tadpoles,
			vertex_mask selected,
			vertex_mask exhausted
		) const noexcept {
			const std::size_t child_tadpoles
				= parent_tadpoles - separated_tadpoles;
			if (child.total_surplus() == 0) {
				// There is only one child allocation.  Its existence together with
				// an active parent was already checked by the (A,D,t) witness test.
				return static_cast<classification_mask_type>(
					child.classification_mask()
					& ~classification_bit(transient_classification::discard)
				);
			}

			const auto child_adjacency = child.support_adjacency();
			if (!child.fixed_active_conditions(
					child_adjacency, child_tadpoles
				)) {
				return 0;
			}

			const std::size_t parent_surplus = total_surplus();
			const std::size_t parent_root_valence
				= std::popcount(parent_adjacency[0])
					+ parent_surplus + parent_tadpoles;
			const bool parent_may_tie
				= allows_tied_root_as_seed_parent(parent_adjacency);
			const std::size_t child_root_valence
				= std::popcount(child_adjacency[0])
					+ child.total_surplus() + child_tadpoles;
			const Int new_vertex = N_VERTICES;

			// The opened tadpoles make the new root bundle have multiplicity
			// separated+1, hence exact excess separated.  It must itself obey the
			// child active bounds.
			const std::size_t new_degree = std::popcount(
				child_adjacency[static_cast<std::size_t>(new_vertex)]
			);
			if (new_degree >= child_root_valence) {
				return 0;
			}
			const std::size_t new_lower = new_degree < 3 ? 3 - new_degree : 0;
			const std::size_t new_upper
				= child_root_valence - new_degree - 1;
			if (separated_tadpoles < new_lower
				|| separated_tadpoles > new_upper) {
				return 0;
			}

			std::size_t lower_sum = 0;
			std::size_t upper_sum = 0;
			vertex_mask parent_neighbors = parent_adjacency[0];
			while (parent_neighbors != 0) {
				const Int vertex = static_cast<Int>(
					std::countr_zero(parent_neighbors)
				);
				parent_neighbors &= parent_neighbors - 1;
				const vertex_mask bit = std::uint64_t{1} << vertex;

				const std::size_t parent_degree = std::popcount(
					parent_adjacency[static_cast<std::size_t>(vertex)]
				);
				std::size_t upper = 0;
				if (!root_bundle_excess_upper_bound(
						parent_degree,
						parent_root_valence,
						parent_may_tie,
						upper
					)) {
					return 0;
				}
				std::size_t lower
					= parent_degree < 3 ? 3 - parent_degree : 0;

				if ((exhausted & bit) != 0) {
					if (lower != 0) {
						return 0;
					}
					continue;
				}
				if ((selected & bit) != 0) {
					lower = std::max<std::size_t>(lower, 1);
				}

				const std::size_t child_degree = std::popcount(
					child_adjacency[static_cast<std::size_t>(vertex)]
				);
				if (child_degree >= child_root_valence) {
					return 0;
				}
				std::size_t child_lower
					= child_degree < 3 ? 3 - child_degree : 0;
				std::size_t child_upper
					= child_root_valence - child_degree - 1;
				if ((selected & bit) != 0) {
					// y_v=x_v-1 for a selected group which remains supported.
					++child_lower;
					++child_upper;
				}
				lower = std::max(lower, child_lower);
				upper = std::min(upper, child_upper);
				if (lower > upper) {
					return 0;
				}
				lower_sum += lower;
				upper_sum += upper;
			}

			const std::size_t parent_root_excess
				= parent_surplus - parent_tadpoles;
			if (parent_root_excess < lower_sum
				|| parent_root_excess > upper_sum) {
				return 0;
			}
			return classification_bit(
				transient_classification::active_transient
			);
		}

		// Decide whether an ACTIVE exact member of this support family realizes
		// the bundle conditions used by a symbolic split.  `zero_excess` means
		// x_v=0 (the selected edge exhausts that group); `positive_excess` means
		// x_v>=1 (the group remains supported after one edge moves).
		bool active_allocation_exists_with_constraints(
			const support_adjacency_array& adjacency,
			std::size_t tadpoles,
			vertex_mask zero_excess,
			vertex_mask positive_excess
		) const noexcept {
			const std::size_t surplus = total_surplus();
			if (tadpoles > surplus || (zero_excess & positive_excess) != 0) {
				return false;
			}
			const vertex_mask root_neighbors = adjacency[0];
			if (((zero_excess | positive_excess) & ~root_neighbors) != 0
				|| !fixed_active_conditions(adjacency, tadpoles)) {
				return false;
			}

			const std::size_t root_valence
				= std::popcount(root_neighbors) + surplus + tadpoles;
			const bool root_may_tie
				= allows_tied_root_as_seed_parent(adjacency);
			std::size_t lower_sum = 0;
			std::size_t upper_sum = 0;
			vertex_mask neighbors = root_neighbors;
			while (neighbors != 0) {
				const Int vertex = static_cast<Int>(std::countr_zero(neighbors));
				neighbors &= neighbors - 1;
				const vertex_mask bit = std::uint64_t{1} << vertex;
				const std::size_t degree = std::popcount(
					adjacency[static_cast<std::size_t>(vertex)]
				);
				std::size_t upper = 0;
				if (!root_bundle_excess_upper_bound(
						degree, root_valence, root_may_tie, upper
					)) {
					return false;
				}
				std::size_t lower = degree < 3 ? 3 - degree : 0;
				if ((zero_excess & bit) != 0) {
					if (lower != 0) {
						return false;
					}
					continue;
				}
				if ((positive_excess & bit) != 0) {
					lower = std::max<std::size_t>(lower, 1);
				}
				if (lower > upper) {
					return false;
				}
				lower_sum += lower;
				upper_sum += upper;
			}
			const std::size_t root_excess = surplus - tadpoles;
			return lower_sum <= root_excess && root_excess <= upper_sum;
		}

		template <typename ExcessAt>
		bool exact_parent_allocation_is_stable(
			const support_adjacency_array& adjacency,
			std::size_t tadpoles,
			ExcessAt&& excess_at
		) const noexcept {
			if (!fixed_active_conditions(adjacency, tadpoles)) {
				return false;
			}
			const std::size_t surplus = total_surplus();
			if (tadpoles > surplus) {
				return false;
			}
			const std::size_t root_valence
				= std::popcount(adjacency[0]) + surplus + tadpoles;
			const bool root_may_tie
				= allows_tied_root_as_seed_parent(adjacency);
			std::size_t excess_sum = 0;
			vertex_mask neighbors = adjacency[0];
			while (neighbors != 0) {
				const Int vertex = static_cast<Int>(std::countr_zero(neighbors));
				neighbors &= neighbors - 1;
				const std::size_t excess = std::invoke(excess_at, vertex);
				excess_sum += excess;
				const std::size_t exact_degree = std::popcount(
					adjacency[static_cast<std::size_t>(vertex)]
				) + excess;
				if (exact_degree > root_valence
					|| (!root_may_tie && exact_degree == root_valence)) {
					return false;
				}
			}
			return excess_sum == surplus - tadpoles;
		}

		template <typename ChildAdjacency>
		static vertex_mask simple_unique_maximum_mask(
			const ChildAdjacency& adjacency
		) noexcept {
			std::size_t maximum = 0;
			Int maximum_vertex = 0;
			std::size_t maximum_count = 0;
			for (Int vertex = 0; vertex <= N_VERTICES; ++vertex) {
				const std::size_t degree = std::popcount(
					adjacency[static_cast<std::size_t>(vertex)]
				);
				if (degree > maximum) {
					maximum = degree;
					maximum_vertex = vertex;
					maximum_count = 1;
				} else if (degree == maximum) {
					++maximum_count;
				}
			}
			return maximum_count == 1
				? (std::uint64_t{1} << maximum_vertex)
				: vertex_mask{0};
		}

		template <typename ChildAdjacency>
		static bool common_bundle_makes_unique_maximum(
			const ChildAdjacency& adjacency,
			Int candidate,
			std::size_t surplus
		) noexcept {
			const std::size_t candidate_valence = std::popcount(
				adjacency[static_cast<std::size_t>(candidate)]
			) + surplus;
			for (Int vertex = 0; vertex <= N_VERTICES; ++vertex) {
				if (vertex == candidate) {
					continue;
				}
				std::size_t valence = std::popcount(
					adjacency[static_cast<std::size_t>(vertex)]
				);
				if (vertex == 0) {
					valence += surplus;
				}
				if (valence >= candidate_valence) {
					return false;
				}
			}
			return true;
		}

		vertex_mask exact_active_root_mask_for_split(
			const support_adjacency_array& parent_adjacency,
			const split_type& child,
			std::size_t parent_tadpoles,
			std::size_t separated_tadpoles,
			vertex_mask selected,
			vertex_mask exhausted,
			classification_mask_type old_root_classifications
		) const noexcept {
			vertex_mask result = 0;
			const classification_mask_type active_bits
				= static_cast<classification_mask_type>(
					classification_bit(transient_classification::active_transient)
					| classification_bit(
						transient_classification::active_admissible
					)
				);
			if ((old_root_classifications & active_bits) != 0) {
				result |= std::uint64_t{1};
			}

			const auto child_adjacency = child.support_adjacency();
			const std::size_t child_surplus = child.total_surplus();
			if (child_surplus == 0) {
				return static_cast<vertex_mask>(
					result | simple_unique_maximum_mask(child_adjacency)
				);
			}

			const vertex_mask positive_selected = selected & ~exhausted;
			const Int new_vertex = N_VERTICES;
			for (Int candidate = 1; candidate <= N_VERTICES; ++candidate) {
				const vertex_mask candidate_bit
					= std::uint64_t{1} << candidate;
				// With positive child surplus, changing the mark is exact only
				// through a surviving common bundle 0--candidate.
				if ((child_adjacency[0] & candidate_bit) == 0
					|| !common_bundle_makes_unique_maximum(
						child_adjacency, candidate, child_surplus
					)) {
					continue;
				}

				bool stable_parent = false;
				if (candidate == new_vertex) {
					// Every remaining defect is the excess on the split edge.
					if (parent_tadpoles == separated_tadpoles
						&& separated_tadpoles == child_surplus) {
						stable_parent = exact_parent_allocation_is_stable(
							parent_adjacency,
							parent_tadpoles,
							[&](Int vertex) -> std::size_t {
								return (positive_selected
									& (std::uint64_t{1} << vertex)) != 0
									? 1
									: 0;
							}
						);
					}
				} else if (parent_tadpoles == 0
					&& separated_tadpoles == 0) {
					// Every remaining defect is excess on the old-root--candidate
					// bundle.  A selected surviving bundle consumed one additional
					// parent excess to create its new support edge.
					stable_parent = exact_parent_allocation_is_stable(
						parent_adjacency,
						0,
						[&](Int vertex) -> std::size_t {
							std::size_t excess = (positive_selected
								& (std::uint64_t{1} << vertex)) != 0
								? 1
								: 0;
							if (vertex == candidate) {
								excess += child_surplus;
							}
							return excess;
						}
					);
				}
				if (stable_parent) {
					result |= candidate_bit;
				}
			}
			return result;
		}

		static bool equal_split_has_witness(
			std::size_t tadpoles,
			std::size_t root_excess,
			vertex_mask selected,
			vertex_mask exhausted,
			vertex_mask root_neighbors,
			std::size_t positive_selected_count
		) noexcept {
			// Exact generator rule: a tied split with any parent tadpole is
			// terminal/discarded.  Simplicity of the retained root side then forces
			// selected non-exhausted groups to have x=1 and every other group x=0.
			if (tadpoles != 0 || root_excess != positive_selected_count) {
				return false;
			}
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				const vertex_mask bit = std::uint64_t{1} << vertex;
				if ((root_neighbors & bit) == 0) {
					continue;
				}
				if ((selected & bit) == 0) {
					// moved count 0 is lexicographically before retained count 1.
					return true;
				}
				if ((exhausted & bit) != 0) {
					// moved count 1 is after retained count 0.
					return false;
				}
				// A selected non-exhausted group has exact size two here, so the
				// two complementary counts are equal and comparison continues.
			}
			return true;
		}

		classification_mask_type classification_mask_for_domain(
			std::size_t tadpoles,
			bool has_fixed_bundle,
			Int fixed_vertex,
			std::size_t fixed_excess
		) const {
			const std::size_t surplus = total_surplus();
			if (tadpoles > surplus) {
				return 0;
			}
			const std::size_t total_excess = surplus - tadpoles;
			const support_adjacency_array adjacency = support_adjacency();
			const vertex_mask root_neighbors = adjacency[0];
			const std::size_t bundle_count = std::popcount(root_neighbors);

			std::size_t variable_excess = total_excess;
			std::size_t variable_count = bundle_count;
			if (has_fixed_bundle) {
				if (fixed_vertex == 0 || fixed_vertex >= N_VERTICES
					|| (root_neighbors
						& (std::uint64_t{1} << fixed_vertex)) == 0
					|| fixed_excess > total_excess) {
					return 0;
				}
				variable_excess -= fixed_excess;
				--variable_count;
			}
			if (variable_count == 0 && variable_excess != 0) {
				return 0;
			}

			if (surplus == 0) {
				return classification_bit(simple_classification(adjacency));
			}

			const bool active = active_allocation_exists(
				adjacency,
				tadpoles,
				has_fixed_bundle,
				fixed_vertex,
				fixed_excess,
				variable_excess,
				variable_count
			);
			classification_mask_type result = active
				? classification_bit(transient_classification::active_transient)
				: classification_mask_type{0};
			if (!active || !every_allocation_is_active(
					adjacency,
					tadpoles,
					has_fixed_bundle,
					fixed_vertex,
					fixed_excess,
					variable_excess,
					variable_count
				)) {
				result = static_cast<classification_mask_type>(
					result | classification_bit(
						transient_classification::discard
					)
				);
			}
			return result;
		}

		static bool support_is_connected(
			const support_adjacency_array& adjacency
		) noexcept {
			vertex_mask visited = std::uint64_t{1};
			vertex_mask frontier = visited;
			while (frontier != 0) {
				const unsigned vertex = std::countr_zero(frontier);
				frontier &= frontier - 1;
				const vertex_mask newly_reached
					= adjacency[vertex] & ~visited;
				visited |= newly_reached;
				frontier |= newly_reached;
			}
			vertex_mask all_vertices;
			if constexpr (N_VERTICES == 64) {
				all_vertices = std::numeric_limits<vertex_mask>::max();
			} else {
				all_vertices = (std::uint64_t{1} << N_VERTICES) - 1;
			}
			return visited == all_vertices;
		}

		static bool support_is_vertex_irreducible(
			const support_adjacency_array& adjacency
		) noexcept {
			if constexpr (N_VERTICES <= 2) {
				return false;
			}
			vertex_mask all_vertices;
			if constexpr (N_VERTICES == 64) {
				all_vertices = std::numeric_limits<vertex_mask>::max();
			} else {
				all_vertices = (std::uint64_t{1} << N_VERTICES) - 1;
			}
			for (Int removed = 0; removed < N_VERTICES; ++removed) {
				const vertex_mask remaining
					= all_vertices & ~(std::uint64_t{1} << removed);
				vertex_mask visited
					= std::uint64_t{1} << std::countr_zero(remaining);
				vertex_mask frontier = visited;
				while (frontier != 0) {
					const unsigned vertex = std::countr_zero(frontier);
					frontier &= frontier - 1;
					const vertex_mask newly_reached
						= adjacency[vertex] & remaining & ~visited;
					visited |= newly_reached;
					frontier |= newly_reached;
				}
				if (visited != remaining) {
					return false;
				}
			}
			return true;
		}

		// Continuing exact states require a strict root maximum.  Equal exact
		// split sides are handled separately as terminal images; there is no
		// seed-specific rerooting or tied-root exception.
		bool allows_tied_root_as_seed_parent(
			const support_adjacency_array& /* adjacency */
		) const noexcept {
			return false;
		}

		static bool root_bundle_excess_upper_bound(
			std::size_t support_degree,
			std::size_t root_valence,
			bool root_may_tie,
			std::size_t& upper
		) noexcept {
			if (support_degree > root_valence
				|| (!root_may_tie && support_degree == root_valence)) {
				return false;
			}
			upper = root_valence - support_degree;
			if (!root_may_tie) {
				--upper;
			}
			return true;
		}

		static transient_classification simple_classification(
			const support_adjacency_array& adjacency
		) noexcept {
			if (!support_is_connected(adjacency)) {
				return transient_classification::discard;
			}
			const std::size_t root_degree = std::popcount(adjacency[0]);
			std::size_t maximum_degree = root_degree;
			std::size_t maximum_count = 1;
			if (root_degree < 3) {
				return transient_classification::discard;
			}
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				const std::size_t degree = std::popcount(
					adjacency[static_cast<std::size_t>(vertex)]
				);
				if (degree < 3) {
					return transient_classification::discard;
				}
				if (degree > maximum_degree) {
					maximum_degree = degree;
					maximum_count = 1;
				} else if (degree == maximum_degree) {
					++maximum_count;
				}
			}
			if (root_degree != maximum_degree) {
				return transient_classification::discard;
			}
			return maximum_count == 1
				? transient_classification::active_admissible
				: transient_classification::terminal_admissible;
		}

		bool fixed_active_conditions(
			const support_adjacency_array& adjacency,
			std::size_t tadpoles
		) const noexcept {
			if (!support_is_connected(adjacency)) {
				return false;
			}
			const std::size_t root_support_degree = std::popcount(adjacency[0]);
			if (root_support_degree < 1) {
				return false;
			}
			const std::size_t root_valence
				= root_support_degree + total_surplus() + tadpoles;
			if (root_valence < 3) {
				return false;
			}
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				const std::size_t degree = std::popcount(
					adjacency[static_cast<std::size_t>(vertex)]
				);
				if (degree == 0) {
					return false;
				}
				if ((adjacency[0] & (std::uint64_t{1} << vertex)) == 0) {
					if (degree < 3 || degree >= root_valence) {
						return false;
					}
				}
			}
			return true;
		}

		bool active_allocation_exists(
			const support_adjacency_array& adjacency,
			std::size_t tadpoles,
			bool has_fixed_bundle,
			Int fixed_vertex,
			std::size_t fixed_excess,
			std::size_t variable_excess,
			std::size_t /* variable_count */
		) const noexcept {
			if (!fixed_active_conditions(adjacency, tadpoles)) {
				return false;
			}
			const std::size_t root_degree = std::popcount(adjacency[0]);
			const std::size_t root_valence
				= root_degree + total_surplus() + tadpoles;
			const bool root_may_tie
				= allows_tied_root_as_seed_parent(adjacency);
			std::size_t lower_sum = 0;
			std::size_t upper_sum = 0;
			vertex_mask neighbors = adjacency[0];
			while (neighbors != 0) {
				const Int vertex = static_cast<Int>(std::countr_zero(neighbors));
				neighbors &= neighbors - 1;
				const std::size_t degree = std::popcount(
					adjacency[static_cast<std::size_t>(vertex)]
				);
				const std::size_t lower = degree < 3 ? 3 - degree : 0;
				std::size_t upper = 0;
				if (!root_bundle_excess_upper_bound(
						degree, root_valence, root_may_tie, upper
					)) {
					return false;
				}
				if (has_fixed_bundle && vertex == fixed_vertex) {
					if (fixed_excess < lower || fixed_excess > upper) {
						return false;
					}
				} else {
					lower_sum += lower;
					upper_sum += upper;
				}
			}
			return lower_sum <= variable_excess
				&& variable_excess <= upper_sum;
		}

		bool every_allocation_is_active(
			const support_adjacency_array& adjacency,
			std::size_t tadpoles,
			bool has_fixed_bundle,
			Int fixed_vertex,
			std::size_t fixed_excess,
			std::size_t variable_excess,
			std::size_t variable_count
		) const noexcept {
			if (!active_allocation_exists(
					adjacency, tadpoles, has_fixed_bundle, fixed_vertex,
					fixed_excess, variable_excess, variable_count
				)) {
				return false;
			}
			if (variable_count <= 1) {
				return true;
			}

			const std::size_t root_degree = std::popcount(adjacency[0]);
			const std::size_t root_valence
				= root_degree + total_surplus() + tadpoles;
			const bool root_may_tie
				= allows_tied_root_as_seed_parent(adjacency);
			vertex_mask neighbors = adjacency[0];
			while (neighbors != 0) {
				const Int vertex = static_cast<Int>(std::countr_zero(neighbors));
				neighbors &= neighbors - 1;
				if (has_fixed_bundle && vertex == fixed_vertex) {
					continue;
				}
				const std::size_t degree = std::popcount(
					adjacency[static_cast<std::size_t>(vertex)]
				);
				const std::size_t lower = degree < 3 ? 3 - degree : 0;
				std::size_t upper = 0;
				if (!root_bundle_excess_upper_bound(
						degree, root_valence, root_may_tie, upper
					)) {
					return false;
				}
				if (lower != 0 || upper < variable_excess) {
					return false;
				}
			}
			return true;
		}

		support_transient_graph permuted_nonroot_unchecked(
			const nonroot_permutation& old_to_new_local
		) const noexcept {
			support_transient_graph result;
			for_each_support_edge([&](Int old_first, Int old_second) {
				const Int new_first = old_first == 0
					? 0
					: static_cast<Int>(old_to_new_local[
						static_cast<std::size_t>(old_first - 1)
					] + 1);
				const Int new_second = old_second == 0
					? 0
					: static_cast<Int>(old_to_new_local[
						static_cast<std::size_t>(old_second - 1)
					] + 1);
				result.set_support_edge_unchecked(new_first, new_second, true);
			});
			return result;
		}

		static void validate_vertex(Int vertex) {
			if (vertex >= N_VERTICES) {
				throw std::out_of_range("support transient vertex");
			}
		}

		static std::size_t support_bit_index(Int first, Int second) {
			validate_vertex(first);
			validate_vertex(second);
			if (first == second) {
				throw std::invalid_argument("simple support has no diagonal edge");
			}
			if (second < first) {
				std::swap(first, second);
			}
			return support_bit_index_unchecked(first, second);
		}

		static constexpr std::size_t support_bit_index_unchecked(
			Int first, Int second
		) noexcept {
			const std::size_t u = static_cast<std::size_t>(first);
			const std::size_t v = static_cast<std::size_t>(second);
			return u * (2 * static_cast<std::size_t>(N_VERTICES) - u - 1) / 2
				+ (v - u - 1);
		}

		void set_support_edge_unchecked(
			Int first, Int second, bool present
		) {
			const std::size_t bit = support_bit_index(first, second);
			const std::uint64_t mask = std::uint64_t{1} << (bit % 64);
			if (present) {
				support_words_[bit / 64] |= mask;
			} else {
				support_words_[bit / 64] &= ~mask;
			}
		}

		static void validate_nonroot_permutation(
			const nonroot_permutation& old_to_new_local
		) {
			std::array<bool, N_NONROOT_VERTICES> seen{};
			for (const Int image : old_to_new_local) {
				const std::size_t index = static_cast<std::size_t>(image);
				if (index >= N_NONROOT_VERTICES || seen[index]) {
					throw std::invalid_argument(
						"invalid support nonroot permutation"
					);
				}
				seen[index] = true;
			}
		}

		support_word_array support_words_{};
};

} // namespace GraphGeneration

namespace std {

template <Int N_VERTICES, Int N_EDGES, typename FieldType>
struct hash<GraphGeneration::support_transient_graph<
	N_VERTICES, N_EDGES, FieldType
>> {
	std::size_t operator()(
		const GraphGeneration::support_transient_graph<
			N_VERTICES, N_EDGES, FieldType
		>& graph
	) const noexcept {
		return graph.hash_value();
	}
};

} // namespace std
