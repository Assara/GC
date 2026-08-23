#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <span>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "CombinatorialUtils.hpp"
#include "TransientGraph.hpp"

namespace GraphGeneration {

template <Int N_VERTICES, Int N_EDGES, typename FieldType>
class rooted_transient_standardizer;

// Compact representation of the transient graphs used by the generator.
//
// Full vertex label 0 is always the distinguished root.  All graph defects
// which are allowed in a transient graph are represented separately:
//
//   * tadpole_count_ is the number of loops at 0;
//   * root_multiplicities_[v - 1] is the exact number of edges 0--v;
//   * core_words_ is the packed upper triangle of the SIMPLE graph induced by
//     vertices 1,...,N_VERTICES-1.
//
// Consequently, this type cannot represent a tadpole or a parallel-edge
// bundle away from the root.  That is intentional: such a graph is discarded
// by the transient generation policy.
template <Int N_VERTICES, Int N_EDGES, typename FieldType = fieldType>
class rooted_transient_graph {
	public:
		static_assert(N_VERTICES > 0);

		using field_type = FieldType;
		using transient_type = transient_graph<N_VERTICES, N_EDGES, FieldType>;
		using hairy_graph_type = typename transient_type::graph_type;
		using hairless_graph_type = typename transient_type::hairless_graph_type;
		using split_type = rooted_transient_graph<
			static_cast<Int>(N_VERTICES + 1),
			static_cast<Int>(N_EDGES + 1),
			FieldType
		>;

		static constexpr Int N_VERTICES_ = N_VERTICES;
		static constexpr Int N_EDGES_ = N_EDGES;
		static constexpr std::size_t N_NONROOT_VERTICES
			= static_cast<std::size_t>(N_VERTICES - 1);
		static constexpr std::size_t CORE_BIT_COUNT
			= N_NONROOT_VERTICES < 2
				? 0
				: N_NONROOT_VERTICES * (N_NONROOT_VERTICES - 1) / 2;
		static constexpr std::size_t CORE_WORD_COUNT = (CORE_BIT_COUNT + 63) / 64;

		using root_multiplicity_array = std::array<Int, N_NONROOT_VERTICES>;
		using core_word_array = std::array<std::uint64_t, CORE_WORD_COUNT>;
		// old_to_new_local[i] gives the new local label of old full vertex i+1.
		using nonroot_permutation = std::array<Int, N_NONROOT_VERTICES>;

		struct analysis_type {
			std::array<Int, N_VERTICES> valences{};
			std::array<Int, N_VERTICES> reduced_valences{};
			Int maximum_valence = 0;
			Int maximum_reduced_valence = 0;
			Int maximum_valence_count = 0;
			bool connected = true;
			bool at_least_trivalent = true;
			bool simple = true;
			// This is true by construction for every structurally valid instance.
			bool defects_at_root = true;
			bool root_is_maximum = false;
			bool root_is_unique_maximum = false;
			bool root_is_reduced_maximum = false;
		};

		constexpr rooted_transient_graph() noexcept = default;

		static rooted_transient_graph rose() {
			static_assert(N_VERTICES == 1);
			rooted_transient_graph result;
			result.tadpole_count_ = N_EDGES;
			return result;
		}

		rooted_transient_graph(
			Int tadpole_count,
			root_multiplicity_array root_multiplicities,
			core_word_array core_words = {}
		)
			: tadpole_count_(tadpole_count),
			  root_multiplicities_(std::move(root_multiplicities)),
			  core_words_(std::move(core_words)) {
			validate();
		}

		explicit rooted_transient_graph(const transient_type& source)
			: rooted_transient_graph(from_transient(source)) {}

		// Relabel the marked vertex to 0 and preserve the relative order of all
		// other labels.  This conversion is label preserving modulo that single
		// deterministic root move; it is not canonicalization.
		static rooted_transient_graph from_transient(const transient_type& source) {
			const Int old_root = source.root();
			if (old_root >= N_VERTICES) {
				throw std::invalid_argument("compact transient source has invalid root");
			}

			std::array<Int, N_VERTICES> old_to_new{};
			old_to_new[static_cast<std::size_t>(old_root)] = 0;
			Int next = 1;
			for (Int old_vertex = 0; old_vertex < N_VERTICES; ++old_vertex) {
				if (old_vertex != old_root) {
					old_to_new[static_cast<std::size_t>(old_vertex)] = next++;
				}
			}

			rooted_transient_graph result;
			for (Int edge = 0; edge < N_EDGES; ++edge) {
				auto [old_u, old_v] = source.graph().getEdge(edge);
				const Int u = old_to_new[static_cast<std::size_t>(old_u)];
				const Int v = old_to_new[static_cast<std::size_t>(old_v)];

				if (u == 0 && v == 0) {
					increment_checked(result.tadpole_count_, "too many root tadpoles");
				} else if (u == 0 || v == 0) {
					const Int other = (u == 0) ? v : u;
					increment_checked(
						result.root_multiplicities_[nonroot_index(other)],
						"root-edge multiplicity overflow"
					);
				} else {
					if (u == v) {
						throw std::invalid_argument(
							"compact transient cannot represent a nonroot tadpole"
						);
					}
					if (result.has_core_edge(u, v)) {
						throw std::invalid_argument(
							"compact transient cannot represent parallel core edges"
						);
					}
					result.set_core_edge_unchecked(u, v, true);
				}
			}
			result.validate();
			return result;
		}

		constexpr Int root() const noexcept { return 0; }

		constexpr Int tadpole_count() const noexcept { return tadpole_count_; }

		constexpr const root_multiplicity_array& root_multiplicities() const noexcept {
			return root_multiplicities_;
		}

		constexpr const core_word_array& core_words() const noexcept {
			return core_words_;
		}

		Int root_multiplicity(Int full_vertex) const {
			return root_multiplicities_[nonroot_index(full_vertex)];
		}

		bool has_core_edge(Int full_u, Int full_v) const {
			const std::size_t bit = core_bit_index(full_u, full_v);
			return (core_words_[bit / 64] & (std::uint64_t{1} << (bit % 64))) != 0;
		}

		// These mutators permit efficient decoding and fixture construction.  A
		// sequence of mutations may temporarily have the wrong edge total; call
		// validate() before using the value as a graph or a hash key.
		constexpr void set_tadpole_count(Int count) noexcept { tadpole_count_ = count; }

		void set_root_multiplicity(Int full_vertex, Int multiplicity) {
			root_multiplicities_[nonroot_index(full_vertex)] = multiplicity;
		}

		void set_core_edge(Int full_u, Int full_v, bool present = true) {
			set_core_edge_unchecked(full_u, full_v, present);
		}

		std::size_t represented_edge_count() const noexcept {
			std::size_t count = tadpole_count_;
			for (const Int multiplicity : root_multiplicities_) {
				count += multiplicity;
			}
			for (const std::uint64_t word : core_words_) {
				count += std::popcount(word);
			}
			return count;
		}

		void validate() const {
			if constexpr (CORE_WORD_COUNT != 0 && CORE_BIT_COUNT % 64 != 0) {
				constexpr std::uint64_t used_mask
					= (std::uint64_t{1} << (CORE_BIT_COUNT % 64)) - 1;
				if ((core_words_.back() & ~used_mask) != 0) {
					throw std::invalid_argument("bits outside compact core triangle");
				}
			}
			if (represented_edge_count() != static_cast<std::size_t>(N_EDGES)) {
				throw std::invalid_argument("compact transient edge count mismatch");
			}
		}

		bool structurally_valid() const noexcept {
			if constexpr (CORE_WORD_COUNT != 0 && CORE_BIT_COUNT % 64 != 0) {
				constexpr std::uint64_t used_mask
					= (std::uint64_t{1} << (CORE_BIT_COUNT % 64)) - 1;
				if ((core_words_.back() & ~used_mask) != 0) {
					return false;
				}
			}
			return represented_edge_count() == static_cast<std::size_t>(N_EDGES);
		}

		hairless_graph_type to_hairless_graph() const {
			validate();
			hairless_graph_type result;
			write_edges(result, 0);
			return result;
		}

		hairless_graph_type without_hair() const { return to_hairless_graph(); }
		hairless_graph_type forget_root() const { return to_hairless_graph(); }

		hairy_graph_type to_hairy_graph() const {
			validate();
			hairy_graph_type result;
			result.half_edges[0] = 0;
			write_edges(result, hairy_graph_type::N_HAIR);
			return result;
		}

		transient_type to_transient() const {
			return transient_type(to_hairy_graph());
		}

		transient_type to_legacy_transient() const {
			return to_transient();
		}

		// Apply a permutation only to full vertices 1,...,N_VERTICES-1.  Vertex
		// 0 remains fixed.  The input values are local labels in 0,...,V-2.
		rooted_transient_graph permuted_nonroot(
			const nonroot_permutation& old_to_new_local
		) const {
			validate_nonroot_permutation(old_to_new_local);
			return permuted_nonroot_unchecked(old_to_new_local);
		}

		template <typename Emit>
		void for_each_core_edge(Emit&& emit) const {
			for (std::size_t word_index = 0;
				 word_index < CORE_WORD_COUNT;
				 ++word_index) {
				std::uint64_t bits = core_words_[word_index];
				while (bits != 0) {
					const std::size_t bit = 64 * word_index
						+ static_cast<std::size_t>(std::countr_zero(bits));
					std::invoke(
						emit,
						CORE_FIRST_VERTICES[bit],
						CORE_SECOND_VERTICES[bit]
					);
					bits &= bits - 1;
				}
			}
		}

		std::array<Int, N_VERTICES> edge_valence_array() const noexcept {
			std::array<Int, N_VERTICES> result{};
			result[0] = static_cast<Int>(2 * tadpole_count_);
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				const Int multiplicity
					= root_multiplicities_[nonroot_index(vertex)];
				result[0] = static_cast<Int>(result[0] + multiplicity);
				result[static_cast<std::size_t>(vertex)] = multiplicity;
			}
			for_each_core_edge([&](Int u, Int v) {
				++result[static_cast<std::size_t>(u)];
				++result[static_cast<std::size_t>(v)];
			});
			return result;
		}

		std::array<Int, N_VERTICES> ordinary_valences() const noexcept {
			return edge_valence_array();
		}

		analysis_type analyze() const {
			validate();
			analysis_type result;
			result.valences = edge_valence_array();
			result.simple = tadpole_count_ == 0;
			result.reduced_valences[0] = 0;

			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				const Int multiplicity
					= root_multiplicities_[nonroot_index(vertex)];
				result.simple &= multiplicity <= 1;
				if (multiplicity != 0) {
					++result.reduced_valences[0];
					++result.reduced_valences[static_cast<std::size_t>(vertex)];
				}
			}
			for_each_core_edge([&](Int u, Int v) {
				++result.reduced_valences[static_cast<std::size_t>(u)];
				++result.reduced_valences[static_cast<std::size_t>(v)];
			});

			for (Int vertex = 0; vertex < N_VERTICES; ++vertex) {
				const std::size_t index = static_cast<std::size_t>(vertex);
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

			std::array<bool, N_VERTICES> visited{};
			std::array<Int, N_VERTICES> stack{};
			std::size_t stack_size = 1;
			std::size_t visited_count = 1;
			visited[0] = true;
			stack[0] = 0;
			while (stack_size != 0) {
				const Int current = stack[--stack_size];
				for (Int candidate = 0; candidate < N_VERTICES; ++candidate) {
					if (candidate == current
						|| visited[static_cast<std::size_t>(candidate)]) {
						continue;
					}
					const bool adjacent = current == 0
						? root_multiplicities_[nonroot_index(candidate)] != 0
						: candidate == 0
							? root_multiplicities_[nonroot_index(current)] != 0
							: has_core_edge(current, candidate);
					if (adjacent) {
						visited[static_cast<std::size_t>(candidate)] = true;
						stack[stack_size++] = candidate;
						++visited_count;
					}
				}
			}
			result.connected = visited_count == static_cast<std::size_t>(N_VERTICES);
			result.root_is_maximum = result.valences[0] == result.maximum_valence;
			result.root_is_unique_maximum
				= result.root_is_maximum && result.maximum_valence_count == 1;
			result.root_is_reduced_maximum
				= result.reduced_valences[0] == result.maximum_reduced_valence;
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
			return state.simple
				? transient_classification::active_admissible
				: transient_classification::active_transient;
		}

		bool can_reduce_to_admissible() const {
			return classification() != transient_classification::discard;
		}

		bool is_admissible_after_removing_hair() const {
			const analysis_type state = analyze();
			return state.simple && state.connected && state.at_least_trivalent;
		}

		// Enumerate only the grouped splits relevant to generation.  A split is
		// described by the number of opened root tadpoles and by a set of
		// root--vertex bundles.  At most one edge is removed from each selected
		// bundle, so the simple nonroot core invariant is preserved directly.
		template <typename Emit>
		std::uint64_t for_each_relevant_root_split(Emit&& emit) const {
			const analysis_type parent = analyze();
			const std::size_t valence = parent.valences[0];
			if (valence < 4) {
				return 0;
			}

			std::vector<Int> bundle_vertices;
			std::vector<Int> bundle_sizes;
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				const Int multiplicity = root_multiplicity(vertex);
				if (multiplicity != 0) {
					bundle_vertices.push_back(vertex);
					bundle_sizes.push_back(multiplicity);
				}
			}
			const std::vector<Int> selection_bounds(bundle_vertices.size(), Int{1});

			Int maximum_untouched_valence = 0;
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				maximum_untouched_valence = std::max(
					maximum_untouched_valence,
					parent.valences[static_cast<std::size_t>(vertex)]
				);
			}

			std::uint64_t emitted = 0;
			for (Int separated_tadpoles = 0;
				 separated_tadpoles <= tadpole_count_;
				 ++separated_tadpoles) {
				combutils::for_each_bounded_count_vector(
					selection_bounds,
					[&](std::span<const Int> selected) {
						std::size_t moved_count = separated_tadpoles;
						for (const Int value : selected) {
							moved_count += value;
						}
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
							if (tadpole_count_ != 0) {
								return;
							}
							for (std::size_t i = 0; i < selected.size(); ++i) {
								if (bundle_sizes[i] - selected[i] > 1) {
									return;
								}
							}
							if (!combutils::is_first_grouped_bipartition(
									selected, bundle_sizes
								)) {
								return;
							}
						}

						split_type child;
						child.set_tadpole_count(static_cast<Int>(
							tadpole_count_ - separated_tadpoles
						));
						for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
							child.set_root_multiplicity(
								vertex, root_multiplicity(vertex)
							);
						}
						for_each_core_edge([&](Int first, Int second) {
							child.set_core_edge(first, second);
						});

						const Int new_vertex = N_VERTICES;
						child.set_root_multiplicity(
							new_vertex,
							static_cast<Int>(separated_tadpoles + 1)
						);
						for (std::size_t i = 0; i < selected.size(); ++i) {
							if (selected[i] == 0) {
								continue;
							}
							const Int vertex = bundle_vertices[i];
							child.set_root_multiplicity(
								vertex,
								static_cast<Int>(child.root_multiplicity(vertex) - 1)
							);
							child.set_core_edge(vertex, new_vertex);
						}

						const transient_classification kind = child.classification();
						if (kind != transient_classification::discard) {
							std::invoke(emit, std::move(child), kind);
							++emitted;
						}
					}
				);
			}
			return emitted;
		}

		signedInt compare(const rooted_transient_graph& other) const noexcept {
			if (tadpole_count_ != other.tadpole_count_) {
				return tadpole_count_ < other.tadpole_count_ ? -1 : 1;
			}
			for (std::size_t i = 0; i < N_NONROOT_VERTICES; ++i) {
				if (root_multiplicities_[i] != other.root_multiplicities_[i]) {
					return root_multiplicities_[i] < other.root_multiplicities_[i]
						? -1 : 1;
				}
			}
			for (std::size_t bit = 0; bit < CORE_BIT_COUNT; ++bit) {
				const bool ours
					= (core_words_[bit / 64] & (std::uint64_t{1} << (bit % 64))) != 0;
				const bool theirs
					= (other.core_words_[bit / 64]
						& (std::uint64_t{1} << (bit % 64))) != 0;
				if (ours != theirs) {
					return ours ? 1 : -1;
				}
			}
			return 0;
		}

		bool operator==(const rooted_transient_graph&) const noexcept = default;
		bool operator<(const rooted_transient_graph& other) const noexcept {
			return compare(other) < 0;
		}

		std::size_t hash_value() const noexcept {
			std::uint64_t hash = 0x9e3779b97f4a7c15ULL;
			auto mix = [&](std::uint64_t value) {
				value += 0x9e3779b97f4a7c15ULL;
				value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ULL;
				value = (value ^ (value >> 27)) * 0x94d049bb133111ebULL;
				value ^= value >> 31;
				hash ^= value + 0x9e3779b97f4a7c15ULL + (hash << 6) + (hash >> 2);
			};
			mix(tadpole_count_);
			for (const Int multiplicity : root_multiplicities_) {
				mix(multiplicity);
			}
			for (const std::uint64_t word : core_words_) {
				mix(word);
			}
			return static_cast<std::size_t>(hash);
		}

	private:
		template <Int, Int, typename>
		friend class rooted_transient_standardizer;

		using core_endpoint_array = std::array<Int, CORE_BIT_COUNT>;

		static consteval core_endpoint_array make_core_first_vertices() {
			core_endpoint_array result{};
			std::size_t bit = 0;
			for (Int first = 1; first < N_VERTICES; ++first) {
				for (Int second = static_cast<Int>(first + 1);
					 second < N_VERTICES;
					 ++second) {
					result[bit++] = first;
				}
			}
			return result;
		}

		static consteval core_endpoint_array make_core_second_vertices() {
			core_endpoint_array result{};
			std::size_t bit = 0;
			for (Int first = 1; first < N_VERTICES; ++first) {
				for (Int second = static_cast<Int>(first + 1);
					 second < N_VERTICES;
					 ++second) {
					result[bit++] = second;
				}
			}
			return result;
		}

		inline static constexpr core_endpoint_array CORE_FIRST_VERTICES
			= make_core_first_vertices();
		inline static constexpr core_endpoint_array CORE_SECOND_VERTICES
			= make_core_second_vertices();

		rooted_transient_graph permuted_nonroot_unchecked(
			const nonroot_permutation& old_to_new_local
		) const noexcept {
			rooted_transient_graph result;
			result.tadpole_count_ = tadpole_count_;
			for (std::size_t old_local = 0;
				 old_local < N_NONROOT_VERTICES;
				 ++old_local) {
				const std::size_t new_local
					= static_cast<std::size_t>(old_to_new_local[old_local]);
				result.root_multiplicities_[new_local]
					= root_multiplicities_[old_local];
			}

			for (std::size_t word_index = 0;
				 word_index < CORE_WORD_COUNT;
				 ++word_index) {
				std::uint64_t bits = core_words_[word_index];
				while (bits != 0) {
					const std::size_t old_bit = 64 * word_index
						+ static_cast<std::size_t>(std::countr_zero(bits));
					std::size_t first = old_to_new_local[
						static_cast<std::size_t>(CORE_FIRST_VERTICES[old_bit] - 1)
					];
					std::size_t second = old_to_new_local[
						static_cast<std::size_t>(CORE_SECOND_VERTICES[old_bit] - 1)
					];
					if (second < first) {
						std::swap(first, second);
					}
					const std::size_t new_bit
						= first * (2 * N_NONROOT_VERTICES - first - 1) / 2
							+ (second - first - 1);
					result.core_words_[new_bit / 64]
						|= std::uint64_t{1} << (new_bit % 64);
					bits &= bits - 1;
				}
			}
			return result;
		}

		static constexpr std::size_t nonroot_index(Int full_vertex) {
			if (full_vertex == 0 || full_vertex >= N_VERTICES) {
				throw std::out_of_range("compact transient nonroot vertex");
			}
			return static_cast<std::size_t>(full_vertex - 1);
		}

		static constexpr std::size_t core_bit_index(Int full_u, Int full_v) {
			std::size_t u = nonroot_index(full_u);
			std::size_t v = nonroot_index(full_v);
			if (u == v) {
				throw std::invalid_argument("compact core has no diagonal bits");
			}
			if (u > v) {
				std::swap(u, v);
			}
			return u * (2 * N_NONROOT_VERTICES - u - 1) / 2 + (v - u - 1);
		}

		void set_core_edge_unchecked(Int full_u, Int full_v, bool present) {
			const std::size_t bit = core_bit_index(full_u, full_v);
			const std::uint64_t mask = std::uint64_t{1} << (bit % 64);
			if (present) {
				core_words_[bit / 64] |= mask;
			} else {
				core_words_[bit / 64] &= ~mask;
			}
		}

		static void validate_nonroot_permutation(
			const nonroot_permutation& old_to_new_local
		) {
			std::array<bool, N_NONROOT_VERTICES> seen{};
			for (const Int image : old_to_new_local) {
				const std::size_t index = static_cast<std::size_t>(image);
				if (index >= N_NONROOT_VERTICES || seen[index]) {
					throw std::invalid_argument("invalid compact nonroot permutation");
				}
				seen[index] = true;
			}
		}

		static void increment_checked(Int& value, const char* message) {
			if (value == std::numeric_limits<Int>::max()) {
				throw std::overflow_error(message);
			}
			++value;
		}

		template <typename GraphType>
		void write_edges(GraphType& result, Int first_edge_position) const {
			Int position = first_edge_position;
			auto write = [&](Int u, Int v) {
				result.half_edges[position++] = u;
				result.half_edges[position++] = v;
			};

			for (Int i = 0; i < tadpole_count_; ++i) {
				write(0, 0);
			}
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				for (Int copy = 0; copy < root_multiplicity(vertex); ++copy) {
					write(0, vertex);
				}
			}
			for_each_core_edge([&](Int u, Int v) { write(u, v); });
		}

		Int tadpole_count_ = 0;
		root_multiplicity_array root_multiplicities_{};
		core_word_array core_words_{};
};

} // namespace GraphGeneration

namespace std {

template <Int N_VERTICES, Int N_EDGES, typename FieldType>
struct hash<GraphGeneration::rooted_transient_graph<
	N_VERTICES, N_EDGES, FieldType
>> {
	std::size_t operator()(
		const GraphGeneration::rooted_transient_graph<
			N_VERTICES, N_EDGES, FieldType
		>& graph
	) const noexcept {
		return graph.hash_value();
	}
};

} // namespace std
