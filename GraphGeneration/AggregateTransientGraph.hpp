#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "GraphGeneration/RootedTransientGraph.hpp"

namespace GraphGeneration {

// Experimental quotient of rooted_transient_graph.
//
// An instance does not denote one multigraph.  It denotes every distribution
// of a fixed number of excess root edges over the displayed root-neighbour
// groups.  Its stored data are only
//
//   * the simple support graph (with the distinguished root fixed at 0), and
//   * the number of root tadpoles.
//
// At fixed N_EDGES the total excess is derived as
//
//   N_EDGES - tadpoles - support_edges.
//
// This deliberately forgets which root-neighbour gets which excess.  Split
// generation below therefore has existential/family semantics: it emits the
// union of all split cases admitted by at least one represented allocation.
// Merging child families can forget correlations from an earlier split.  The
// type is kept separate from the exact production representation so that this
// quotient can be tested without changing the working generator.
template <Int N_VERTICES, Int N_EDGES, typename FieldType = fieldType>
class aggregate_transient_graph {
	public:
		static_assert(N_VERTICES > 0);
		static_assert(N_VERTICES <= 64,
			"aggregate split descriptors use one 64-bit vertex mask");

		using field_type = FieldType;
		using exact_type = rooted_transient_graph<N_VERTICES, N_EDGES, FieldType>;
		using hairless_graph_type = Graph<
			N_VERTICES, N_EDGES, 0, 0, 0, 0, FieldType
		>;
		using split_type = aggregate_transient_graph<
			static_cast<Int>(N_VERTICES + 1),
			static_cast<Int>(N_EDGES + 1),
			FieldType
		>;

		static constexpr Int N_VERTICES_ = N_VERTICES;
		static constexpr Int N_EDGES_ = N_EDGES;
		static constexpr std::size_t SUPPORT_BIT_COUNT
			= static_cast<std::size_t>(N_VERTICES)
				* static_cast<std::size_t>(N_VERTICES - 1) / 2;
		static constexpr std::size_t SUPPORT_WORD_COUNT
			= (SUPPORT_BIT_COUNT + 63) / 64;
		static constexpr std::size_t N_NONROOT_VERTICES
			= static_cast<std::size_t>(N_VERTICES - 1);

		using support_word_array = std::array<std::uint64_t, SUPPORT_WORD_COUNT>;
		// old_to_new_local[i] is the new local label of old full vertex i+1.
		using nonroot_permutation = std::array<Int, N_NONROOT_VERTICES>;
		using permutation = nonroot_permutation;
		using excess_array = std::array<Int, N_NONROOT_VERTICES>;

		struct symbolic_split_descriptor {
			Int separated_tadpoles = 0;
			// Bits are indexed by full vertex labels.  Bit zero is unused.
			std::uint64_t selected_groups = 0;
			std::uint64_t exhausted_groups = 0;
		};

		struct analysis_type {
			std::array<Int, N_VERTICES> support_degrees{};
			Int root_full_valence = 0;
			Int maximum_support_degree = 0;
			Int excess_edges = 0;
			bool structurally_valid = false;
			bool connected = false;
			bool simple = false;
			bool root_is_reduced_maximum = false;
			bool has_active_allocation = false;
		};

		constexpr aggregate_transient_graph() noexcept = default;

		static aggregate_transient_graph rose() {
			static_assert(N_VERTICES == 1);
			aggregate_transient_graph result;
			result.tadpole_count_ = N_EDGES;
			return result;
		}

		aggregate_transient_graph(Int tadpoles, support_word_array support = {})
			: tadpole_count_(tadpoles), support_words_(std::move(support)) {
			validate();
		}

		explicit aggregate_transient_graph(const exact_type& exact)
			: aggregate_transient_graph(from_exact(exact)) {}

		static aggregate_transient_graph from_exact(const exact_type& exact) {
			exact.validate();
			aggregate_transient_graph result;
			result.tadpole_count_ = exact.tadpole_count();
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				if (exact.root_multiplicity(vertex) != 0) {
					result.set_support_edge_unchecked(0, vertex, true);
				}
			}
			exact.for_each_core_edge([&](Int first, Int second) {
				result.set_support_edge_unchecked(first, second, true);
			});
			result.validate();
			return result;
		}

		constexpr Int root() const noexcept { return 0; }
		constexpr Int tadpole_count() const noexcept { return tadpole_count_; }
		constexpr const support_word_array& support_words() const noexcept {
			return support_words_;
		}

		constexpr void set_tadpole_count(Int count) noexcept {
			tadpole_count_ = count;
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

		std::size_t support_edge_count() const noexcept {
			std::size_t result = 0;
			for (const std::uint64_t word : support_words_) {
				result += std::popcount(word);
			}
			return result;
		}

		Int excess_edge_count() const noexcept {
			const std::size_t used
				= static_cast<std::size_t>(tadpole_count_) + support_edge_count();
			if (used > static_cast<std::size_t>(N_EDGES)) {
				return std::numeric_limits<Int>::max();
			}
			return static_cast<Int>(static_cast<std::size_t>(N_EDGES) - used);
		}

		Int root_multiplicity(Int vertex) const {
			return has_root_support(vertex) ? Int{1} : Int{0};
		}

		Int support_degree(Int vertex) const {
			if (vertex >= N_VERTICES) {
				throw std::out_of_range("aggregate support vertex");
			}
			Int result = 0;
			for (Int other = 0; other < N_VERTICES; ++other) {
				if (other != vertex && has_support_edge(vertex, other)) {
					++result;
				}
			}
			return result;
		}

		template <typename Emit>
		void for_each_support_edge(Emit&& emit) const {
			for (Int first = 0; first < N_VERTICES; ++first) {
				for (Int second = static_cast<Int>(first + 1);
					 second < N_VERTICES; ++second) {
					if (has_support_edge(first, second)) {
						std::invoke(emit, first, second);
					}
				}
			}
		}

		template <typename Emit>
		void for_each_core_support_edge(Emit&& emit) const {
			for (Int first = 1; first < N_VERTICES; ++first) {
				for (Int second = static_cast<Int>(first + 1);
					 second < N_VERTICES; ++second) {
					if (has_support_edge(first, second)) {
						std::invoke(emit, first, second);
					}
				}
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
			return static_cast<std::size_t>(tadpole_count_) + support_edge_count()
				<= static_cast<std::size_t>(N_EDGES);
		}

		void validate() const {
			if (!structurally_valid()) {
				throw std::invalid_argument("invalid aggregate transient state");
			}
		}

		bool is_simple_family() const noexcept {
			return structurally_valid() && tadpole_count_ == 0
				&& excess_edge_count() == 0;
		}

		bool support_connected() const noexcept {
			if constexpr (N_VERTICES == 1) {
				return true;
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
					if (has_support_edge(current, candidate)) {
						visited[static_cast<std::size_t>(candidate)] = true;
						stack[stack_size++] = candidate;
						++visited_count;
					}
				}
			}
			return visited_count == static_cast<std::size_t>(N_VERTICES);
		}

		// Is there an allocation of the total excess over root groups that is
		// active and also realizes the requested split statuses?  selected means
		// one edge is removed from that group; exhausted means its multiplicity
		// was exactly one, so its root-support edge disappears in the child.
		bool active_allocation_feasible(
			std::uint64_t selected_groups = 0,
			std::uint64_t exhausted_groups = 0
		) const noexcept {
			if (!structurally_valid()
				|| (exhausted_groups & ~selected_groups) != 0
				|| (selected_groups & ~root_neighbor_mask()) != 0
				|| !support_connected()) {
				return false;
			}

			const Int root_degree = aggregate_root_valence();
			if (root_degree < 3) {
				return false;
			}

			std::size_t lower_sum = 0;
			std::size_t upper_sum = 0;
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				const Int degree = support_degree(vertex);
				const std::uint64_t bit = std::uint64_t{1} << vertex;
				if (!has_root_support(vertex)) {
					if (degree < 3 || degree >= root_degree) {
						return false;
					}
					continue;
				}

				Int lower = degree < 3 ? static_cast<Int>(3 - degree) : Int{0};
				if ((selected_groups & bit) != 0
					&& (exhausted_groups & bit) == 0) {
					lower = std::max(lower, Int{1});
				}
				const int upper = static_cast<int>(root_degree)
					- 1 - static_cast<int>(degree);
				if (upper < 0) {
					return false;
				}
				if ((exhausted_groups & bit) != 0) {
					if (lower != 0) {
						return false;
					}
					continue;
				}
				if (static_cast<int>(lower) > upper) {
					return false;
				}
				lower_sum += lower;
				upper_sum += static_cast<std::size_t>(upper);
			}

			const std::size_t excess = excess_edge_count();
			return lower_sum <= excess && excess <= upper_sum;
		}

		transient_classification classification() const noexcept {
			if (!structurally_valid() || !support_connected()) {
				return transient_classification::discard;
			}
			if (!is_simple_family()) {
				return active_allocation_feasible()
					? transient_classification::active_transient
					: transient_classification::discard;
			}

			std::array<Int, N_VERTICES> degrees{};
			Int maximum = 0;
			Int maximum_count = 0;
			for (Int vertex = 0; vertex < N_VERTICES; ++vertex) {
				const Int degree = support_degree(vertex);
				degrees[static_cast<std::size_t>(vertex)] = degree;
				if (degree < 3) {
					return transient_classification::discard;
				}
				if (degree > maximum) {
					maximum = degree;
					maximum_count = 1;
				} else if (degree == maximum) {
					++maximum_count;
				}
			}
			if (degrees[0] != maximum) {
				return transient_classification::discard;
			}
			return maximum_count == 1
				? transient_classification::active_admissible
				: transient_classification::terminal_admissible;
		}

		std::uint8_t classification_mask() const noexcept {
			return static_cast<std::uint8_t>(
				std::uint8_t{1}
				<< static_cast<std::uint8_t>(classification())
			);
		}

		bool allocation_feasible_for(transient_classification kind) const noexcept {
			return classification() == kind;
		}

		bool can_reduce_to_admissible() const noexcept {
			return classification() != transient_classification::discard;
		}

		analysis_type analyze() const noexcept {
			analysis_type result;
			result.structurally_valid = structurally_valid();
			if (!result.structurally_valid) {
				return result;
			}
			result.connected = support_connected();
			result.excess_edges = excess_edge_count();
			result.simple = is_simple_family();
			for (Int vertex = 0; vertex < N_VERTICES; ++vertex) {
				result.support_degrees[static_cast<std::size_t>(vertex)]
					= support_degree(vertex);
				result.maximum_support_degree = std::max(
					result.maximum_support_degree,
					result.support_degrees[static_cast<std::size_t>(vertex)]
				);
			}
			result.root_full_valence = aggregate_root_valence();
			result.root_is_reduced_maximum
				= result.support_degrees[0] == result.maximum_support_degree;
			result.has_active_allocation = active_allocation_feasible();
			return result;
		}

		hairless_graph_type to_hairless_graph() const {
			if (!is_simple_family()) {
				throw std::invalid_argument(
					"only an exact simple aggregate state has one hairless graph"
				);
			}
			hairless_graph_type result{};
			Int edge = 0;
			for_each_support_edge([&](Int first, Int second) {
				result.setEdge(edge++, first, second);
			});
			if (edge != N_EDGES) {
				throw std::logic_error("simple aggregate edge count mismatch");
			}
			return result;
		}

		hairless_graph_type without_hair() const { return to_hairless_graph(); }

		exact_type materialize(const excess_array& extras) const {
			validate();
			std::size_t total = 0;
			exact_type result;
			result.set_tadpole_count(tadpole_count_);
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				const Int extra = extras[static_cast<std::size_t>(vertex - 1)];
				if (!has_root_support(vertex) && extra != 0) {
					throw std::invalid_argument("excess on an absent root group");
				}
				total += extra;
				result.set_root_multiplicity(
					vertex,
					has_root_support(vertex) ? static_cast<Int>(extra + 1) : Int{0}
				);
			}
			if (total != excess_edge_count()) {
				throw std::invalid_argument("aggregate excess allocation mismatch");
			}
			for_each_core_support_edge([&](Int first, Int second) {
				result.set_core_edge(first, second);
			});
			result.validate();
			return result;
		}

		aggregate_transient_graph permuted_nonroot(
			const nonroot_permutation& old_to_new_local
		) const {
			validate_nonroot_permutation(old_to_new_local);
			aggregate_transient_graph result;
			result.tadpole_count_ = tadpole_count_;
			for_each_support_edge([&](Int old_first, Int old_second) {
				const Int first = old_first == 0 ? 0 : static_cast<Int>(
					old_to_new_local[static_cast<std::size_t>(old_first - 1)] + 1
				);
				const Int second = old_second == 0 ? 0 : static_cast<Int>(
					old_to_new_local[static_cast<std::size_t>(old_second - 1)] + 1
				);
				result.set_support_edge_unchecked(first, second, true);
			});
			return result;
		}

		// Enumerate the union of grouped split cases over all active exact
		// allocations represented by this family.  The child is again the bare
		// aggregate family: constraints such as the newly-created root bundle
		// having excess exactly separated_tadpoles are intentionally forgotten.
		template <typename Emit>
		std::uint64_t for_each_relevant_root_split(Emit&& emit) const {
			if (!active_allocation_feasible()) {
				return 0;
			}
			const Int root_valence = aggregate_root_valence();
			if (root_valence < 4) {
				return 0;
			}

			std::vector<Int> neighbors;
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				if (has_root_support(vertex)) {
					neighbors.push_back(vertex);
				}
			}

			std::uint64_t emitted = 0;
			for (Int separated = 0; separated <= tadpole_count_; ++separated) {
				auto visit_group_state
					= [&](std::uint64_t selected, std::uint64_t exhausted) {
						const std::size_t moved = static_cast<std::size_t>(separated)
							+ std::popcount(selected);
						const std::size_t retained
							= static_cast<std::size_t>(root_valence) - moved;
						if (moved < 2 || retained < 2 || moved > retained
							|| !active_allocation_feasible(selected, exhausted)) {
							return;
						}

						split_type child;
						child.tadpole_count_ = static_cast<Int>(
							tadpole_count_ - separated
						);
						for_each_support_edge([&](Int first, Int second) {
							child.set_support_edge_unchecked(first, second, true);
						});
						for (const Int vertex : neighbors) {
							const std::uint64_t bit = std::uint64_t{1} << vertex;
							if ((exhausted & bit) != 0) {
								child.set_support_edge_unchecked(0, vertex, false);
							}
							if ((selected & bit) != 0) {
								child.set_support_edge_unchecked(
									vertex, N_VERTICES, true
								);
							}
						}
						child.set_support_edge_unchecked(0, N_VERTICES, true);

						const transient_classification kind = child.classification();
						if (kind == transient_classification::discard) {
							return;
						}
						if (moved == retained && !child.is_simple_family()) {
							return;
						}

						symbolic_split_descriptor descriptor{
							separated, selected, exhausted
						};
						invoke_split_callback(
							emit, std::move(child), kind, descriptor
						);
						++emitted;
					};
				enumerate_group_states(
					neighbors, 0, 0, 0, visit_group_state
				);
			}
			return emitted;
		}

		signedInt compare(const aggregate_transient_graph& other) const noexcept {
			if (tadpole_count_ != other.tadpole_count_) {
				return tadpole_count_ < other.tadpole_count_ ? -1 : 1;
			}
			for (std::size_t i = 0; i < SUPPORT_WORD_COUNT; ++i) {
				if (support_words_[i] != other.support_words_[i]) {
					return support_words_[i] < other.support_words_[i] ? -1 : 1;
				}
			}
			return 0;
		}

		bool operator==(const aggregate_transient_graph&) const noexcept = default;

		std::size_t hash_value() const noexcept {
			std::size_t result = std::hash<Int>{}(tadpole_count_);
			for (const std::uint64_t word : support_words_) {
				result ^= std::hash<std::uint64_t>{}(word)
					+ 0x9e3779b97f4a7c15ULL + (result << 6) + (result >> 2);
			}
			return result;
		}

	private:
		template <Int, Int, typename>
		friend class aggregate_transient_graph;

		static std::size_t support_bit_index(Int first, Int second) {
			if (first >= N_VERTICES || second >= N_VERTICES || first == second) {
				throw std::out_of_range("aggregate support edge");
			}
			if (second < first) {
				std::swap(first, second);
			}
			const std::size_t u = first;
			const std::size_t v = second;
			const std::size_t n = N_VERTICES;
			return u * (2 * n - u - 1) / 2 + (v - u - 1);
		}

		void set_support_edge_unchecked(Int first, Int second, bool present) {
			const std::size_t bit = support_bit_index(first, second);
			const std::uint64_t mask = std::uint64_t{1} << (bit % 64);
			if (present) {
				support_words_[bit / 64] |= mask;
			} else {
				support_words_[bit / 64] &= ~mask;
			}
		}

		std::uint64_t root_neighbor_mask() const noexcept {
			std::uint64_t result = 0;
			for (Int vertex = 1; vertex < N_VERTICES; ++vertex) {
				if (has_root_support(vertex)) {
					result |= std::uint64_t{1} << vertex;
				}
			}
			return result;
		}

		Int aggregate_root_valence() const noexcept {
			return static_cast<Int>(
				2 * static_cast<unsigned>(tadpole_count_)
				+ static_cast<unsigned>(support_degree(0))
				+ static_cast<unsigned>(excess_edge_count())
			);
		}

		void validate_nonroot_permutation(
			const nonroot_permutation& permutation_value
		) const {
			std::array<bool, N_NONROOT_VERTICES> seen{};
			for (const Int image : permutation_value) {
				if (image >= N_NONROOT_VERTICES
					|| seen[static_cast<std::size_t>(image)]) {
					throw std::invalid_argument("invalid aggregate nonroot permutation");
				}
				seen[static_cast<std::size_t>(image)] = true;
			}
		}

		template <typename Visit>
		static void enumerate_group_states(
			const std::vector<Int>& neighbors,
			std::size_t index,
			std::uint64_t selected,
			std::uint64_t exhausted,
			Visit& visit
		) {
			if (index == neighbors.size()) {
				std::invoke(visit, selected, exhausted);
				return;
			}
			const std::uint64_t bit = std::uint64_t{1} << neighbors[index];
			enumerate_group_states(
				neighbors, index + 1, selected, exhausted, visit
			);
			enumerate_group_states(
				neighbors, index + 1, selected | bit, exhausted, visit
			);
			enumerate_group_states(
				neighbors, index + 1, selected | bit, exhausted | bit, visit
			);
		}

		template <typename Emit>
		static void invoke_split_callback(
			Emit& emit,
			split_type child,
			transient_classification kind,
			const symbolic_split_descriptor& descriptor
		) {
			if constexpr (std::is_invocable_v<
				Emit&, split_type, transient_classification,
				const symbolic_split_descriptor&
			>) {
				std::invoke(emit, std::move(child), kind, descriptor);
			} else {
				std::invoke(emit, std::move(child), kind);
			}
		}

		Int tadpole_count_ = 0;
		support_word_array support_words_{};
};

} // namespace GraphGeneration

namespace std {

template <Int N_VERTICES, Int N_EDGES, typename FieldType>
struct hash<GraphGeneration::aggregate_transient_graph<
	N_VERTICES, N_EDGES, FieldType
>> {
	std::size_t operator()(
		const GraphGeneration::aggregate_transient_graph<
			N_VERTICES, N_EDGES, FieldType
		>& graph
	) const noexcept {
		return graph.hash_value();
	}
};

} // namespace std
