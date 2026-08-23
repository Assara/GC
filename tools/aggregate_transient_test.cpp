#include <algorithm>
#include <array>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>

#include "GraphGeneration/AggregateTransientGraph.hpp"
#include "GraphGeneration/AggregateTransientStandardizer.hpp"
#include "GraphGeneration/RootedTransientStandardizer.hpp"

namespace {

int failures = 0;

void check(bool condition, std::string_view description) {
	if (!condition) {
		std::cerr << "FAIL: " << description << '\n';
		++failures;
	}
}

template <typename Aggregate, typename Emit>
void for_each_exact_allocation(const Aggregate& aggregate, Emit&& emit) {
	typename Aggregate::excess_array extras{};
	std::array<Int, Aggregate::N_NONROOT_VERTICES> root_neighbors{};
	std::size_t root_neighbor_count = 0;
	for (Int vertex = 1; vertex < Aggregate::N_VERTICES_; ++vertex) {
		if (aggregate.has_root_support(vertex)) {
			root_neighbors[root_neighbor_count++] = vertex;
		}
	}

	const auto enumerate = [&](auto&& self, std::size_t neighbor,
		Int remaining) -> void {
		if (neighbor == root_neighbor_count) {
			if (remaining == 0) {
				std::invoke(emit, aggregate.materialize(extras));
			}
			return;
		}

		const Int vertex = root_neighbors[neighbor];
		const std::size_t index = static_cast<std::size_t>(vertex - 1);
		if (neighbor + 1 == root_neighbor_count) {
			extras[index] = remaining;
			self(self, neighbor + 1, Int{0});
			extras[index] = 0;
			return;
		}
		for (Int value = 0; value <= remaining; ++value) {
			extras[index] = value;
			self(self, neighbor + 1, static_cast<Int>(remaining - value));
		}
		extras[index] = 0;
	};

	if (root_neighbor_count == 0) {
		if (aggregate.excess_edge_count() == 0) {
			std::invoke(emit, aggregate.materialize(extras));
		}
		return;
	}
	enumerate(enumerate, 0, aggregate.excess_edge_count());
}

template <typename Aggregate>
struct brute_force_result {
	Aggregate canonical_graph;
	std::uint64_t automorphism_count = 0;
};

template <typename Aggregate>
brute_force_result<Aggregate> brute_force_standardization(
	const Aggregate& input
) {
	typename Aggregate::nonroot_permutation permutation{};
	std::iota(permutation.begin(), permutation.end(), Int{0});

	brute_force_result<Aggregate> result;
	bool have_candidate = false;
	do {
		Aggregate candidate = input.permuted_nonroot(permutation);
		if (!have_candidate || result.canonical_graph.compare(candidate) < 0) {
			result.canonical_graph = std::move(candidate);
			result.automorphism_count = 1;
			have_candidate = true;
		} else if (result.canonical_graph == candidate) {
			++result.automorphism_count;
		}
	} while (std::next_permutation(permutation.begin(), permutation.end()));
	return result;
}

template <typename Aggregate>
void check_standardizer(const Aggregate& graph, std::uint64_t expected_aut) {
	using Standardizer = GraphGeneration::aggregate_transient_standardizer<
		Aggregate::N_VERTICES_, Aggregate::N_EDGES_,
		typename Aggregate::field_type
	>;
	Standardizer standardizer;
	const auto brute = brute_force_standardization(graph);
	const auto actual = standardizer.standardize_with_info(graph);

	check(actual.automorphism_count == expected_aut,
		"aggregate standardizer returns the expected rooted automorphism order");
	check(actual.automorphism_count == brute.automorphism_count,
		"aggregate standardizer agrees with exhaustive automorphism counting");
	check(standardizer.standardize_no_sign(brute.canonical_graph)
		== actual.canonical_graph,
		"aggregate standardization is constant on the exhaustive orbit");
	check(standardizer.standardize_no_sign(actual.canonical_graph)
		== actual.canonical_graph,
		"aggregate standardization is idempotent");

	typename Aggregate::nonroot_permutation permutation{};
	std::iota(permutation.begin(), permutation.end(), Int{0});
	do {
		check(standardizer.standardize_no_sign(
			graph.permuted_nonroot(permutation)
		) == actual.canonical_graph,
			"every nonroot relabeling has the same aggregate representative");
	} while (std::next_permutation(permutation.begin(), permutation.end()));
}

void test_derived_excess_and_collision() {
	using Exact = GraphGeneration::rooted_transient_graph<3, 8>;
	using Aggregate = GraphGeneration::aggregate_transient_graph<3, 8>;

	Exact two_four;
	two_four.set_tadpole_count(1);
	two_four.set_root_multiplicity(1, 2);
	two_four.set_root_multiplicity(2, 4);
	two_four.set_core_edge(1, 2);
	two_four.validate();

	Exact three_three;
	three_three.set_tadpole_count(1);
	three_three.set_root_multiplicity(1, 3);
	three_three.set_root_multiplicity(2, 3);
	three_three.set_core_edge(1, 2);
	three_three.validate();

	const Aggregate collapsed_two_four(two_four);
	const Aggregate collapsed_three_three(three_three);
	check(collapsed_two_four.support_edge_count() == 3,
		"the weighted rooted triangle collapses to three support edges");
	check(collapsed_two_four.excess_edge_count() == 4,
		"aggregate excess is derived from E - tadpoles - support edges");
	check(collapsed_two_four.tadpole_count() == 1,
		"collapse preserves the root tadpole count");
	check(collapsed_two_four == collapsed_three_three,
		"the {2,4} and {3,3} root bundles intentionally share one family");

	Aggregate::excess_array extras{};
	extras[0] = 1;
	extras[1] = 3;
	check(collapsed_two_four.materialize(extras) == two_four,
		"a requested excess allocation materializes the {2,4} graph exactly");
	extras[0] = 2;
	extras[1] = 2;
	check(collapsed_two_four.materialize(extras) == three_three,
		"a requested excess allocation materializes the {3,3} graph exactly");

	std::uint64_t allocation_count = 0;
	for_each_exact_allocation(collapsed_two_four, [&](const Exact&) {
		++allocation_count;
	});
	check(allocation_count == 5,
		"four excess edges over two root groups have five labeled allocations");
	check_standardizer(collapsed_two_four, 2);
}

void test_canonical_relabeling_and_automorphisms() {
	using Aggregate = GraphGeneration::aggregate_transient_graph<4, 8>;
	Aggregate rooted_k4_family;
	rooted_k4_family.set_tadpole_count(1);
	for (Int first = 0; first < 4; ++first) {
		for (Int second = static_cast<Int>(first + 1); second < 4; ++second) {
			rooted_k4_family.set_support_edge(first, second);
		}
	}
	rooted_k4_family.validate();
	check(rooted_k4_family.excess_edge_count() == 1,
		"the rooted K4 family carries one derived excess edge");
	check_standardizer(rooted_k4_family, 6);
}

template <typename Aggregate>
void check_symbolic_child_union(
	const Aggregate& input,
	std::string_view description
) {
	using Exact = typename Aggregate::exact_type;
	using ExactChild = typename Exact::split_type;
	using AggregateChild = typename Aggregate::split_type;
	using ExactStandardizer = GraphGeneration::rooted_transient_standardizer<
		Exact::N_VERTICES_, Exact::N_EDGES_, typename Exact::field_type
	>;
	using ChildAggregateStandardizer
		= GraphGeneration::aggregate_transient_standardizer<
			AggregateChild::N_VERTICES_, AggregateChild::N_EDGES_,
			typename AggregateChild::field_type
		>;

	ExactStandardizer exact_standardizer;
	ChildAggregateStandardizer child_standardizer;
	std::unordered_set<Exact> exact_parents;
	std::unordered_set<AggregateChild> exact_image;
	std::unordered_set<AggregateChild> symbolic_image;

	for_each_exact_allocation(input, [&](Exact exact) {
		if (GraphGeneration::is_active(exact.classification())) {
			exact_parents.insert(exact_standardizer.standardize_no_sign(exact));
		}
	});
	for (const Exact& exact : exact_parents) {
		exact.for_each_relevant_root_split(
			[&](ExactChild child, GraphGeneration::transient_classification) {
				exact_image.insert(child_standardizer.standardize_no_sign(
					AggregateChild(child)
				));
			}
		);
	}

	const auto canonical_input
		= GraphGeneration::aggregate_transient_standardizer<
			Aggregate::N_VERTICES_, Aggregate::N_EDGES_,
			typename Aggregate::field_type
		>{}.standardize_no_sign(input);
	canonical_input.for_each_relevant_root_split(
		[&](AggregateChild child, GraphGeneration::transient_classification,
			const typename Aggregate::symbolic_split_descriptor&) {
			symbolic_image.insert(
				child_standardizer.standardize_no_sign(child)
			);
		}
	);

	check(symbolic_image == exact_image, description);
}

void test_symbolic_split_union() {
	using Collision = GraphGeneration::aggregate_transient_graph<3, 8>;
	Collision collision;
	collision.set_tadpole_count(1);
	collision.set_support_edge(0, 1);
	collision.set_support_edge(0, 2);
	collision.set_support_edge(1, 2);
	collision.validate();
	check_symbolic_child_union(collision,
		"symbolic splitting equals the exact allocation union for the collision family");

	using Asymmetric = GraphGeneration::aggregate_transient_graph<4, 9>;
	Asymmetric asymmetric;
	asymmetric.set_tadpole_count(1);
	asymmetric.set_support_edge(0, 1);
	asymmetric.set_support_edge(0, 2);
	asymmetric.set_support_edge(0, 3);
	asymmetric.set_support_edge(1, 2);
	asymmetric.set_support_edge(2, 3);
	asymmetric.validate();
	check_symbolic_child_union(asymmetric,
		"symbolic splitting equals the exact allocation union on an asymmetric core");
}

template <Int LOOP, Int VERTICES>
void compare_in_memory_generation(
	const std::unordered_set<GraphGeneration::rooted_transient_graph<
		VERTICES, static_cast<Int>(LOOP + VERTICES - 1)
	>>& exact_frontier,
	const std::unordered_set<GraphGeneration::aggregate_transient_graph<
		VERTICES, static_cast<Int>(LOOP + VERTICES - 1)
	>>& aggregate_frontier
) {
	static constexpr Int EDGES = static_cast<Int>(LOOP + VERTICES - 1);
	using Exact = GraphGeneration::rooted_transient_graph<VERTICES, EDGES>;
	using Aggregate = GraphGeneration::aggregate_transient_graph<VERTICES, EDGES>;
	using ExactChild = typename Exact::split_type;
	using AggregateChild = typename Aggregate::split_type;
	using HairlessChild = typename ExactChild::hairless_graph_type;
	using ExactChildStandardizer = GraphGeneration::rooted_transient_standardizer<
		ExactChild::N_VERTICES_, ExactChild::N_EDGES_
	>;
	using AggregateChildStandardizer
		= GraphGeneration::aggregate_transient_standardizer<
			AggregateChild::N_VERTICES_, AggregateChild::N_EDGES_
		>;
	using HairlessStandardizer = GraphStandardizer<
		HairlessChild::N_VERTICES_, HairlessChild::N_EDGES_, 0, 0, 0, 0,
		typename HairlessChild::Field
	>;

	ExactChildStandardizer exact_standardizer;
	AggregateChildStandardizer aggregate_standardizer;
	HairlessStandardizer hairless_standardizer;
	std::unordered_set<ExactChild> exact_next;
	std::unordered_set<AggregateChild> aggregate_next;
	std::unordered_set<AggregateChild> collapsed_exact_next;
	std::unordered_set<HairlessChild> exact_final;
	std::unordered_set<HairlessChild> aggregate_final;

	for (const Exact& parent : exact_frontier) {
		parent.for_each_relevant_root_split(
			[&](ExactChild child, GraphGeneration::transient_classification kind) {
				child = exact_standardizer.standardize_no_sign(child);
				if (GraphGeneration::is_active(kind)) {
					exact_next.insert(child);
				}
				if (GraphGeneration::is_admissible(kind)) {
					exact_final.insert(hairless_standardizer.standardize_no_sign(
						child.to_hairless_graph()
					));
				}
			}
		);
	}
	for (const ExactChild& child : exact_next) {
		collapsed_exact_next.insert(aggregate_standardizer.standardize_no_sign(
			AggregateChild(child)
		));
	}

	for (const Aggregate& parent : aggregate_frontier) {
		parent.for_each_relevant_root_split(
			[&](AggregateChild child,
				GraphGeneration::transient_classification kind) {
				child = aggregate_standardizer.standardize_no_sign(child);
				if (GraphGeneration::is_active(kind)) {
					aggregate_next.insert(child);
				}
				if (GraphGeneration::is_admissible(kind)) {
					aggregate_final.insert(
						hairless_standardizer.standardize_no_sign(
							child.to_hairless_graph()
						)
					);
				}
			}
		);
	}

	const std::string prefix = "loop " + std::to_string(LOOP)
		+ ", V=" + std::to_string(static_cast<unsigned>(VERTICES + 1));
	check(aggregate_next == collapsed_exact_next,
		prefix + ": aggregate active frontier equals the collapsed exact frontier");
	check(aggregate_final == exact_final,
		prefix + ": aggregate and exact final simple sets agree");

	if constexpr (VERTICES + 1 < 2 * LOOP - 2) {
		compare_in_memory_generation<LOOP, static_cast<Int>(VERTICES + 1)>(
			exact_next, aggregate_next
		);
	}
}

template <Int LOOP>
void test_loop_generation() {
	using Exact = GraphGeneration::rooted_transient_graph<1, LOOP>;
	using Aggregate = GraphGeneration::aggregate_transient_graph<1, LOOP>;
	std::unordered_set<Exact> exact{Exact::rose()};
	std::unordered_set<Aggregate> aggregate{Aggregate::rose()};
	compare_in_memory_generation<LOOP, 1>(exact, aggregate);
}

} // namespace

int main() {
	test_derived_excess_and_collision();
	test_canonical_relabeling_and_automorphisms();
	test_symbolic_split_union();
	test_loop_generation<3>();
	test_loop_generation<4>();
	test_loop_generation<5>();
	test_loop_generation<6>();

	if (failures != 0) {
		std::cerr << failures << " aggregate transient test(s) failed\n";
		return 1;
	}
	std::cout << "all aggregate transient tests passed\n";
	return 0;
}
