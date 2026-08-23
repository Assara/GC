#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <functional>
#include <iostream>
#include <numeric>
#include <string_view>

#include "GraphGeneration/SupportTransientGraph.hpp"
#include "GraphGeneration/SupportTransientStandardizer.hpp"
#include "GraphGeneration/UnrootedSupportTransientGraph.hpp"
#include "GraphGeneration/UnrootedSupportTransientStandardizer.hpp"

namespace {

int failures = 0;

void check(bool condition, std::string_view description) {
	if (!condition) {
		std::cerr << "FAIL: " << description << '\n';
		++failures;
	}
}

template <typename Support>
void check_standardizer(const Support& graph, std::uint64_t expected_aut) {
	using Standardizer = GraphGeneration::support_transient_standardizer<
		Support::N_VERTICES_, Support::N_EDGES_, typename Support::field_type
	>;
	Standardizer standardizer;
	const auto result = standardizer.standardize_with_info(graph);

	typename Support::nonroot_permutation permutation{};
	std::iota(permutation.begin(), permutation.end(), Int{0});
	std::uint64_t automorphisms = 0;
	do {
		const Support relabeled = graph.permuted_nonroot(permutation);
		automorphisms += relabeled == graph;
		check(standardizer.standardize_no_sign(relabeled)
			== result.canonical_graph,
			"every nonroot relabeling has the same support representative");
	} while (std::next_permutation(permutation.begin(), permutation.end()));

	check(result.automorphism_count == expected_aut,
		"support standardizer returns the expected rooted automorphism order");
	check(result.automorphism_count == automorphisms,
		"support standardizer agrees with exhaustive automorphism counting");
	check(standardizer.standardize_no_sign(result.canonical_graph)
		== result.canonical_graph,
		"support standardization is idempotent");
}

void test_triangle_seed() {
	// At loop order five K3 has E=L+2=7 and four hidden edges.
	using Seed = GraphGeneration::support_transient_graph<3, 7>;
	const Seed seed = Seed::triangle();
	check(seed.support_edge_count() == 3,
		"the K3 seed has the three triangle support edges");
	check(seed.has_support_edge(0, 1)
			&& seed.has_support_edge(0, 2)
			&& seed.has_support_edge(1, 2),
		"the K3 seed support is a triangle");
	check(seed.total_surplus() == 4,
		"the loop-five K3 seed has four hidden edges");
	check(seed.has_minimum_support_degree_two(),
		"every K3 support vertex is bivalent");
	check(seed.has_repairable_support_valences(),
		"the K3 bivalent vertices are repairable from the surplus");
	check(seed.has_active_allocation(),
		"the K3 family has an active exact allocation");
	check_standardizer(seed, 2);
}

void test_repairable_bivalent_support_vertices() {
	using EdgeSeed = GraphGeneration::support_transient_graph<2, 6>;
	EdgeSeed edge;
	edge.set_support_edge(0, 1);
	edge.validate();
	check(!edge.has_minimum_support_degree_two(),
		"the two-vertex seed has support leaves");
	check(edge.has_repairable_support_valences(),
		"two hidden bundle copies repair the nonroot support leaf");
	check(edge.has_active_allocation(),
		"the one-edge support seed has an active trivalent allocation");
	bool saw_leaf_split = false;
	edge.for_each_symbolic_relevant_root_split_impl(
		[&](EdgeSeed::symbolic_split_type split) {
			saw_leaf_split |= split.selected_groups == 0;
		},
		false
	);
	check(saw_leaf_split,
		"a tadpole-funded split may create a support leaf with no moved group");

	using Support = GraphGeneration::support_transient_graph<4, 8>;
	Support k4_minus_edge = Support::wheel3();
	k4_minus_edge.set_support_edge(2, 3, false);
	k4_minus_edge.validate();
	check(k4_minus_edge.has_minimum_support_degree_two(),
		"K4 minus an edge still has support minimum degree two");
	check(!k4_minus_edge.has_minimum_support_degree_three(),
		"K4 minus an edge has bivalent support vertices");
	check(k4_minus_edge.has_repairable_support_valences(),
		"two root-adjacent bivalent vertices are repaired by two surplus edges");
	check(k4_minus_edge.has_active_allocation(),
		"the repairable support has an active exact allocation");
	check(k4_minus_edge.for_each_symbolic_relevant_root_split(
		[](Support::symbolic_split_type) {}) != 0,
		"a repairable support can emit children");

	using Insufficient = GraphGeneration::support_transient_graph<4, 6>;
	Insufficient too_little_surplus = Insufficient::wheel3();
	too_little_surplus.set_support_edge(2, 3, false);
	too_little_surplus.validate();
	check(!too_little_surplus.has_repairable_support_valences(),
		"one surplus edge cannot repair two bivalent root neighbours");
	check(!too_little_surplus.has_active_allocation(),
		"insufficient surplus gives no active exact allocation");
	check(too_little_surplus.for_each_symbolic_relevant_root_split(
		[](Insufficient::symbolic_split_type) {}) == 0,
		"insufficient surplus emits no children");

	using Away = GraphGeneration::support_transient_graph<5, 9>;
	Away away;
	for (Int rim = 1; rim < 5; ++rim) {
		away.set_support_edge(0, rim);
	}
	away.set_support_edge(1, 2);
	away.set_support_edge(2, 3);
	away.set_support_edge(3, 4);
	away.set_support_edge(4, 1);
	away.set_support_edge(0, 3, false);
	away.validate();
	check(away.support_degree(3) == 2 && !away.has_root_support(3),
		"the off-root fixture has a bivalent vertex away from the root");
	check(!away.has_repairable_support_valences(),
		"surplus cannot repair a bivalent vertex away from the root");
}

void test_symbolic_splits_and_active_roots() {
	// Loop eight: parent K3 has E=10 and seven hidden edges. This is large
	// enough to exercise both the old-root and changed-root witnesses without
	// running a graph enumeration.
	using Seed = GraphGeneration::support_transient_graph<3, 10>;
	using Child = typename Seed::split_type;
	using ChildStandardizer = GraphGeneration::support_transient_standardizer<
		Child::N_VERTICES_, Child::N_EDGES_, typename Child::field_type
	>;
	const Seed seed = Seed::triangle();
	ChildStandardizer standardizer;
	std::uint64_t children = 0;
	bool saw_active_root = false;
	seed.for_each_symbolic_relevant_root_split(
		[&](Seed::symbolic_split_type split) {
			++children;
			check(split.child.has_repairable_support_valences(),
				"every emitted child has repairable support valences");
			check(std::popcount(split.selected_groups) >= 1
					|| split.separated_tadpoles >= 2,
				"a support-leaf split is stabilized by at least two hidden incidences");
			check(std::popcount(seed.root_support_mask())
					- std::popcount(split.exhausted_groups) >= 2,
				"the retained split vertex keeps at least two old support groups");
			const std::size_t selected_support
				= std::popcount(split.selected_groups);
			const std::size_t retained_support
				= std::popcount(seed.root_support_mask())
					- std::popcount(split.exhausted_groups);
			check(selected_support <= retained_support,
				"the diverging side never has more support groups than the root");
			if (selected_support == retained_support
				&& split.active_root_mask != 0) {
				check((split.selected_groups & ~split.exhausted_groups) != 0,
					"a continuing support tie uses a virtual double edge");
			}
			check((split.active_root_mask & ~std::uint64_t{1}) == 0,
				"a continuing child always retains vertex zero as root");
			check(split.has_simple_image == split.child.has_simple_expansion(),
				"a feasible simple image is exactly a zero-surplus child");

			if ((split.active_root_mask & std::uint64_t{1}) != 0) {
				saw_active_root = true;
				check(split.child.has_repairable_support_valences(),
					"an active child remains support-repairable");
				check(split.child.has_active_allocation(),
					"an exact root witness projects to an active family");
				const Child canonical
					= standardizer.standardize_no_sign(split.child);
				check(standardizer.standardize_no_sign(canonical) == canonical,
					"a rooted child canonicalization is idempotent");
			}

		}
	);
	check(children != 0, "the loop-eight K3 seed emits symbolic children");
	check(saw_active_root, "some K3 child has an exact active-root witness");
}

void test_unrooted_standardizer_and_partition() {
	using Graph = GraphGeneration::unrooted_support_transient_graph<4, 8>;
	using Support = typename Graph::support_type;
	using Standardizer
		= GraphGeneration::unrooted_support_transient_standardizer<4, 8>;
	Support support = Support::wheel3();
	support.set_support_edge(2, 3, false);
	const Graph graph(support);
	Standardizer standardizer;
	const auto result = standardizer.standardize_with_info(graph);
	check(result.maximum_valence == 3 && result.maximum_valence_count == 2,
		"unrooted standardization recovers maximum support valence and tie count");
	check(result.canonical_graph.storage_partition() == std::pair<Int, Int>{3, 2},
		"the canonical graph reports its (maximum,tie) storage partition");
	const auto canonical_degrees
		= result.canonical_graph.support_valences();
	check(std::is_sorted(canonical_degrees.begin(), canonical_degrees.end()),
		"literal initial valence grouping places maximum vertices at the end");

	std::array<Int, 4> permutation{0, 1, 2, 3};
	do {
		const Graph relabeled = graph.permuted(permutation);
		check(standardizer.standardize_no_sign(relabeled)
				== result.canonical_graph,
			"every full vertex relabeling has one unrooted support representative");
	} while (std::next_permutation(permutation.begin(), permutation.end()));
}

void test_ordered_special_and_maximum_increasing_splits() {
	using Support = GraphGeneration::support_transient_graph<5, 13>;
	using Graph = GraphGeneration::unrooted_support_transient_graph<5, 13>;
	Support graph;
	for (Int first = 0; first < 5; ++first) {
		for (Int second = static_cast<Int>(first + 1); second < 5; ++second) {
			if (first != 0 || second != 4) {
				graph.set_support_edge(first, second);
			}
		}
	}
	graph.validate();
	check(graph.support_degree(0) == 3,
		"the special-phase fixture has a nonuniversal root");
	check(graph.has_active_allocation(),
		"the special-phase fixture has an active allocation");
	bool saw_nonuniversal_root_leaf_split = false;
	Graph(graph).for_each_allowed_root_split(
		[&](Graph::symbolic_split_type split, Int root) {
			saw_nonuniversal_root_leaf_split
				|= root == 0 && split.selected_groups == 0;
		}
	);
	check(saw_nonuniversal_root_leaf_split,
		"a global universal vertex opens the special-split phase for every feasible root");

	using NonuniversalSupport
		= GraphGeneration::support_transient_graph<6, 12>;
	using NonuniversalGraph
		= GraphGeneration::unrooted_support_transient_graph<6, 12>;
	NonuniversalSupport k33;
	for (Int left = 0; left < 3; ++left) {
		for (Int right = 3; right < 6; ++right) {
			k33.set_support_edge(left, right);
		}
	}
	k33.validate();
	check(NonuniversalGraph(k33).maximum_support_valence() == 3,
		"the nonuniversal fixture has maximum support valence three");
	std::uint64_t ordinary_children = 0;
	NonuniversalGraph(k33).for_each_allowed_root_split(
		[&](NonuniversalGraph::symbolic_split_type split, Int) {
			++ordinary_children;
			check(split.selected_groups != 0,
				"a nonuniversal parent emits no support-leaf split");
			check(split.child.maximum_support_valence() <= 3,
				"a nonuniversal parent emits no maximum-increasing split");
		}
	);
	check(ordinary_children != 0,
		"the nonuniversal fixture still emits ordered ordinary splits");
}

} // namespace

int main() {
	const auto run = [](std::string_view name, auto&& test) {
		try {
			std::invoke(test);
		} catch (const std::exception& exception) {
			std::cerr << "FAIL: " << name << " threw: "
			          << exception.what() << '\n';
			++failures;
		}
	};
	run("K3 seed", test_triangle_seed);
	run("repairable bivalent support vertices",
		test_repairable_bivalent_support_vertices);
	run("symbolic splits and active roots",
		test_symbolic_splits_and_active_roots);
	run("ordered special and maximum-increasing splits",
		test_ordered_special_and_maximum_increasing_splits);
	run("unrooted standardizer and partition",
		test_unrooted_standardizer_and_partition);

	if (failures != 0) {
		std::cerr << failures << " support transient test(s) failed\n";
		return 1;
	}
	std::cout << "all support transient tests passed\n";
	return 0;
}
