#include <algorithm>
#include <array>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include "GraphGeneration/TransientGraph.hpp"

namespace {

int failures = 0;

void check(bool condition, std::string_view description) {
	if (!condition) {
		std::cerr << "FAIL: " << description << '\n';
		++failures;
	}
}

template <typename GraphType>
GraphType graph_with_edges(
	std::initializer_list<std::pair<Int, Int>> edges
) {
	if (edges.size() != GraphType::N_EDGES_) {
		throw std::invalid_argument("wrong edge count in transient graph test");
	}
	GraphType graph;
	Int edge_index = 0;
	for (const auto [first, second] : edges) {
		graph.setEdge(edge_index++, first, second);
	}
	return graph;
}

template <typename GraphType>
using UnsignedStandardizer = GraphStandardizer<
	GraphType::N_VERTICES_,
	GraphType::N_EDGES_,
	GraphType::N_OUT_HAIR_,
	GraphType::N_IN_HAIR_,
	GraphType::C_,
	GraphType::D_,
	typename GraphType::Field
>;

template <typename Transient>
struct relevant_split_sets {
	using Split = typename Transient::split_type;
	using Hairy = typename Split::graph_type;
	using Hairless = typename Split::hairless_graph_type;

	std::unordered_set<Hairy> active;
	std::unordered_set<Hairless> admissible;
	std::uint64_t callbacks = 0;
};

template <typename Transient>
void insert_relevant_split(
	relevant_split_sets<Transient>& result,
	typename Transient::split_type child,
	GraphGeneration::transient_classification kind
) {
	using Split = typename Transient::split_type;
	using Hairy = typename Split::graph_type;
	using Hairless = typename Split::hairless_graph_type;
	static UnsignedStandardizer<Hairy> hairy_standardizer;
	static UnsignedStandardizer<Hairless> hairless_standardizer;

	++result.callbacks;
	if (GraphGeneration::is_active(kind)) {
		result.active.insert(
			hairy_standardizer.standardize_no_sign(child.graph())
		);
	}
	if (GraphGeneration::is_admissible(kind)) {
		result.admissible.insert(
			hairless_standardizer.standardize_no_sign(child.without_hair())
		);
	}
}

template <typename Transient>
relevant_split_sets<Transient> exhaustive_relevant_split_sets(
	const Transient& parent
) {
	relevant_split_sets<Transient> result;
	parent.for_each_root_split([&](typename Transient::split_type child) {
		const auto kind = child.classification();
		if (kind != GraphGeneration::transient_classification::discard) {
			insert_relevant_split(result, std::move(child), kind);
		}
	});
	return result;
}

template <typename Transient>
relevant_split_sets<Transient> optimized_relevant_split_sets(
	const Transient& parent
) {
	relevant_split_sets<Transient> result;
	parent.for_each_relevant_root_split(
		[&](
			typename Transient::split_type child,
			GraphGeneration::transient_classification kind
		) {
			insert_relevant_split(result, std::move(child), kind);
		}
	);
	return result;
}

void test_bounded_count_combinatorics() {
	const std::array<Int, 2> bounds{1, 2};
	std::vector<std::array<Int, 2>> values;
	const std::uint64_t count = combutils::for_each_bounded_count_vector(
		bounds,
		[&](std::span<const Int> value) {
			values.push_back({value[0], value[1]});
		}
	);
	check(count == 6 && values.size() == 6,
		"bounded count vectors enumerate their Cartesian product");
	check(values.front() == std::array<Int, 2>{0, 0}
		&& values.back() == std::array<Int, 2>{1, 2},
		"bounded count vectors have deterministic odometer order");

	const std::array<Int, 2> group_sizes{1, 2};
	const std::array<Int, 2> first{0, 1};
	const std::array<Int, 2> second{1, 1};
	check(combutils::is_first_grouped_bipartition(first, group_sizes),
		"the first grouped bipartition orientation is retained");
	check(!combutils::is_first_grouped_bipartition(second, group_sizes),
		"the complementary grouped bipartition orientation is rejected");
}

void test_hair_is_only_a_marker() {
	using Transient = GraphGeneration::transient_graph<4, 6>;
	using Hairless = Transient::hairless_graph_type;
	const Hairless k4 = graph_with_edges<Hairless>({
		{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
	});
	const Transient transient(k4, 0);

	check(transient.root() == 0, "the hair endpoint is the root");
	check(
		transient.edge_incidence_positions(0) == std::vector<Int>({1, 3, 5}),
		"ordinary incidence positions exclude the hair at position zero"
	);
	check(
		transient.ordinary_valences() == std::array<Int, 4>({3, 3, 3, 3}),
		"ordinary valence excludes the marker hair"
	);
	check(
		transient.graph().valence_array() == std::array<Int, 4>({4, 3, 3, 3}),
		"the backing Graph still sees its actual hair"
	);
	check(
		transient.classification()
			== GraphGeneration::transient_classification::terminal_admissible,
		"a simple graph with tied maximum valence is final but terminal"
	);
	check(transient.without_hair() == k4, "removing the hair preserves raw labels");
	check(
		transient.remove_hair_if_admissible().value() == k4,
		"an admissible underlying graph can be extracted"
	);
}

void test_direct_split_moves_the_marker_to_the_larger_side() {
	using Transient = GraphGeneration::transient_graph<7, 6>;
	using Hairless = Transient::hairless_graph_type;
	const Hairless star = graph_with_edges<Hairless>({
		{0, 1}, {0, 2}, {0, 3}, {0, 4}, {0, 5}, {0, 6}
	});
	const Transient transient(star, 0);

	const std::array<Int, 2> move_two{1, 2};
	const auto split_two = transient.split_root(move_two);
	check(
		split_two.ordinary_valences()[0] == 5
			&& split_two.ordinary_valences()[7] == 3,
		"a 6-valent root split with two moved incidences is 5--3"
	);
	check(split_two.root() == 0, "the marker stays on the 5-valent old side");
	check(
		split_two.graph().getEdge(6) == std::pair<Int, Int>{0, 7},
		"the new splitting edge joins the two split vertices"
	);

	const std::array<Int, 4> move_four{1, 2, 3, 4};
	const auto split_four = transient.split_root(move_four);
	check(
		split_four.ordinary_valences()[0] == 3
			&& split_four.ordinary_valences()[7] == 5,
		"moving four incidences reverses the 5--3 split"
	);
	check(split_four.root() == 7, "the marker moves to the 5-valent new side");

	const std::array<Int, 3> move_three{1, 2, 3};
	const auto split_three = transient.split_root(move_three);
	check(
		split_three.ordinary_valences()[0] == 4
			&& split_three.ordinary_valences()[7] == 4,
		"three moved incidences give the equal 4--4 split"
	);
	check(split_three.root() == 0, "an equal split has a deterministic marker side");
}

void test_unoriented_split_count() {
	using Transient = GraphGeneration::transient_graph<1, 3>;
	const Transient rose;
	std::uint64_t old_roots = 0;
	std::uint64_t new_roots = 0;
	const std::uint64_t count = rose.for_each_root_split([&](const auto& child) {
		if (child.root() == 0) {
			++old_roots;
		} else if (child.root() == 1) {
			++new_roots;
		} else {
			check(false, "a split root is one of the two split vertices");
		}
		const auto incidences = child.edge_incidence_positions(child.root());
		check(
			incidences.end() == std::find(incidences.begin(), incidences.end(), Int{0}),
			"the artificial hair is never an ordinary split incidence"
		);
	});
	check(count == 25, "a 6-valent root has 25 stable unoriented splits");
	check(old_roots == 20 && new_roots == 5, "the larger-side root rule is deterministic");

	std::uint64_t relevant_callbacks = 0;
	rose.for_each_relevant_root_split(
		[&](const auto& child, GraphGeneration::transient_classification kind) {
			++relevant_callbacks;
			check(kind != GraphGeneration::transient_classification::discard,
				"the filtered split enumerator never emits discard states");
			check(child.classification() == kind,
				"the emitted classification describes the child");
		}
	);
	check(relevant_callbacks == 1,
		"grouped generation emits one relevant split of the three-tadpole rose");

	const auto exhaustive = exhaustive_relevant_split_sets(rose);
	const auto optimized = optimized_relevant_split_sets(rose);
	check(exhaustive.active == optimized.active
		&& exhaustive.admissible == optimized.admissible,
		"grouped rose splits preserve every relevant isomorphism class");
	check(exhaustive.callbacks == 12 && optimized.callbacks == 1,
		"grouped rose splits discard choices of tadpole and endpoint");
}

void test_ten_tadpole_grouped_split_count() {
	using Transient = GraphGeneration::transient_graph<1, 10>;
	const Transient rose;
	std::array<bool, 11> root_loop_counts{};
	std::uint64_t callbacks = 0;
	const std::uint64_t returned = rose.for_each_relevant_root_split(
		[&](const auto& child, GraphGeneration::transient_classification kind) {
			++callbacks;
			check(kind == GraphGeneration::transient_classification::active_transient,
				"every relevant ten-tadpole rose child remains transient");
			check(child.root() == 0,
				"the larger side of a ten-tadpole rose split remains the root");

			Int root_loops = 0;
			Int cross_edges = 0;
			bool reached_retained_tadpoles = false;
			for (Int edge = 0; edge < child.N_EDGES_; ++edge) {
				const auto [first, second] = child.graph().getEdge(edge);
				root_loops += first == 0 && second == 0;
				cross_edges += first != second;
				if (edge < 10) {
					if (first == 0 && second == 0) {
						reached_retained_tadpoles = true;
					} else {
						check(!reached_retained_tadpoles && first == 0 && second == 1,
							"the lowest tadpoles are separated in one fixed direction");
					}
				}
			}
			root_loop_counts[static_cast<std::size_t>(root_loops)] = true;
			check(cross_edges == static_cast<Int>(11 - root_loops),
				"each separated tadpole adds one edge to the split-edge bundle");
		}
	);
	check(returned == 8 && callbacks == 8,
		"ten tadpoles require only the eight stable unequal grouped splits");
	for (Int loops = 1; loops <= 8; ++loops) {
		check(root_loop_counts[static_cast<std::size_t>(loops)],
			"ten-tadpole grouped splits cover every surviving separation count");
	}
}

void test_parallel_bundle_grouped_splits() {
	using Transient = GraphGeneration::transient_graph<5, 9>;
	using Hairless = Transient::hairless_graph_type;

	const auto check_fixture = [&]
		(const Hairless& graph, Int lower_edge, Int upper_edge,
		 std::string_view description) {
		const Transient parent(graph, 0);
		check(parent.classification()
			== GraphGeneration::transient_classification::active_transient,
			"parallel-bundle fixture is an active transient");
		const auto exhaustive = exhaustive_relevant_split_sets(parent);
		const auto optimized = optimized_relevant_split_sets(parent);
		check(exhaustive.active == optimized.active
			&& exhaustive.admissible == optimized.admissible,
			description);
		check(exhaustive.callbacks == 6 && optimized.callbacks == 3,
			"a double-edge group moves only its lowest edge");
		parent.for_each_relevant_root_split(
			[&](const auto& child, GraphGeneration::transient_classification) {
				const auto lower = child.graph().getEdge(lower_edge);
				const auto upper = child.graph().getEdge(upper_edge);
				const bool lower_moved = lower.first == 5 || lower.second == 5;
				const bool upper_moved = upper.first == 5 || upper.second == 5;
				check(lower_moved && !upper_moved,
					"only the lowest-indexed edge of a parallel group may move");
			}
		);
	};

	check_fixture(graph_with_edges<Hairless>({
		{0, 1}, {0, 2}, {0, 2}, {0, 3}, {0, 4},
		{1, 2}, {1, 4}, {2, 3}, {3, 4}
	}), 1, 2,
		"grouped splits preserve classes when the double edge follows a singleton");

	check_fixture(graph_with_edges<Hairless>({
		{0, 1}, {0, 1}, {0, 2}, {0, 3}, {0, 4},
		{1, 2}, {1, 4}, {2, 3}, {3, 4}
	}), 0, 1,
		"grouped splits preserve classes when the first root group is double");
}

void test_admissibility_filters() {
	using Classification = GraphGeneration::transient_classification;

	{
		using Transient = GraphGeneration::transient_graph<3, 6>;
		using Hairless = Transient::hairless_graph_type;
		const Transient graph(graph_with_edges<Hairless>({
			{0, 0}, {0, 1}, {0, 1}, {0, 2}, {0, 2}, {1, 2}
		}), 0);
		const auto state = graph.analyze();
		check(state.valences == std::array<Int, 3>({6, 3, 3}),
			"a tadpole contributes two to ordinary valence");
		check(state.reduced_valences == std::array<Int, 3>({2, 2, 2}),
			"reduced valence removes loops and collapses parallel bundles");
		check(graph.classification() == Classification::active_transient,
			"root defects with a reduced-valence tie remain transient");
		check(!graph.remove_hair_if_admissible().has_value(),
			"a nonsimple transient is not a final admissible graph");
	}

	{
		using Transient = GraphGeneration::transient_graph<3, 7>;
		using Hairless = Transient::hairless_graph_type;
		const Transient off_root_loop(graph_with_edges<Hairless>({
			{0, 1}, {0, 1}, {0, 1}, {0, 2}, {0, 2}, {0, 2}, {1, 1}
		}), 0);
		check(!off_root_loop.analyze().defects_at_root,
			"a tadpole away from the root is a permanent defect");
		check(off_root_loop.classification() == Classification::discard,
			"a transient with an off-root tadpole is discarded");
	}

	{
		using Transient = GraphGeneration::transient_graph<3, 8>;
		using Hairless = Transient::hairless_graph_type;
		const Transient off_root_bundle(graph_with_edges<Hairless>({
			{0, 1}, {0, 1}, {0, 1}, {0, 2}, {0, 2}, {0, 2}, {1, 2}, {1, 2}
		}), 0);
		check(!off_root_bundle.analyze().defects_at_root,
			"a parallel bundle away from the root is a permanent defect");
		check(off_root_bundle.classification() == Classification::discard,
			"a transient with an off-root parallel bundle is discarded");
	}

	{
		using Transient = GraphGeneration::transient_graph<5, 14>;
		using Hairless = Transient::hairless_graph_type;
		const Transient reduced_failure(graph_with_edges<Hairless>({
			{0, 1}, {0, 1}, {0, 1}, {0, 1},
			{0, 2}, {0, 2}, {0, 2}, {0, 2},
			{1, 2}, {1, 3}, {1, 4}, {2, 3}, {2, 4}, {3, 4}
		}), 0);
		const auto state = reduced_failure.analyze();
		check(state.root_is_unique_maximum && state.defects_at_root,
			"the reduced-valence fixture passes the earlier transient tests");
		check(!state.root_is_reduced_maximum,
			"the root can lose maximality after defects are discounted");
		check(reduced_failure.classification() == Classification::discard,
			"a reduced-valence failure is discarded");
	}

	{
		using Transient = GraphGeneration::transient_graph<2, 3>;
		using Hairless = Transient::hairless_graph_type;
		const Transient theta(graph_with_edges<Hairless>({
			{0, 1}, {0, 1}, {0, 1}
		}), 0);
		check(theta.classification() == Classification::discard,
			"a nonsimple graph with tied maximum valence is terminal and discarded");
	}

	{
		using Transient = GraphGeneration::transient_graph<5, 8>;
		using Hairless = Transient::hairless_graph_type;
		const Hairless wheel = graph_with_edges<Hairless>({
			{0, 1}, {0, 2}, {0, 3}, {0, 4},
			{1, 2}, {2, 3}, {3, 4}, {4, 1}
		});
		check(Transient(wheel, 0).classification() == Classification::active_admissible,
			"a simple graph with a unique maximum is final and remains active");
		check(Transient(wheel, 1).classification() == Classification::discard,
			"a split result whose marker misses the maximum is pruned");
	}
}

void test_invalid_split_descriptions_are_rejected() {
	using Transient = GraphGeneration::transient_graph<1, 3>;
	const Transient rose;
	bool rejected_zero = false;
	try {
		const std::array<Int, 2> invalid{0, 1};
		(void)rose.split_root(invalid);
	} catch (const std::invalid_argument&) {
		rejected_zero = true;
	}
	check(rejected_zero, "the fixed incidence cannot be moved");

	bool rejected_duplicate = false;
	try {
		const std::array<Int, 2> invalid{1, 1};
		(void)rose.split_root(invalid);
	} catch (const std::invalid_argument&) {
		rejected_duplicate = true;
	}
	check(rejected_duplicate, "an incidence cannot be moved twice");
}

} // namespace

int main() {
	test_bounded_count_combinatorics();
	test_hair_is_only_a_marker();
	test_direct_split_moves_the_marker_to_the_larger_side();
	test_unoriented_split_count();
	test_ten_tadpole_grouped_split_count();
	test_parallel_bundle_grouped_splits();
	test_admissibility_filters();
	test_invalid_split_descriptions_are_rejected();

	if (failures != 0) {
		std::cerr << failures << " transient graph test(s) failed\n";
		return 1;
	}
	std::cout << "all transient graph tests passed\n";
	return 0;
}
