#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <initializer_list>
#include <iostream>
#include <numeric>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <utility>

#include "GraphGeneration/RootedTransientGraph.hpp"
#include "GraphGeneration/RootedTransientStandardizer.hpp"

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
		throw std::invalid_argument("wrong edge count in rooted transient test");
	}
	GraphType graph;
	Int edge = 0;
	for (const auto [first, second] : edges) {
		graph.setEdge(edge++, first, second);
	}
	return graph;
}

template <typename GraphType>
using unsigned_standardizer = GraphStandardizer<
	GraphType::N_VERTICES_,
	GraphType::N_EDGES_,
	GraphType::N_OUT_HAIR_,
	GraphType::N_IN_HAIR_,
	GraphType::C_,
	GraphType::D_,
	typename GraphType::Field
>;

template <typename RootedGraph>
struct brute_force_result {
	RootedGraph canonical_graph;
	std::uint64_t automorphism_count = 0;
	std::uint64_t permutation_count = 0;
};

template <typename RootedGraph>
brute_force_result<RootedGraph> brute_force_standardization(
	const RootedGraph& input
) {
	typename RootedGraph::nonroot_permutation permutation{};
	std::iota(permutation.begin(), permutation.end(), Int{0});

	brute_force_result<RootedGraph> result;
	bool have_candidate = false;
	do {
		RootedGraph candidate = input.permuted_nonroot(permutation);
		++result.permutation_count;
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

template <typename RootedGraph>
void check_standardizer_against_brute_force(const RootedGraph& graph) {
	using Standardizer = GraphGeneration::rooted_transient_standardizer<
		RootedGraph::N_VERTICES_,
		RootedGraph::N_EDGES_,
		typename RootedGraph::field_type
	>;
	Standardizer standardizer;
	const auto exhaustive_orbit = brute_force_standardization(graph);
	const auto actual = standardizer.standardize_with_info(graph);

	// Standardize4 selects the maximal refinement trace and only then compares
	// its final survivors.  Its representative therefore need not be the global
	// lexicographic maximum of the entire orbit.  The exhaustive maximum is an
	// independently chosen orbit member and must canonicalize identically.
	check(standardizer.standardize_no_sign(
		exhaustive_orbit.canonical_graph
	) == actual.canonical_graph,
		"rooted standardization is constant on the exhaustive permutation orbit");
	check(actual.automorphism_count == exhaustive_orbit.automorphism_count,
		"rooted standardizer returns the exact root-fixing automorphism count");
	check(actual.final_attempt_count >= actual.automorphism_count,
		"final refinement attempts contain every rooted automorphism");
	check(standardizer.standardize_no_sign(actual.canonical_graph)
		== actual.canonical_graph,
		"rooted standardization is idempotent");

	typename RootedGraph::nonroot_permutation permutation{};
	std::iota(permutation.begin(), permutation.end(), Int{0});
	std::uint64_t canonical_preimages = 0;
	do {
		const RootedGraph relabeled = graph.permuted_nonroot(permutation);
		canonical_preimages += relabeled == actual.canonical_graph;
		check(standardizer.standardize_no_sign(relabeled)
			== actual.canonical_graph,
			"every nonroot relabeling has the same rooted representative");
	} while (std::next_permutation(permutation.begin(), permutation.end()));
	check(canonical_preimages == actual.automorphism_count,
		"exhaustive orbit fibers confirm the rooted automorphism count");
}

template <typename RootedGraph>
void check_analysis_matches_legacy(const RootedGraph& graph) {
	const auto rooted = graph.analyze();
	const auto legacy = graph.to_transient().analyze();
	check(rooted.valences == legacy.valences,
		"compact and hairy ordinary valences agree");
	check(rooted.reduced_valences == legacy.reduced_valences,
		"compact and hairy reduced valences agree");
	check(rooted.maximum_valence == legacy.maximum_valence
		&& rooted.maximum_reduced_valence == legacy.maximum_reduced_valence
		&& rooted.maximum_valence_count == legacy.maximum_valence_count,
		"compact and hairy valence extrema agree");
	check(rooted.connected == legacy.connected
		&& rooted.at_least_trivalent == legacy.at_least_trivalent
		&& rooted.simple == legacy.simple,
		"compact and hairy basic classifications agree");
	check(rooted.root_is_maximum == legacy.root_is_maximum
		&& rooted.root_is_unique_maximum == legacy.root_is_unique_maximum
		&& rooted.root_is_reduced_maximum == legacy.root_is_reduced_maximum,
		"compact and hairy root conditions agree");
	check(graph.classification() == graph.to_transient().classification(),
		"compact and hairy transient classifications agree");
}

template <typename LegacyTransient>
void check_conversion_and_old_canonical_oracle(
	const LegacyTransient& legacy
) {
	using Rooted = GraphGeneration::rooted_transient_graph<
		LegacyTransient::N_VERTICES_,
		LegacyTransient::N_EDGES_,
		typename LegacyTransient::graph_type::Field
	>;
	using Hairy = typename LegacyTransient::graph_type;
	using RootedStandardizer = GraphGeneration::rooted_transient_standardizer<
		Rooted::N_VERTICES_, Rooted::N_EDGES_, typename Rooted::field_type
	>;

	const Rooted compact(legacy);
	compact.validate();
	check(Rooted(compact.to_transient()) == compact,
		"compact-to-hairy conversion round trips exactly");
	check_analysis_matches_legacy(compact);

	unsigned_standardizer<Hairy> old_standardizer;
	RootedStandardizer rooted_standardizer;
	const Hairy old_canonical = old_standardizer.standardize_no_sign(
		legacy.graph()
	);
	const Rooted compact_canonical
		= rooted_standardizer.standardize_no_sign(compact);
	const Hairy expanded_canonical = old_standardizer.standardize_no_sign(
		compact_canonical.to_hairy_graph()
	);
	check(expanded_canonical == old_canonical,
		"compact canonicalization preserves the legacy rooted isomorphism class");
}

void test_named_graphs() {
	{
		using Rooted = GraphGeneration::rooted_transient_graph<1, 4>;
		Rooted rose;
		rose.set_tadpole_count(4);
		rose.validate();
		check(rose.tadpole_count() == 4 && rose.represented_edge_count() == 4,
			"the compact rose stores all edges as root tadpoles");
		check_standardizer_against_brute_force(rose);
		check_conversion_and_old_canonical_oracle(rose.to_transient());
	}

	{
		using Legacy = GraphGeneration::transient_graph<4, 6>;
		using Hairless = typename Legacy::hairless_graph_type;
		const Legacy rooted_k4(graph_with_edges<Hairless>({
			{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
		}), 2);
		using Rooted = GraphGeneration::rooted_transient_graph<4, 6>;
		using Standardizer
			= GraphGeneration::rooted_transient_standardizer<4, 6>;
		const Rooted compact(rooted_k4);
		const auto result = Standardizer{}.standardize_with_info(compact);
		check(result.automorphism_count == 6,
			"K4 with its root fixed has six vertex automorphisms");
		check_standardizer_against_brute_force(compact);
		check_conversion_and_old_canonical_oracle(rooted_k4);
	}

	{
		using Legacy = GraphGeneration::transient_graph<5, 8>;
		using Hairless = typename Legacy::hairless_graph_type;
		const Legacy rooted_wheel(graph_with_edges<Hairless>({
			{0, 1}, {0, 2}, {0, 3}, {0, 4},
			{1, 2}, {2, 3}, {3, 4}, {4, 1}
		}), 0);
		using Rooted = GraphGeneration::rooted_transient_graph<5, 8>;
		using Standardizer
			= GraphGeneration::rooted_transient_standardizer<5, 8>;
		const Rooted compact(rooted_wheel);
		const auto result = Standardizer{}.standardize_with_info(compact);
		check(result.automorphism_count == 8,
			"the rooted four-spoke wheel retains its dihedral automorphisms");
		check_standardizer_against_brute_force(compact);
		check_conversion_and_old_canonical_oracle(rooted_wheel);
	}

	{
		using Legacy = GraphGeneration::transient_graph<3, 6>;
		using Hairless = typename Legacy::hairless_graph_type;
		const Legacy weighted_root(graph_with_edges<Hairless>({
			{2, 2}, {2, 0}, {2, 0}, {2, 1}, {2, 1}, {0, 1}
		}), 2);
		using Rooted = GraphGeneration::rooted_transient_graph<3, 6>;
		const Rooted compact(weighted_root);
		check(compact.tadpole_count() == 1
			&& compact.root_multiplicity(1) == 2
			&& compact.root_multiplicity(2) == 2
			&& compact.has_core_edge(1, 2),
			"conversion moves an arbitrarily labeled source root to vertex zero");
		check_standardizer_against_brute_force(compact);
		check_conversion_and_old_canonical_oracle(weighted_root);
	}
}

void test_exhaustive_small_weighted_cores() {
	using Rooted = GraphGeneration::rooted_transient_graph<4, 5>;
	constexpr std::array<std::pair<Int, Int>, 3> core_edges{{
		{1, 2}, {1, 3}, {2, 3}
	}};

	std::uint64_t graph_count = 0;
	for (unsigned mask = 0; mask < (1U << core_edges.size()); ++mask) {
		const Int core_count = static_cast<Int>(std::popcount(mask));
		const Int root_edge_count = static_cast<Int>(Rooted::N_EDGES_ - core_count);
		for (Int tadpoles = 0; tadpoles <= root_edge_count; ++tadpoles) {
			for (Int first = 0;
				 first <= static_cast<Int>(root_edge_count - tadpoles);
				 ++first) {
				for (Int second = 0;
					 second <= static_cast<Int>(root_edge_count - tadpoles - first);
					 ++second) {
					const Int third = static_cast<Int>(
						root_edge_count - tadpoles - first - second
					);
					Rooted graph;
					graph.set_tadpole_count(tadpoles);
					graph.set_root_multiplicity(1, first);
					graph.set_root_multiplicity(2, second);
					graph.set_root_multiplicity(3, third);
					for (std::size_t bit = 0; bit < core_edges.size(); ++bit) {
						if ((mask >> bit) & 1U) {
							graph.set_core_edge(
								core_edges[bit].first, core_edges[bit].second
							);
						}
					}
					graph.validate();
					check_standardizer_against_brute_force(graph);
					check_analysis_matches_legacy(graph);

					using Hairy = typename Rooted::hairy_graph_type;
					unsigned_standardizer<Hairy> old_standardizer;
					GraphGeneration::rooted_transient_standardizer<4, 5>
						rooted_standardizer;
					const Hairy old_canonical = old_standardizer.standardize_no_sign(
						graph.to_hairy_graph()
					);
					const Hairy compact_canonical = old_standardizer.standardize_no_sign(
						rooted_standardizer.standardize_no_sign(graph).to_hairy_graph()
					);
					check(old_canonical == compact_canonical,
						"exhaustive compact representative agrees with the hairy oracle");
					++graph_count;
				}
			}
		}
	}
	check(graph_count == 231,
		"the exhaustive V=4,E=5 weighted-core fixture count is stable");
}

template <typename LegacyTransient>
void check_relevant_split_sets_match(const LegacyTransient& legacy_parent) {
	using RootedParent = GraphGeneration::rooted_transient_graph<
		LegacyTransient::N_VERTICES_,
		LegacyTransient::N_EDGES_,
		typename LegacyTransient::graph_type::Field
	>;
	using RootedChild = typename RootedParent::split_type;
	using Standardizer = GraphGeneration::rooted_transient_standardizer<
		RootedChild::N_VERTICES_,
		RootedChild::N_EDGES_,
		typename RootedChild::field_type
	>;
	using Classification = GraphGeneration::transient_classification;
	using ClassMap = std::unordered_map<RootedChild, Classification>;

	Standardizer standardizer;
	ClassMap legacy_classes;
	ClassMap compact_classes;
	bool classification_conflict = false;
	const auto insert = [&](ClassMap& classes, RootedChild child,
		Classification classification) {
		child = standardizer.standardize_no_sign(child);
		const auto [iterator, inserted]
			= classes.try_emplace(std::move(child), classification);
		classification_conflict |= !inserted
			&& iterator->second != classification;
	};

	const std::uint64_t legacy_callbacks
		= legacy_parent.for_each_relevant_root_split(
			[&](typename LegacyTransient::split_type child,
				Classification classification) {
				insert(legacy_classes, RootedChild(child), classification);
			}
		);
	const RootedParent compact_parent(legacy_parent);
	const std::uint64_t compact_callbacks
		= compact_parent.for_each_relevant_root_split(
			[&](RootedChild child, Classification classification) {
				insert(compact_classes, std::move(child), classification);
			}
		);

	check(!classification_conflict,
		"duplicate split representatives have one invariant classification");
	check(compact_classes == legacy_classes,
		"compact splitting preserves every relevant rooted isomorphism class");
	check(compact_callbacks == legacy_callbacks,
		"compact and legacy grouped splitting emit the same number of children");
}

void test_split_equivalence() {
	check_relevant_split_sets_match(
		GraphGeneration::transient_graph<1, 10>{}
	);

	using WeightedLegacy = GraphGeneration::transient_graph<3, 6>;
	using WeightedHairless = typename WeightedLegacy::hairless_graph_type;
	check_relevant_split_sets_match(WeightedLegacy(
		graph_with_edges<WeightedHairless>({
			{2, 2}, {2, 0}, {2, 0}, {2, 1}, {2, 1}, {0, 1}
		}), 2
	));

	using BundleLegacy = GraphGeneration::transient_graph<5, 9>;
	using BundleHairless = typename BundleLegacy::hairless_graph_type;
	check_relevant_split_sets_match(BundleLegacy(
		graph_with_edges<BundleHairless>({
			{0, 1}, {0, 2}, {0, 2}, {0, 3}, {0, 4},
			{1, 2}, {1, 4}, {2, 3}, {3, 4}
		}), 0
	));
}

void test_off_root_defects_are_rejected() {
	using Legacy = GraphGeneration::transient_graph<3, 7>;
	using Hairless = typename Legacy::hairless_graph_type;
	const Legacy off_root_loop(graph_with_edges<Hairless>({
		{0, 1}, {0, 1}, {0, 1}, {0, 2}, {0, 2}, {0, 2}, {1, 1}
	}), 0);
	bool rejected = false;
	try {
		(void)GraphGeneration::rooted_transient_graph<3, 7>(off_root_loop);
	} catch (const std::invalid_argument&) {
		rejected = true;
	}
	check(rejected,
		"the compact representation rejects tadpoles away from its root");
}

} // namespace

int main() {
	test_named_graphs();
	test_exhaustive_small_weighted_cores();
	test_split_equivalence();
	test_off_root_defects_are_rejected();

	if (failures != 0) {
		std::cerr << failures << " rooted transient test(s) failed\n";
		return 1;
	}
	std::cout << "all rooted transient tests passed\n";
	return 0;
}
