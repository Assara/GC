#include <cstdlib>
#include <iostream>

#include "examplegraphs.hpp"
#include "GraphDirections.hpp"
#include "GraphIsomorphism.hpp"

template <typename GCType>
bool check_odd_even_contraction_split(const GCType& input, const char* label) {
	GCType gc = input;
	auto d_total = gc.d_contraction();
	auto d_even = gc.d_even_contraction();
	auto d_odd = gc.d_odd_contraction();
	auto sum = d_even;
	sum += d_odd;
	sum.standardize_all();
	d_total.standardize_all();

	const bool ok = (sum.data() == d_total.data());
	std::cout << label
		  << ": total=" << d_total.size()
		  << " even=" << d_even.size()
		  << " odd=" << d_odd.size()
		  << " -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

static OddGraphdegZero<10> w9_term_1() {
	OddGraphdegZero<10> g;
	int e = 0;
	for (auto [u, v] : std::initializer_list<std::pair<int, int>>{
		{0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,5}, {2,4}, {2,6},
		{3,5}, {3,7}, {4,6}, {4,8}, {5,7}, {6,9}, {7,8}, {7,9}, {8,9}
	}) {
		g.setEdge(e++, u, v);
	}
	return g;
}

static OddGraphdegZero<10> w9_term_2() {
	OddGraphdegZero<10> g;
	int e = 0;
	for (auto [u, v] : std::initializer_list<std::pair<int, int>>{
		{0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,5}, {2,4}, {2,6},
		{3,5}, {4,7}, {5,6}, {5,8}, {6,8}, {6,9}, {7,8}, {7,9}, {8,9}
	}) {
		g.setEdge(e++, u, v);
	}
	return g;
}

static OddGraphdegZero<10> w9_term_3() {
	OddGraphdegZero<10> g;
	int e = 0;
	for (auto [u, v] : std::initializer_list<std::pair<int, int>>{
		{0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,5}, {2,4}, {2,6},
		{3,7}, {3,8}, {4,6}, {5,6}, {5,7}, {5,9}, {6,9}, {7,8}, {8,9}
	}) {
		g.setEdge(e++, u, v);
	}
	return g;
}

template <typename GraphType>
bool check_graph_isomorphism_action(const char* label) {
	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	GraphIsomorphism<3, 3> iso;
	iso.vertex_perm = {1, 2, 0};
	iso.edge_perm = {2, 0, 1};
	iso.edge_flip = {false, true, false};

	GraphType image = iso.permute(source);

	GraphType expected;
	expected.setEdge(0, 0, 2);
	expected.setEdge(1, 1, 0);
	expected.setEdge(2, 1, 2);

	const bool ok = (image == expected);
	std::cout << label << ": graph isomorphism action -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_graph_isomorphism_composition(const char* label) {
	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	GraphIsomorphism<3, 3> a;
	a.vertex_perm = {1, 2, 0};
	a.edge_perm = {2, 0, 1};
	a.edge_flip = {false, true, false};

	GraphIsomorphism<3, 3> b;
	b.vertex_perm = {2, 0, 1};
	b.edge_perm = {1, 2, 0};
	b.edge_flip = {true, false, true};

	const GraphType stepwise = b.permute(a.permute(source));
	const GraphIsomorphism<3, 3> ab = a.compose(b);
	const GraphType composed = ab.permute(source);

	const bool ok = (stepwise == composed);
	std::cout << label << ": graph isomorphism composition -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_graph_isomorphism_inverse(const char* label) {
	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	GraphIsomorphism<3, 3> iso;
	iso.vertex_perm = {1, 2, 0};
	iso.edge_perm = {2, 0, 1};
	iso.edge_flip = {false, true, false};

	const auto inv = iso.inverse();
	const GraphType recovered = inv.permute(iso.permute(source));
	const GraphType recovered_other_side = iso.permute(inv.permute(source));

	const GraphIsomorphism<3, 3> id_left = iso.compose(inv);
	const GraphIsomorphism<3, 3> id_right = inv.compose(iso);
	const GraphType identity_left_image = id_left.permute(source);
	const GraphType identity_right_image = id_right.permute(source);

	const bool ok =
		(recovered == source) &&
		(recovered_other_side == source) &&
		(identity_left_image == source) &&
		(identity_right_image == source);

	std::cout << label << ": graph isomorphism inverse -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_minimizing_isomorphisms(const char* label, const GraphType& source, bigInt expected_count) {
	GraphStandardizer<
		GraphType::N_VERTICES_,
		GraphType::N_EDGES_,
		GraphType::N_OUT_HAIR_,
		GraphType::N_IN_HAIR_,
		GraphType::C_,
		GraphType::D_,
		fieldType
	> standardizer;

	const auto minimizers = standardizer.minimizing_isomorphisms(source);
	bool ok = !minimizers.empty() && minimizers.size() == expected_count;

	if (!minimizers.empty()) {
		const GraphType target = minimizers.front().permute(source);
		for (const auto& iso : minimizers) {
			if (!(iso.permute(source) == target)) {
				ok = false;
				break;
			}
		}
	}

	const bool count_matches = standardizer.automorphism_group_size(source) == minimizers.size();
	ok &= count_matches;

	std::cout << label
		  << ": minimizing isomorphisms -> "
		  << (ok ? "ok" : "failed")
		  << " (count=" << minimizers.size() << ")\n";
	return ok;
}

template <typename GraphType>
bool check_minimizing_isomorphisms_match_standardization(const char* label, const GraphType& source, bigInt expected_count) {
	GraphStandardizer<
		GraphType::N_VERTICES_,
		GraphType::N_EDGES_,
		GraphType::N_OUT_HAIR_,
		GraphType::N_IN_HAIR_,
		GraphType::C_,
		GraphType::D_,
		fieldType
	> standardizer;

	const auto minimizers = standardizer.minimizing_isomorphisms(source);
	typename GraphType::Basis standardized_input(source, fieldType{1});
	GraphType::std(standardized_input);
	const GraphType standardized_graph = standardized_input.getValue();

	bool ok = minimizers.size() == expected_count;
	for (const auto& iso : minimizers) {
		if (!(iso.permute(source) == standardized_graph)) {
			ok = false;
			break;
		}
	}

	std::cout << label
		  << ": minimizing isomorphisms match std -> "
		  << (ok ? "ok" : "failed")
		  << " (count=" << minimizers.size() << ")\n";
	return ok;
}

template <typename GraphType>
bool check_graph_directions_shape(const char* label) {
	GraphDirections<GraphType> directions;
	directions.fill(false);
	directions[0] = true;
	directions[GraphType::SIZE - 1] = true;

	const bool ok =
		(directions.size() == GraphType::SIZE) &&
		directions[0] &&
		directions[GraphType::SIZE - 1];

	std::cout << label << ": graph directions shape -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_graph_directions_isomorphism_action(const char* label) {
	GraphDirections<GraphType> source;
	source.fill(false);
	source[0] = true;
	source[GraphType::N_HAIR + 0] = true;
	source[GraphType::N_HAIR + 3] = true;

	GraphIsomorphism<3, 3> iso;
	iso.vertex_perm = {1, 2, 0};
	iso.edge_perm = {2, 0, 1};
	iso.edge_flip = {false, true, false};

	const auto image = iso.permute(source);

	GraphDirections<GraphType> expected;
	expected.fill(false);
	expected[GraphType::N_HAIR + 0] = true;
	expected[GraphType::N_HAIR + 4] = true;

	const bool ok = (image == expected);
	std::cout << label << ": graph directions isomorphism action -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_graph_directions_order(const char* label) {
	GraphDirections<GraphType> a;
	GraphDirections<GraphType> b;
	GraphDirections<GraphType> c;

	a.fill(false);
	b.fill(false);
	c.fill(false);

	b[GraphType::SIZE - 1] = true;
	c[0] = true;

	const bool ok =
		(a.compare(a) == 0) &&
		(a < b) &&
		(b < c) &&
		(a < c) &&
		!(b < a) &&
		!(c < b);

	std::cout << label << ": graph directions order -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

int main() {
	bool ok = true;
	OddGCdegZero<10> combined(w9_term_1());
	combined += OddGCdegZero<10>(w9_term_2());
	combined += OddGCdegZero<10>(w9_term_3());

	ok &= check_odd_even_contraction_split(OddGCdegZero<10>(w9_term_1()), "w9_term_1");
	ok &= check_odd_even_contraction_split(OddGCdegZero<10>(w9_term_2()), "w9_term_2");
	ok &= check_odd_even_contraction_split(OddGCdegZero<10>(w9_term_3()), "w9_term_3");
	ok &= check_odd_even_contraction_split(combined, "w9_sum");
	ok &= check_graph_isomorphism_action<OddLoopGraphType<3>>("triangle_iso");
	ok &= check_graph_isomorphism_composition<OddLoopGraphType<3>>("triangle_iso_comp");
	ok &= check_graph_isomorphism_inverse<OddLoopGraphType<3>>("triangle_iso_inv");
	ok &= check_graph_directions_shape<OddLoopGraphType<3>>("triangle_directions");
	ok &= check_graph_directions_isomorphism_action<OddLoopGraphType<3>>("triangle_directions_iso");
	ok &= check_graph_directions_order<OddLoopGraphType<3>>("triangle_directions_order");
	ok &= check_minimizing_isomorphisms("triangle_minimizers", loop_graph<3>(), 6);
	ok &= check_minimizing_isomorphisms_match_standardization("wheel_11_minimizers", wheel_graph<11>(), 22);

	if (!ok) {
		return EXIT_FAILURE;
	}

	std::cout << "all tests passed\n";
	return EXIT_SUCCESS;
}
