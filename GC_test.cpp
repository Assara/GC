#include <cstdlib>
#include <iostream>

#include "examplegraphs.hpp"
#include "OCGraph_test.hpp"

template <typename GraphType>
GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_> make_graph_isomorphism_from_vertex_perm(
	const GraphType& graph,
	const std::array<Int, GraphType::N_VERTICES_>& vertex_perm
) {
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	IsoType iso;
	iso.vertex_perm = vertex_perm;

	for (Int source_edge = 0; source_edge < GraphType::N_EDGES_; ++source_edge) {
		auto [u, v] = graph.getEdge(source_edge);
		u = vertex_perm[u];
		v = vertex_perm[v];

		const Int target_edge = graph.find_edge_index(u, v);
		iso.edge_perm[source_edge] = target_edge;

		const auto [tu, tv] = graph.getEdge(target_edge);
		iso.edge_flip[source_edge] = (tu != u || tv != v);
	}

	return iso;
}

template <Int N>
bool check_wheel_ocgraph_delta_e_isomorphism_equivariance(const char* label) {
	using GraphType = OddGraphdegZero<N + 1>;
	using DirectionType = GraphDirections<GraphType>;
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	const GraphType source = wheel_graph<N>();

	DirectionType directions;
	directions.fill(true);

	std::vector<IsoType> isomorphisms;

	std::array<Int, GraphType::N_VERTICES_> rotation{};
	rotation[0] = 0;
	for (Int v = 1; v < GraphType::N_VERTICES_; ++v) {
		rotation[v] = (v == N) ? 1 : static_cast<Int>(v + 1);
	}
	isomorphisms.push_back(make_graph_isomorphism_from_vertex_perm(source, rotation));

	std::array<Int, GraphType::N_VERTICES_> reflection{};
	reflection[0] = 0;
	reflection[1] = 1;
	for (Int v = 2; v < GraphType::N_VERTICES_; ++v) {
		reflection[v] = static_cast<Int>(GraphType::N_VERTICES_ + 1 - v);
	}
	isomorphisms.push_back(make_graph_isomorphism_from_vertex_perm(source, reflection));

	return check_ocgraph_delta_e_isomorphism_equivariance_on_data(label, source, directions, isomorphisms);
}

template <typename GraphType>
std::vector<GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>> sample_graph_isomorphisms(
	const GraphType&
) {
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	std::vector<IsoType> isomorphisms;

	IsoType a;
	for (Int v = 0; v < GraphType::N_VERTICES_; ++v) {
		a.vertex_perm[v] = static_cast<Int>((v + 1) % GraphType::N_VERTICES_);
	}
	for (Int e = 0; e < GraphType::N_EDGES_; ++e) {
		a.edge_perm[e] = static_cast<Int>((e + 1) % GraphType::N_EDGES_);
		a.edge_flip[e] = (e % 2) == 0;
	}
	isomorphisms.push_back(a);

	IsoType b;
	for (Int v = 0; v < GraphType::N_VERTICES_; ++v) {
		b.vertex_perm[v] = static_cast<Int>((GraphType::N_VERTICES_ - 1 - v));
	}
	for (Int e = 0; e < GraphType::N_EDGES_; ++e) {
		b.edge_perm[e] = static_cast<Int>((GraphType::N_EDGES_ - 1 - e));
		b.edge_flip[e] = (e % 3) == 1;
	}
	isomorphisms.push_back(b);

	IsoType c;
	for (Int v = 0; v < GraphType::N_VERTICES_; ++v) {
		c.vertex_perm[v] = static_cast<Int>((v + 2) % GraphType::N_VERTICES_);
	}
	for (Int e = 0; e < GraphType::N_EDGES_; ++e) {
		c.edge_perm[e] = static_cast<Int>((e + 3) % GraphType::N_EDGES_);
		c.edge_flip[e] = (e % 2) == 1;
	}
	isomorphisms.push_back(c);

	return isomorphisms;
}

template <Int N>
GraphIsomorphism<OddGraphdegZero<N + 1>::N_VERTICES_, OddGraphdegZero<N + 1>::N_EDGES_>
wheel_swap_spoke_and_rim_labels_isomorphism() {
	using GraphType = OddGraphdegZero<N + 1>;
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	IsoType iso;
	for (Int edge = 0; edge < N; ++edge) {
		iso.edge_perm[edge] = static_cast<Int>(N + edge);
		iso.edge_perm[N + edge] = edge;
	}
	return iso;
}

template <Int N>
std::vector<GraphIsomorphism<OddGraphdegZero<N + 1>::N_VERTICES_, OddGraphdegZero<N + 1>::N_EDGES_>>
wheel_OCG_F0_test_isomorphisms() {
	using GraphType = OddGraphdegZero<N + 1>;

	auto isomorphisms = sample_graph_isomorphisms<GraphType>(wheel_graph<N>());
	isomorphisms.push_back(wheel_swap_spoke_and_rim_labels_isomorphism<N>());
	return isomorphisms;
}

template <Int N>
bool check_wheel_OCG_F0_exact_equivariance(const char* label) {
	using GraphType = OddGraphdegZero<N + 1>;
	const GraphType source = wheel_graph<N>();
	return check_OCG_F0_exact_equivariance_on_graph(label, source, wheel_OCG_F0_test_isomorphisms<N>());
}

bool check_w9_OCG_F0_unique_direction_data(const char* label) {
	using GraphType = OddGraphdegZero<10>;
	const GraphType source = wheel_graph<9>();
	return check_OCG_F0_unique_direction_on_graph(label, source, wheel_OCG_F0_test_isomorphisms<9>());
}

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
	ok &= check_graph_isomorphism_graph_sign<OddLoopGraphType<3>>("triangle_iso_graph_sign");
	ok &= check_graph_directions_shape<OddLoopGraphType<3>>("triangle_directions");
	ok &= check_graph_directions_isomorphism_action<OddLoopGraphType<3>>("triangle_directions_iso");
	ok &= check_graph_directions_order<OddLoopGraphType<3>>("triangle_directions_order");
	ok &= check_graph_directions_standardize<OddLoopGraphType<3>>("triangle_directions_std");
	ok &= check_ocgraph_standardize_and_sort<OddLoopGraphType<3>>("triangle_ocgraph_std");
	ok &= check_ocgraph_filters_non_covering_directions<OddLoopGraphType<3>>("triangle_ocgraph_filter");
	ok &= check_ocgraph_total_sign_cancellation("double_edge_ocgraph_total_sign_cancel");
	ok &= check_ocgraph_OCG_F0_shifted_degree("triangle_ocgraph_OCG_F0_shift");
	ok &= check_wheel_ocgraph_delta_e_isomorphism_equivariance<5>("wheel5_ocgraph_delta_e_iso");
	ok &= check_wheel_ocgraph_delta_e_isomorphism_equivariance<7>("wheel7_ocgraph_delta_e_iso");
	ok &= check_wheel_OCG_F0_exact_equivariance<5>("wheel5_OCG_F0_exact");
	ok &= check_wheel_OCG_F0_exact_equivariance<7>("wheel7_OCG_F0_exact");
	ok &= check_wheel_OCG_F0_exact_equivariance<9>("wheel9_OCG_F0_exact");
	ok &= check_w9_OCG_F0_unique_direction_data("wheel9_OCG_F0_direction_data");
	ok &= OCGraph_test<OddLoopGraphType<3>>::object();
	ok &= check_minimizing_isomorphisms("triangle_minimizers", loop_graph<3>(), 6);
	ok &= check_minimizing_isomorphisms_match_standardization("wheel_11_minimizers", wheel_graph<11>(), 22);

	if (!ok) {
		return EXIT_FAILURE;
	}

	std::cout << "all tests passed\n";
	return EXIT_SUCCESS;
}
