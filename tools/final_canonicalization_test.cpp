#include <array>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>

#include "GraphGeneration/FinalCanonicalization.hpp"

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
		throw std::invalid_argument("wrong edge count in canonicalization test");
	}
	GraphType graph;
	Int edge_index = 0;
	for (const auto [first, second] : edges) {
		graph.setEdge(edge_index++, first, second);
	}
	return graph;
}

template <signedInt C, signedInt D, typename GraphType>
bool signed_standardization_survives(const GraphType& graph) {
	using SignedGraph = typename GraphType::template RebindDegree<C, D>;
	using Standardizer = GraphStandardizer<
		SignedGraph::N_VERTICES_,
		SignedGraph::N_EDGES_,
		0,
		0,
		C,
		D,
		typename SignedGraph::Field
	>;
	SignedGraph signed_graph;
	signed_graph.half_edges = graph.half_edges;
	const typename SignedGraph::Basis input(signed_graph);
	return Standardizer{}.standardize4(input).getCoefficient()
		!= typename SignedGraph::Field{};
}

template <typename GraphType>
void check_against_signed_standardizers(
	const GraphType& graph,
	const GraphGeneration::final_canonicalization_result<GraphType>& result,
	std::string_view label
) {
	check(
		result.survives_odd_edges()
			== signed_standardization_survives<0, 1>(graph),
		std::string(label).append(" odd-edge flag matches signed standardization")
	);
	check(
		result.survives_even_edges_odd_vertices()
			== signed_standardization_survives<0, 0>(graph),
		std::string(label).append(" odd-vertex flag matches signed standardization")
	);
}

void test_k4() {
	using GraphType = Graph<4, 6, 0, 0, 0, 0, fieldType>;
	const GraphType graph = graph_with_edges<GraphType>({
		{0, 1}, {0, 2}, {0, 3}, {1, 2}, {1, 3}, {2, 3}
	});
	const auto result = GraphGeneration::canonicalize_final_graph(graph);

	check(result.automorphism_order == 24, "K4 has automorphism order 24");
	check(result.survives_odd_edges(), "K4 survives the odd-edge convention");
	check(!result.survives_even_edges_odd_vertices(),
		"K4 vanishes when vertices are odd");
	check_against_signed_standardizers(graph, result, "K4");

	using Standardizer = GraphStandardizer<4, 6, 0, 0, 0, 0, fieldType>;
	check(result.canonical_graph == Standardizer{}.standardize_no_sign(graph),
		"K4 uses the unsigned standardizer's canonical representative");
}

void test_triangle() {
	using GraphType = Graph<3, 3, 0, 0, 0, 0, fieldType>;
	const GraphType graph = graph_with_edges<GraphType>({
		{0, 1}, {1, 2}, {2, 0}
	});
	const auto result = GraphGeneration::canonicalize_final_graph(graph);

	check(result.automorphism_order == 6, "the triangle has automorphism order 6");
	check(!result.survives_odd_edges(),
		"the triangle has an odd induced edge automorphism");
	check(!result.survives_even_edges_odd_vertices(),
		"the triangle has an odd vertex automorphism");
	check_against_signed_standardizers(graph, result, "triangle");
}

void test_asymmetric_graph_and_relabeling() {
	using GraphType = Graph<6, 6, 0, 0, 0, 0, fieldType>;
	GraphType graph = graph_with_edges<GraphType>({
		{0, 2}, {0, 3}, {0, 5}, {1, 2}, {1, 4}, {2, 3}
	});
	const auto result = GraphGeneration::canonicalize_final_graph(graph);

	check(result.automorphism_order == 1,
		"the asymmetric fixture has trivial automorphism group");
	check(result.survives_odd_edges() && result.survives_even_edges_odd_vertices(),
		"a graph with trivial automorphism group survives both conventions");
	check_against_signed_standardizers(graph, result, "asymmetric graph");

	const Permutation<6> relabel(std::array<Int, 6>({4, 1, 5, 0, 2, 3}));
	(void)graph.permuteVertices(relabel);
	const auto relabeled_result = GraphGeneration::canonicalize_final_graph(graph);
	check(relabeled_result.canonical_graph == result.canonical_graph,
		"canonical graph is invariant under vertex relabeling");
	check(relabeled_result.automorphism_order == result.automorphism_order,
		"automorphism order is invariant under vertex relabeling");
	check(relabeled_result.survival == result.survival,
		"both survival flags are invariant under vertex relabeling");
}

void test_nonsimple_graphs_are_rejected() {
	using GraphType = Graph<2, 2, 0, 0, 0, 0, fieldType>;
	const GraphType parallel = graph_with_edges<GraphType>({{0, 1}, {0, 1}});
	bool rejected_parallel = false;
	try {
		(void)GraphGeneration::canonicalize_final_graph(parallel);
	} catch (const std::invalid_argument&) {
		rejected_parallel = true;
	}
	check(rejected_parallel, "parallel edges are rejected before final canonicalization");

	const GraphType tadpole = graph_with_edges<GraphType>({{0, 0}, {0, 1}});
	bool rejected_tadpole = false;
	try {
		(void)GraphGeneration::canonicalize_final_graph(tadpole);
	} catch (const std::invalid_argument&) {
		rejected_tadpole = true;
	}
	check(rejected_tadpole, "tadpoles are rejected before final canonicalization");
}

} // namespace

int main() {
	test_k4();
	test_triangle();
	test_asymmetric_graph_and_relabeling();
	test_nonsimple_graphs_are_rejected();

	if (failures != 0) {
		std::cerr << failures << " final canonicalization test(s) failed\n";
		return 1;
	}
	std::cout << "final canonicalization tests passed\n";
	return 0;
}
