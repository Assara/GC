#include <cstdlib>
#include <iostream>

#include "examplegraphs.hpp"

template <typename GraphType>
bool check_contract_preserve_order_matches(const GraphType& graph, const char* label) {
	typename GraphType::Basis standardized_input(graph, fieldType{1});
	GraphType::std(standardized_input);
	const GraphType standardized_graph = standardized_input.getValue();

	bool ok = true;
	for (Int i = 0; i < GraphType::N_EDGES_; ++i) {
		auto standard = standardized_graph.contract_edge(i, fieldType{1});
		auto preserve = standardized_graph.contract_preserve_order(i, fieldType{1});
		using ContGraph = typename GraphType::ContGraph;
		ContGraph::std(standard);
		ContGraph::std(preserve);

		const bool edge_ok =
			standard.getCoefficient() == preserve.getCoefficient() &&
			standard.getValue() == preserve.getValue();
		if (!edge_ok) {
			std::cout << label << ": contract_preserve_order mismatch at edge " << int(i) << '\n';
			ok = false;
		}
	}

	std::cout << label << ": preserve-order contraction -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
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
	ok &= check_contract_preserve_order_matches(w9_term_1(), "w9_term_1");
	ok &= check_contract_preserve_order_matches(w9_term_2(), "w9_term_2");
	ok &= check_contract_preserve_order_matches(w9_term_3(), "w9_term_3");
	ok &= check_contract_preserve_order_matches(wheel_graph<5>(), "wheel_5");

	if (!ok) {
		return EXIT_FAILURE;
	}

	std::cout << "all tests passed\n";
	return EXIT_SUCCESS;
}
