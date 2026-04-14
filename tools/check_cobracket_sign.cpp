#include <iostream>
#include <vector>

#include "DegreeZeroGraphWedge.hpp"
#include "examplegraphs.hpp"

namespace {

using FullGraph = OddGraphdegZero<13>;
using G5 = OddGraphdegZero<6>;
using G7 = OddGraphdegZero<8>;

bool check_single_factor_sign() {
	using OneWedge = DegreeZeroGraphWedge<48, 1>;
	using TwoWedge = DegreeZeroGraphWedge<48, 2>;

	std::vector<std::vector<bool>> blocks{{false}, {false, false}};
	FullGraph graph = iterated_V_graph_lcr_extension<12>(blocks);
	auto cobracket = OneWedge::from_graph(graph).cobracket();

	std::vector<bool> left{false};
	std::vector<bool> left_left{false, false};
	TwoWedge target = TwoWedge::from_graphs(V_graph<5>(left), V_graph<7>(left_left));
	BasisElement<TwoWedge, fieldType> be(target, fieldType{1});
	const auto canon_target = TwoWedge::canonized(be);

	if (cobracket.size() != 1) {
		std::cout << "single-factor check failed: expected 1 cobracket term, got "
			  << cobracket.size() << std::endl;
		return false;
	}

	const auto& term = cobracket.raw_elements().front();
	if (!(term.getValue() == canon_target.getValue())) {
		std::cout << "single-factor check failed: wedge value mismatch" << std::endl;
		return false;
	}
	if (term.getCoefficient() != canon_target.getCoefficient()) {
		std::cout << "single-factor check failed: sign mismatch" << std::endl;
		return false;
	}

	return true;
}

bool check_two_factor_antisymmetry() {
	using TwoWedge = DegreeZeroGraphWedge<68, 2>;

	std::vector<bool> left{false};
	G5 x = V_graph<5>(left);

	std::vector<std::vector<bool>> blocks{{false}, {false, false}};
	FullGraph g = iterated_V_graph_lcr_extension<12>(blocks);

	BasisElement<TwoWedge, fieldType> be_xg(TwoWedge::from_graphs(x, g), fieldType{1});
	BasisElement<TwoWedge, fieldType> be_gx(TwoWedge::from_graphs(g, x), fieldType{1});
	const auto canon_xg = TwoWedge::canonized(be_xg);
	const auto canon_gx = TwoWedge::canonized(be_gx);

	if (!(canon_xg.getValue() == canon_gx.getValue())) {
		std::cout << "two-factor check failed: canonical wedge values differ" << std::endl;
		return false;
	}
	if (canon_xg.getCoefficient() != fieldType{-1} * canon_gx.getCoefficient()) {
		std::cout << "two-factor check failed: wedge antisymmetry sign mismatch" << std::endl;
		return false;
	}

	auto d_xg = canon_xg.getValue().cobracket();
	d_xg.scalar_multiply(canon_xg.getCoefficient());
	d_xg.standardize_and_sort();

	auto d_gx = canon_gx.getValue().cobracket();
	d_gx.scalar_multiply(canon_gx.getCoefficient());
	d_gx.standardize_and_sort();

	auto neg_d_xg = d_xg;
	neg_d_xg.scalar_multiply(fieldType{-1});
	neg_d_xg.standardize_and_sort();

	if (!(d_gx == neg_d_xg)) {
		std::cout << "two-factor check failed: cobracket is not antisymmetric under wedge swap" << std::endl;
		return false;
	}

	return true;
}

} // namespace

int main() {
	const bool ok_single = check_single_factor_sign();
	const bool ok_two_factor = check_two_factor_antisymmetry();

	std::cout << "single-factor sign check: " << (ok_single ? "pass" : "fail") << std::endl;
	std::cout << "two-factor antisymmetry check: " << (ok_two_factor ? "pass" : "fail") << std::endl;

	return (ok_single && ok_two_factor) ? 0 : 1;
}
