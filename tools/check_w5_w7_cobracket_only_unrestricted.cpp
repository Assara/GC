#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "coalgebra_search.hpp"
#include "coalgebra_utils.hpp"
#include "examplegraphs.hpp"

namespace {

using W5GC = OddGCdegZero<6>;
using W7GC = OddGCdegZero<8>;
using SearchGC = coalgebra_utils::GraInsertLowValentGC<W5GC, W7GC>;
using GraphType = SearchGC::GraphType;
using WedgeType = coalgebra_search::TwoFactorWedge<SearchGC, W5GC, W7GC>;
using WedgeLinComb = typename WedgeType::LinComb;
using L = SearchGC::L;

W5GC build_w5() {
	std::vector<bool> Lseq{false};
	return W5GC(V_graph<5>(Lseq));
}

W7GC build_w7() {
	std::vector<bool> LL{false, false};
	std::vector<bool> LR{false, true};
	auto w7 = W7GC(V_graph<7>(LL));
	w7 += W7GC(V_graph<7>(LR));
	w7.standardize_all();
	w7.sort_elements();
	return w7;
}

std::unordered_set<GraphType> build_unrestricted_lie_support(const W5GC& w5, const W7GC& w7) {
	const auto lie = coalgebra_utils::gra_lie_unrestricted(w5, w7);
	std::unordered_set<GraphType> support;
	support.reserve(lie.data().size());
	for (const auto& be : lie.data()) {
		support.insert(be.getValue());
	}
	return support;
}

std::unordered_map<GraphType, WedgeLinComb> build_cobracket_map(const std::unordered_set<GraphType>& support) {
	std::unordered_map<GraphType, WedgeLinComb> map;
	map.reserve(support.size());
	std::size_t processed = 0;
	for (const auto& graph : support) {
		auto cobracket = coalgebra_search::cobracket_of_gc<SearchGC, W5GC, W7GC>(SearchGC(graph, AssumeBasisOrderTag{}));
		cobracket.standardize_and_sort();
		map.emplace(graph, std::move(cobracket));
		++processed;
		if (processed % 50 == 0) {
			std::cout << "built cobracket rows: " << processed << " / " << support.size() << std::endl;
		}
	}
	return map;
}

SearchGC make_solution_gc(L solution) {
	solution.standardize_and_sort();
	return SearchGC(solution);
}

} // namespace

int main() {
	const auto w5 = build_w5();
	const auto w7 = build_w7();
	const auto support = build_unrestricted_lie_support(w5, w7);
	std::cout << "unrestricted lie support size: " << support.size() << std::endl;

	const auto target = coalgebra_search::wedge_target<W5GC, W7GC, WedgeType>(w5, w7);
	std::cout << "target wedge terms: " << target.size() << std::endl;

	const auto cobracket_map = build_cobracket_map(support);
	std::cout << "built cobracket map with " << cobracket_map.size() << " domain graphs" << std::endl;

	VectorSpace::wiedemann_primitive_finder<WedgeType, GraphType, fieldType> solver(cobracket_map);
	auto solution_opt = solver.find_primitive_or_empty(target);
	if (solution_opt.has_value() == false) {
		std::cout << "no solution to Delta C = W5 ^ W7 in unrestricted lie support" << std::endl;
		return 0;
	}

	auto solution_gc = make_solution_gc(*solution_opt);
	auto actual = coalgebra_search::cobracket_of_gc<SearchGC, W5GC, W7GC>(solution_gc);
	std::cout << "found solution to Delta-side only" << std::endl;
	std::cout << "solution support size: " << solution_gc.data().size() << std::endl;
	std::cout << "actual cobracket equals target: " << ((actual == target) ? "yes" : "no") << std::endl;

	auto d_solution = solution_gc.d_contraction();
	d_solution.standardize_all();
	d_solution.sort_elements();
	std::cout << "d(solution) terms: " << d_solution.data().size() << std::endl;

	return 0;
}
