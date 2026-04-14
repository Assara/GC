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
using ContL = SearchGC::ContL;
using L = SearchGC::L;
using WedgeType = coalgebra_search::TwoFactorWedge<SearchGC, W5GC, W7GC>;

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

std::unordered_set<GraphType> build_lie_support(const SearchGC& lie) {
	std::unordered_set<GraphType> support;
	support.reserve(lie.data().size());
	for (const auto& be : lie.data()) {
		support.insert(be.getValue());
	}
	return support;
}

std::unordered_map<GraphType, ContL> build_restricted_differential_map(
	const std::unordered_set<GraphType>& support
) {
	std::unordered_map<GraphType, ContL> differential_map;
	differential_map.reserve(support.size());

	for (const auto& graph : support) {
		ContL d_graph = SearchGC(graph, AssumeBasisOrderTag{}).d_contraction().data();
		d_graph.standardize_all();
		d_graph.sort_elements();
		differential_map.emplace(graph, std::move(d_graph));
	}

	return differential_map;
}

SearchGC make_cycle_from_solution(const GraphType& pivot_graph, L solution) {
	solution.append_in_basis_order(pivot_graph, fieldType{1});
	solution.standardize_and_sort();
	return SearchGC(solution);
}

template <typename LinCombType>
std::size_t overlap_count(const LinCombType& left, const LinCombType& right) {
	const auto& a = left.raw_elements();
	const auto& b = right.raw_elements();
	std::size_t i = 0;
	std::size_t j = 0;
	std::size_t overlap = 0;

	while (i < a.size() && j < b.size()) {
		const auto cmp = a[i].getValue().compare(b[j].getValue());
		if (cmp < 0) {
			++i;
		} else if (cmp > 0) {
			++j;
		} else {
			++overlap;
			++i;
			++j;
		}
	}

	return overlap;
}

} // namespace

int main() {
	const auto w5 = build_w5();
	const auto w7 = build_w7();
	const auto lie = coalgebra_utils::gra_lie_unrestricted(w5, w7);
	const auto support = build_lie_support(lie);
	const auto differential_map = build_restricted_differential_map(support);

	std::vector<GraphType> ordered_support(support.begin(), support.end());
	std::sort(ordered_support.begin(), ordered_support.end());

	std::cout << "unrestricted lie support size: " << ordered_support.size() << std::endl;

	std::size_t pivot_index = 0;
	for (const auto& pivot_graph : ordered_support) {
		auto target_it = differential_map.find(pivot_graph);
		if (target_it == differential_map.end()) {
			++pivot_index;
			continue;
		}

		ContL target = target_it->second;
		target.scalar_multiply(fieldType{-1});

		std::unordered_map<GraphType, ContL> restricted_map;
		restricted_map.reserve(differential_map.size() - 1);
		for (const auto& [graph, differential] : differential_map) {
			if (graph == pivot_graph) {
				continue;
			}
			restricted_map.emplace(graph, differential);
		}

		std::cout << "trying pivot " << pivot_index << " / " << ordered_support.size() << std::endl;
		VectorSpace::wiedemann_primitive_finder<SearchGC::ContGraphType, GraphType, fieldType> solver(restricted_map);
		auto solution_opt = solver.find_primitive_or_empty(target);
		if (!solution_opt.has_value()) {
			++pivot_index;
			continue;
		}

		auto cycle = make_cycle_from_solution(pivot_graph, *solution_opt);
		auto d_cycle = cycle.d_contraction();
		d_cycle.standardize_all();
		d_cycle.sort_elements();

		if (d_cycle.size() != 0) {
			std::cout << "candidate from pivot " << pivot_index
				  << " did not verify; d(cycle) terms = " << d_cycle.size() << std::endl;
			++pivot_index;
			continue;
		}

		const auto cycle_cobracket = coalgebra_search::cobracket_of_gc<SearchGC, W5GC, W7GC>(cycle);
		const auto target_wedge = coalgebra_search::wedge_target<W5GC, W7GC, WedgeType>(w5, w7);

		std::cout << "found unrestricted lie-support cycle" << std::endl;
		std::cout << "pivot index: " << pivot_index << std::endl;
		std::cout << "cycle support size: " << cycle.data().size() << std::endl;
		std::cout << "cobracket support size: " << cycle_cobracket.size() << std::endl;
		std::cout << "target wedge support size: " << target_wedge.size() << std::endl;
		std::cout << "cobracket equals target wedge: "
			  << ((cycle_cobracket == target_wedge) ? "yes" : "no") << std::endl;
		std::cout << "common wedge terms with target: "
			  << overlap_count(cycle_cobracket, target_wedge) << std::endl;

		for (const auto& be : cycle_cobracket) {
			std::cout << "cobracket coeff: " << be.getCoefficient() << std::endl;
		}
		for (const auto& be : target_wedge) {
			std::cout << "target coeff: " << be.getCoefficient() << std::endl;
		}

		return 0;
	}

	std::cout << "no unrestricted lie-support cycle found" << std::endl;
	return 0;
}
