#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "coalgebra_utils.hpp"
#include "examplegraphs.hpp"

namespace {

using W5GC = OddGCdegZero<6>;
using W7GC = OddGCdegZero<8>;
using InsertGC = coalgebra_utils::GraInsertLowValentGC<W5GC, W7GC>;
using GraphType = InsertGC::GraphType;
using ContL = InsertGC::ContL;
using L = InsertGC::L;

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

std::unordered_set<GraphType> build_insertion_support(const W5GC& w5, const W7GC& w7) {
	std::unordered_set<GraphType> support;

	for (const auto& left_be : w5.data()) {
		for (const auto& right_be : w7.data()) {
			auto forward = coalgebra_utils::gra_insertion_support(left_be.getValue(), right_be.getValue());
			support.insert(forward.begin(), forward.end());

			auto backward = coalgebra_utils::gra_insertion_support(right_be.getValue(), left_be.getValue());
			support.insert(backward.begin(), backward.end());
		}
	}

	return support;
}

std::unordered_map<GraphType, ContL> build_restricted_differential_map(
	const std::unordered_set<GraphType>& support
) {
	std::unordered_map<GraphType, ContL> differential_map;
	differential_map.reserve(support.size());

	for (const auto& graph : support) {
		ContL d_graph = InsertGC(graph, AssumeBasisOrderTag{}).d_contraction().data();
		d_graph.standardize_all();
		d_graph.sort_elements();
		differential_map.emplace(graph, std::move(d_graph));
	}

	return differential_map;
}

InsertGC make_cycle_from_solution(const GraphType& pivot_graph, L solution) {
	solution.append_in_basis_order(pivot_graph, fieldType{1});
	solution.standardize_and_sort();
	return InsertGC(solution);
}

} // namespace

int main() {
	const auto w5 = build_w5();
	const auto w7 = build_w7();
	const auto support = build_insertion_support(w5, w7);
	const auto differential_map = build_restricted_differential_map(support);

	std::vector<GraphType> ordered_support(support.begin(), support.end());
	std::sort(ordered_support.begin(), ordered_support.end());

	std::cout << "standardised insertion support size: " << ordered_support.size() << '\n';

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

		std::cout << "trying pivot " << pivot_index << " / " << ordered_support.size() << '\n';
		VectorSpace::wiedemann_primitive_finder<InsertGC::ContGraphType, GraphType, fieldType> solver(restricted_map);
		auto solution_opt = solver.find_primitive_or_empty(target);
		if (!solution_opt.has_value()) {
			++pivot_index;
			continue;
		}

		auto cycle = make_cycle_from_solution(pivot_graph, *solution_opt);
		auto d_cycle = cycle.d_contraction();
		d_cycle.standardize_all();
		d_cycle.sort_elements();

		if (d_cycle.size() == 0) {
			std::cout << "found non-zero cycle in standardised insertion support\n";
			std::cout << "pivot index: " << pivot_index << '\n';
			std::cout << "cycle support size: " << cycle.data().size() << '\n';
			return 0;
		}

		std::cout << "candidate from pivot " << pivot_index
			  << " did not verify; d(cycle) terms = " << d_cycle.size() << '\n';
		++pivot_index;
	}

	std::cout << "no non-zero cycle found inside the standardised insertion support\n";
	return 0;
}
