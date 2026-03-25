#include <iostream>
#include <unordered_map>
#include <vector>

#include "coalgebra_utils.hpp"
#include "examplegraphs.hpp"

namespace {

using W5GC = OddGCdegZero<6>;
using W7GC = OddGCdegZero<8>;
using LieGC = coalgebra_utils::GraInsertLowValentGC<W5GC, W7GC>;
using GraphType = LieGC::GraphType;
using ContL = LieGC::ContL;
using L = LieGC::L;

LieGC build_full_lie_bracket() {
	std::vector<bool> Lseq{false};
	std::vector<bool> LL{false, false};
	std::vector<bool> LR{false, true};

	auto w5 = W5GC(V_graph<5>(Lseq));
	auto w7 = W7GC(V_graph<7>(LL));
	w7 += W7GC(V_graph<7>(LR));

	auto lie = coalgebra_utils::gra_lie(w5, w7);
	lie.standardize_all();
	lie.sort_elements();
	return lie;
}

std::unordered_map<GraphType, ContL> build_restricted_differential_map(const LieGC& lie) {
	std::unordered_map<GraphType, ContL> differential_map;
	differential_map.reserve(lie.data().size());

	for (const auto& be : lie.data()) {
		ContL d_graph = LieGC(be.getValue(), AssumeBasisOrderTag{}).d_contraction().data();
		d_graph.standardize_all();
		d_graph.sort_elements();
		differential_map.emplace(be.getValue(), std::move(d_graph));
	}

	return differential_map;
}

LieGC make_cycle_from_solution(const GraphType& pivot_graph, L solution) {
	solution.append_in_basis_order(pivot_graph, fieldType{1});
	solution.standardize_and_sort();
	return LieGC(solution);
}

} // namespace

int main() {
	const LieGC lie = build_full_lie_bracket();
	const auto differential_map = build_restricted_differential_map(lie);

	std::cout << "full lie support size: " << lie.data().size() << '\n';

	std::size_t pivot_index = 0;
	for (const auto& pivot_be : lie.data()) {
		const GraphType& pivot_graph = pivot_be.getValue();

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

		std::cout << "trying pivot " << pivot_index << " / " << lie.data().size() << '\n';
		VectorSpace::wiedemann_primitive_finder<LieGC::ContGraphType, GraphType, fieldType> solver(restricted_map);
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
			std::cout << "found non-zero cycle in lie support\n";
			std::cout << "pivot index: " << pivot_index << '\n';
			std::cout << "cycle support size: " << cycle.data().size() << '\n';
			return 0;
		}

		std::cout << "candidate from pivot " << pivot_index
			  << " did not verify; d(cycle) terms = " << d_cycle.size() << '\n';
		++pivot_index;
	}

	std::cout << "no non-zero cycle found inside the full lie support\n";
	return 0;
}
