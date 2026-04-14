#include <iostream>
#include <unordered_map>
#include <unordered_set>
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

std::unordered_set<GraphType> build_lie_support(const LieGC& lie) {
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

	std::size_t processed = 0;
	for (const auto& graph : support) {
		ContL d_graph = LieGC(graph, AssumeBasisOrderTag{}).d_contraction().data();
		d_graph.standardize_all();
		d_graph.sort_elements();
		differential_map.emplace(graph, std::move(d_graph));

		++processed;
		if (processed % 50 == 0) {
			std::cout << "built differential rows: " << processed
				  << " / " << support.size() << std::endl;
		}
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
	const auto w5 = build_w5();
	const auto w7 = build_w7();
	std::cout << "built W5 and W7 reps" << std::endl;

	const auto pre_5_7 = coalgebra_utils::gra_pre_lie_unrestricted(w5, w7);
	std::cout << "built unrestricted pre_5_7" << std::endl;
	const auto pre_7_5 = coalgebra_utils::gra_pre_lie_unrestricted(w7, w5);
	std::cout << "built unrestricted pre_7_5" << std::endl;
	const auto lie = coalgebra_utils::gra_lie_unrestricted(w5, w7);
	std::cout << "built unrestricted lie" << std::endl;
	const auto support = build_lie_support(lie);
	std::cout << "built unrestricted lie support" << std::endl;
	const auto differential_map = build_restricted_differential_map(support);
	std::cout << "built restricted differential map" << std::endl;

	std::vector<GraphType> ordered_support(support.begin(), support.end());
	std::sort(ordered_support.begin(), ordered_support.end());

	std::cout << "unrestricted pre_5_7 size: " << pre_5_7.data().size() << '\n';
	std::cout << "unrestricted pre_7_5 size: " << pre_7_5.data().size() << '\n';
	std::cout << "unrestricted lie size: " << lie.data().size() << '\n';
	std::cout << "unrestricted lie support size: " << ordered_support.size() << '\n';

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
			std::cout << "found non-zero cycle in unrestricted lie support" << std::endl;
			std::cout << "pivot index: " << pivot_index << '\n';
			std::cout << "cycle support size: " << cycle.data().size() << '\n';
			return 0;
		}

		std::cout << "candidate from pivot " << pivot_index
			  << " did not verify; d(cycle) terms = " << d_cycle.size() << std::endl;
		++pivot_index;
	}

	std::cout << "no non-zero cycle found inside the unrestricted lie support" << std::endl;
	return 0;
}
