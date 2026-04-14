#include <iostream>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "DegreeZeroGraphWedge.hpp"
#include "coalgebra_utils.hpp"
#include "examplegraphs.hpp"

namespace {

using W5GC = OddGCdegZero<6>;
using W7GC = OddGCdegZero<8>;
using LieGC = coalgebra_utils::GraInsertLowValentGC<W5GC, W7GC>;
using GraphType = LieGC::GraphType;
using ContGraphType = LieGC::ContGraphType;
using ContGC = LieGC::ContGC;
using ContL = LieGC::ContL;
using L = LieGC::L;
using OneWedge = DegreeZeroGraphWedge<48, 1>;

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

LieGC build_full_lie_bracket() {
	const auto w5 = build_w5();
	const auto w7 = build_w7();

	auto lie = coalgebra_utils::gra_lie(w5, w7);
	lie.standardize_all();
	lie.sort_elements();
	return lie;
}

GraphType build_lcr_pivot() {
	std::vector<bool> first{false};
	std::vector<bool> second{false, false};
	std::vector<std::vector<bool>> blocks{first, second};
	auto pivot = iterated_V_graph_lcr_extension<12>(blocks);
	BasisElement<GraphType, fieldType> be(pivot, fieldType{1});
	return GraphType::canonized(be).getValue();
}

bool cobracket_vanishes(const GraphType& graph) {
	auto wedge = OneWedge::from_graph(graph);
	auto cobracket = wedge.cobracket();
	cobracket.standardize_and_sort();
	return cobracket.size() == 0;
}

std::unordered_set<GraphType> build_zero_cobracket_split_candidates(const ContGC& differential) {
	std::unordered_set<GraphType> support;
	std::size_t processed_remainder_terms = 0;
	std::size_t processed_split_terms = 0;

	for (const auto& be : differential.data()) {
		const ContGraphType& graph = be.getValue();
		const auto valences = graph.valence_array();
		++processed_remainder_terms;
		if (processed_remainder_terms % 50 == 0) {
			std::cout << "processed remainder terms: " << processed_remainder_terms
				  << ", zero-cobracket support so far: " << support.size() << std::endl;
		}

		for (Int vertex = 0; vertex < ContGraphType::N_VERTICES_; ++vertex) {
			if (valences[vertex] != 5 && valences[vertex] != 6) {
				continue;
			}

			typename ContGC::SplitL local_splits;
			graph.split_vertex_even(vertex, graph.adjacent(vertex), local_splits, fieldType{1});
			local_splits.standardize_all();
			local_splits.sort_elements();

			for (const auto& split_be : local_splits) {
				if (split_be.getCoefficient() == fieldType{}) {
					continue;
				}
				++processed_split_terms;
				if (processed_split_terms % 500 == 0) {
					std::cout << "processed split terms: " << processed_split_terms
						  << ", zero-cobracket support so far: " << support.size() << std::endl;
				}

				if (cobracket_vanishes(split_be.getValue())) {
					support.insert(split_be.getValue());
				}
			}
		}
	}

	return support;
}

std::unordered_map<GraphType, ContL> build_differential_map(const std::unordered_set<GraphType>& support) {
	std::unordered_map<GraphType, ContL> differential_map;
	differential_map.reserve(support.size());
	std::size_t processed = 0;

	for (const auto& graph : support) {
		ContL differential = LieGC(graph, AssumeBasisOrderTag{}).d_contraction().data();
		differential.standardize_all();
		differential.sort_elements();
		differential_map.emplace(graph, std::move(differential));
		++processed;
		if (processed % 100 == 0) {
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

std::optional<LieGC> find_cycle_for_pivot(
	const GraphType& pivot_graph,
	const std::unordered_map<GraphType, ContL>& differential_map
) {
	LieGC pivot_gc(pivot_graph, AssumeBasisOrderTag{});
	auto target_gc = pivot_gc.d_contraction();
	target_gc.standardize_all();
	target_gc.sort_elements();
	ContL target = target_gc.data();
	target.scalar_multiply(fieldType{-1});

	std::unordered_map<GraphType, ContL> restricted_map;
	restricted_map.reserve(differential_map.size());
	for (const auto& [graph, differential] : differential_map) {
		if (graph == pivot_graph) {
			continue;
		}
		restricted_map.emplace(graph, differential);
	}

	std::cout << "starting solver for fixed LCR pivot" << std::endl;
	VectorSpace::wiedemann_primitive_finder<ContGraphType, GraphType, fieldType> solver(restricted_map);
	auto solution_opt = solver.find_primitive_or_empty(target);
	if (!solution_opt.has_value()) {
		return std::nullopt;
	}

	auto cycle = make_cycle_from_solution(pivot_graph, *solution_opt);
	auto d_cycle = cycle.d_contraction();
	d_cycle.standardize_all();
	d_cycle.sort_elements();

	if (d_cycle.size() == 0) {
		return cycle;
	}

	std::cout << "candidate did not verify; d(cycle) terms = " << d_cycle.size() << std::endl;
	return std::nullopt;
}

} // namespace

int main() {
	LieGC lie = build_full_lie_bracket();
	std::cout << "built lie bracket" << std::endl;

	auto d_lie = lie.d_contraction();
	d_lie.standardize_all();
	d_lie.sort_elements();
	std::cout << "built and standardised d(lie) with " << d_lie.data().size() << " terms" << std::endl;

	const auto extra_support = build_zero_cobracket_split_candidates(d_lie);
	std::cout << "finished zero-cobracket split harvest" << std::endl;

	std::unordered_set<GraphType> support;
	support.reserve(lie.data().size() + extra_support.size() + 1);
	for (const auto& be : lie.data()) {
		support.insert(be.getValue());
	}
	support.insert(extra_support.begin(), extra_support.end());

	const GraphType pivot_graph = build_lcr_pivot();
	if (!support.contains(pivot_graph)) {
		std::cerr << "fixed LCR pivot is not in the augmented support" << std::endl;
		return 1;
	}

	const auto differential_map = build_differential_map(support);
	std::cout << "finished augmented differential map" << std::endl;
	std::cout << "lie support size: " << lie.data().size() << std::endl;
	std::cout << "extra zero-cobracket split candidates: " << extra_support.size() << std::endl;
	std::cout << "augmented support size: " << support.size() << std::endl;

	auto cycle_opt = find_cycle_for_pivot(pivot_graph, differential_map);
	if (!cycle_opt.has_value()) {
		std::cout << "no cycle found with the fixed LCR pivot" << std::endl;
		return 0;
	}

	const auto& cycle = *cycle_opt;
	std::size_t lie_terms = 0;
	for (const auto& be : cycle.data()) {
		if (std::ranges::find_if(
			lie.data().raw_elements(),
			[&](const auto& lie_be) { return lie_be.getValue() == be.getValue(); }
		) != lie.data().raw_elements().end()) {
			++lie_terms;
		}
	}

	std::cout << "found cycle with fixed LCR pivot" << std::endl;
	std::cout << "cycle support size: " << cycle.data().size() << std::endl;
	std::cout << "cycle terms coming from the original lie support: " << lie_terms << std::endl;
	std::cout << "cycle terms outside the original lie support: " << (cycle.data().size() - lie_terms) << std::endl;

	return 0;
}
