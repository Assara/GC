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

std::unordered_map<GraphType, ContL> build_differential_map(
	const std::unordered_set<GraphType>& support
) {
	std::unordered_map<GraphType, ContL> differential_map;
	differential_map.reserve(support.size());
	std::size_t processed = 0;

	for (const auto& graph : support) {
		ContL differential = LieGC(graph, AssumeBasisOrderTag{}).d_contraction().data();
		differential.standardize_all();
		differential.sort_elements();
		differential_map.emplace(graph, std::move(differential));
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

std::optional<LieGC> find_cycle_with_pivots(
	const std::vector<GraphType>& pivots,
	const std::unordered_map<GraphType, ContL>& differential_map
) {
	std::size_t pivot_index = 0;
	for (const auto& pivot_graph : pivots) {
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

		std::cout << "trying pivot " << pivot_index << " / " << pivots.size() << std::endl;
		VectorSpace::wiedemann_primitive_finder<ContGraphType, GraphType, fieldType> solver(restricted_map);
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
			return cycle;
		}

		std::cout << "candidate from pivot " << pivot_index
			  << " did not verify; d(cycle) terms = " << d_cycle.size() << std::endl;
		++pivot_index;
	}

	return std::nullopt;
}

} // namespace

int main() {
	LieGC lie = build_full_lie_bracket();

	std::unordered_set<GraphType> original_support;
	original_support.reserve(lie.data().size());
	std::vector<GraphType> original_pivots;
	original_pivots.reserve(lie.data().size());
	for (const auto& be : lie.data()) {
		original_support.insert(be.getValue());
		original_pivots.push_back(be.getValue());
	}
	std::cout << "built lie bracket support with " << original_support.size() << " terms" << std::endl;

	auto d_lie = lie.d_contraction();
	d_lie.standardize_all();
	d_lie.sort_elements();
	std::cout << "built and standardised d(lie) with " << d_lie.data().size() << " terms" << std::endl;

	const auto extra_support = build_zero_cobracket_split_candidates(d_lie);
	std::cout << "finished zero-cobracket split harvest" << std::endl;

	std::unordered_set<GraphType> augmented_support = original_support;
	augmented_support.insert(extra_support.begin(), extra_support.end());

	const auto differential_map = build_differential_map(augmented_support);
	std::cout << "finished augmented differential map" << std::endl;

	std::cout << "original lie support size: " << original_support.size() << '\n';
	std::cout << "standardised d(lie) terms: " << d_lie.data().size() << '\n';
	std::cout << "extra zero-cobracket split candidates: " << extra_support.size() << '\n';
	std::cout << "augmented support size: " << augmented_support.size() << '\n';

	auto cycle_opt = find_cycle_with_pivots(original_pivots, differential_map);
	if (!cycle_opt.has_value()) {
		std::cout << "no non-zero cycle found using the augmented map with an original lie-support pivot\n";
		return 0;
	}

	auto cycle = *cycle_opt;
	std::size_t overlap_with_original = 0;
	for (const auto& be : cycle.data()) {
		if (original_support.contains(be.getValue())) {
			++overlap_with_original;
		}
	}

	std::cout << "found non-zero cycle in augmented support\n";
	std::cout << "cycle support size: " << cycle.data().size() << '\n';
	std::cout << "cycle terms from original lie support: " << overlap_with_original << '\n';
	std::cout << "cycle terms from added zero-cobracket support: "
		  << (cycle.data().size() - overlap_with_original) << '\n';

	return 0;
}
