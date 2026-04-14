#include <iostream>
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
using ContGraphType = SearchGC::ContGraphType;

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

bool is_three_four_valent(const GraphType& graph) {
	for (const Int valence : graph.valence_array()) {
		if (valence != 3 && valence != 4) {
			return false;
		}
	}
	return true;
}

std::unordered_set<GraphType> build_restricted_lie_support(const W5GC& w5, const W7GC& w7) {
	const auto lie = coalgebra_utils::gra_lie(w5, w7);

	std::unordered_set<GraphType> support;
	support.reserve(lie.data().size());
	for (const auto& be : lie.data()) {
		support.insert(be.getValue());
	}
	return support;
}

std::unordered_set<GraphType> build_three_four_split_support(
	const std::unordered_set<GraphType>& restricted_support
) {
	std::unordered_set<GraphType> split_support;

	std::size_t processed_domain = 0;
	for (const auto& graph : restricted_support) {
		auto differential = SearchGC(graph, AssumeBasisOrderTag{}).d_contraction();
		differential.standardize_all();
		differential.sort_elements();

		for (const auto& be : differential.data()) {
			const ContGraphType& h = be.getValue();

			for (Int vertex = 0; vertex < ContGraphType::N_VERTICES_; ++vertex) {
				typename SearchGC::ContGC::SplitL local_splits;
				h.split_vertex_even(vertex, h.adjacent(vertex), local_splits, fieldType{1});
				local_splits.standardize_all();
				local_splits.sort_elements();

				for (const auto& split_be : local_splits) {
					if (split_be.getCoefficient() == fieldType{0}) {
						continue;
					}
					if (is_three_four_valent(split_be.getValue())) {
						split_support.insert(split_be.getValue());
					}
				}
			}
		}

		++processed_domain;
		if (processed_domain % 25 == 0) {
			std::cout << "processed restricted support graphs: " << processed_domain
				  << " / " << restricted_support.size()
				  << ", split support so far: " << split_support.size() << std::endl;
		}
	}

	return split_support;
}

} // namespace

int main() {
	const auto w5 = build_w5();
	const auto w7 = build_w7();

	const auto restricted_support = build_restricted_lie_support(w5, w7);
	std::cout << "restricted lie support size: " << restricted_support.size() << std::endl;

	const auto split_support = build_three_four_split_support(restricted_support);
	std::cout << "three-four split support size: " << split_support.size() << std::endl;

	std::unordered_set<GraphType> augmented_support = restricted_support;
	augmented_support.insert(split_support.begin(), split_support.end());
	std::cout << "augmented search space size: " << augmented_support.size() << std::endl;

	auto cycle_opt = coalgebra_search::find_cobracket_lift_in_support<SearchGC>(
		augmented_support,
		w5,
		w7,
		true
	);

	if (!cycle_opt.has_value()) {
		std::cout << "no cycle with target cobracket W5 ^ W7 found in restricted-split augmented space" << std::endl;
		return 0;
	}

	const auto& cycle = *cycle_opt;
	std::cout << "found cycle with target cobracket W5 ^ W7 in restricted-split augmented space" << std::endl;
	std::cout << "cycle support size: " << cycle.data().size() << std::endl;

	return 0;
}
