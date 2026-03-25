#include <iostream>
#include <unordered_set>
#include <vector>

#include "examplegraphs.hpp"

namespace {

using FullGC = GC<13, 24, 0, 0, 0, 1>;
using ContGC = GC<12, 23, 0, 0, 0, 1>;
using FullGraph = typename FullGC::GraphType;
using ContGraph = typename ContGC::GraphType;
using FullLinComb = typename FullGC::L;

std::vector<std::vector<bool>> all_sequences(Int length) {
	std::vector<std::vector<bool>> result;
	const std::size_t count = std::size_t{1} << length;
	result.reserve(count);

	for (std::size_t mask = 0; mask < count; ++mask) {
		std::vector<bool> sequence;
		sequence.reserve(length);
		for (Int i = 0; i < length; ++i) {
			sequence.push_back(((mask >> i) & 1U) != 0U);
		}
		result.push_back(std::move(sequence));
	}

	return result;
}

FullGC build_sum(Int first_length, Int second_length) {
	FullGC result;
	const auto first_sequences = all_sequences(first_length);
	const auto second_sequences = all_sequences(second_length);

	for (const auto& first : first_sequences) {
		for (const auto& second : second_sequences) {
			std::vector<std::vector<bool>> blocks{first, second};
			result += FullGC(iterated_V_graph_lcr_extension<12>(blocks));
		}
	}

	result.standardize_all();
	result.sort_elements();
	return result;
}

template <typename GraphType>
bool is_three_four_valent(const GraphType& graph) {
	for (const Int valence : graph.valence_array()) {
		if (valence != 3 && valence != 4) {
			return false;
		}
	}
	return true;
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

FullGC build_three_four_split_candidates(const ContGC& remainder) {
	FullLinComb candidates;

	for (const auto& be : remainder.data()) {
		const ContGraph& graph = be.getValue();
		const auto valences = graph.valence_array();
		for (Int v = 0; v < ContGraph::N_VERTICES_; ++v) {
			if (valences[v] != 5 && valences[v] != 6) {
				continue;
			}

			typename ContGC::SplitL local_splits;
			graph.split_vertex_even(v, graph.adjacent(v), local_splits, be.getCoefficient());
			local_splits.standardize_all();

			for (const auto& split_be : local_splits) {
				if (split_be.getCoefficient() == fieldType{}) {
					continue;
				}
				if (!is_three_four_valent(split_be.getValue())) {
					continue;
				}
				candidates.append_in_basis_order(split_be);
			}
		}
	}

	candidates.standardize_and_sort();
	return FullGC(std::move(candidates.raw_elements_nonconst()));
}

FullGC filter_not_in_original(const FullGC& candidates, const FullGC& original) {
	std::unordered_set<FullGraph> original_support;
	original_support.reserve(original.data().size());
	for (const auto& be : original.data()) {
		original_support.insert(be.getValue());
	}

	FullLinComb filtered;
	for (const auto& be : candidates.data()) {
		if (!original_support.contains(be.getValue())) {
			filtered.append_in_basis_order(be);
		}
	}
	filtered.standardize_and_sort();
	return FullGC(std::move(filtered.raw_elements_nonconst()));
}

} // namespace

int main() {
	const FullGC sum_1_then_2 = build_sum(1, 2);
	const ContGC d_1_then_2 = [&]() {
		auto tmp = sum_1_then_2;
		return tmp.d_contraction();
	}();

	const FullGC split_candidates = build_three_four_split_candidates(d_1_then_2);
	const FullGC new_candidates = filter_not_in_original(split_candidates, sum_1_then_2);
	const ContGC d_new_candidates = [&]() {
		auto tmp = new_candidates;
		return tmp.d_contraction();
	}();

	auto residual_minus = d_1_then_2.data();
	residual_minus = residual_minus.add_scaled(d_new_candidates.data(), fieldType{-1});
	auto residual_plus = d_1_then_2.data();
	residual_plus = residual_plus.add_scaled(d_new_candidates.data(), fieldType{1});

	std::cout << "original sum terms: " << sum_1_then_2.data().size() << '\n';
	std::cout << "original d_contraction terms: " << d_1_then_2.data().size() << '\n';
	std::cout << "3/4-valent split candidates terms: " << split_candidates.data().size() << '\n';
	std::cout << "new candidates not already in original sum: " << new_candidates.data().size() << '\n';
	std::cout << "d(new candidates) terms: " << d_new_candidates.data().size() << '\n';
	std::cout << "common terms between original d and d(new candidates): "
		  << overlap_count(d_1_then_2.data(), d_new_candidates.data()) << '\n';
	std::cout << "size of d(original) - d(new candidates): " << residual_minus.size() << '\n';
	std::cout << "size of d(original) + d(new candidates): " << residual_plus.size() << '\n';

	return 0;
}
