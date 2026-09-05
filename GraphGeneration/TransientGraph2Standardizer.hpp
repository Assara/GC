#pragma once

#include "GraphGeneration/TransientGraph2.hpp"
#include "GraphStandardizer.hpp"

namespace GraphGeneration {

template <
	Int N_VERTICES,
	Int N_EDGES
>
class transient_graph2_standardizer {
	public:
		using transient_type = transient_graph2<
			N_VERTICES,
			N_EDGES
		>;
		using graph_type = typename transient_type::graph_type;

		transient_type standardize_no_sign(const transient_type& input) const {
			const auto canonical = standardizer_.standardize_no_sign(input.graph());
			const auto valences = canonical.valence_array();
			// Sort the canonical labels, preserving their order within each valence.
			// This is deterministic for an isomorphism class and gives the pipeline
			// its required increasing-valence order without changing the old standardizer.
			std::array<Int, N_VERTICES> order{};
			for (Int v = 0; v < N_VERTICES; ++v) order[v] = v;
			std::stable_sort(order.begin(), order.end(), [&](Int a, Int b) {
				return valences[a] < valences[b];
			});
			Permutation<N_VERTICES> labels;
			for (Int v = 0; v < N_VERTICES; ++v) labels[order[v]] = v;
			graph_type result;
			result.assignPermutedDirectedSortedEdgesNoSign(canonical, labels);
			return transient_type(std::move(result));
		}

	private:
		GraphStandardizer<
			N_VERTICES,
			N_EDGES,
			0,
			0,
			0,
			0,
			fieldType
		> standardizer_;
};

} // namespace GraphGeneration
