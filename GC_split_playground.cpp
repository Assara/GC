#include <iostream>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <optional>

#include "examplegraphs.hpp"

using namespace std;

template <typename GCType>
void print_split_grading(const typename GCType::SplitGC& split_gc) {
	auto graded = split_gc.graded_by_odd_vertex_pairs();
	for (std::size_t i = 0; i < graded.size(); ++i) {
		if (graded[i].size() == 0) {
			continue;
		}
		std::cout << "odd-pair grade " << i << ": size=" << graded[i].size() << '\n';
	}
}

template <typename GCType>
typename GCType::SplitL filter_split_grades_above(
	const typename GCType::SplitL& split_l,
	signedInt cutoff_grade
) {
	typename GCType::SplitL filtered;
	for (const auto& be : split_l) {
		if (be.getValue().n_odd_pairs() > cutoff_grade) {
			filtered.append_in_basis_order(be);
		}
	}
	return filtered;
}

template <typename GCType>
typename GCType::L filter_graphs_on_exact_odd_pairs(
	const typename GCType::L& graph_l,
	signedInt target_grade
) {
	typename GCType::L filtered;
	for (const auto& be : graph_l) {
		if (be.getValue().n_odd_pairs() == target_grade) {
			filtered.append_in_basis_order(be);
		}
	}
	return filtered;
}

template <typename GCType>
std::vector<std::unordered_map<typename GCType::GraphType, typename GCType::L>>
make_staged_contracted_solver_maps(
	const std::vector<std::unordered_map<typename GCType::GraphType, typename GCType::SplitL>>& staged_delta_maps,
	signedInt filtered_out_grade,
	signedInt target_grade
) {
	using GraphType = typename GCType::GraphType;

	std::vector<std::unordered_map<GraphType, typename GCType::L>> solver_maps;
	solver_maps.reserve(staged_delta_maps.size());

	for (const auto& stage_delta_map : staged_delta_maps) {
		std::unordered_map<GraphType, typename GCType::L> solver_stage;
		solver_stage.reserve(stage_delta_map.size());

		for (const auto& [graph, delta_graph] : stage_delta_map) {
			auto filtered_delta = filter_split_grades_above<GCType>(delta_graph, filtered_out_grade);
			if (filtered_delta.size() == 0) {
				continue;
			}

			auto filtered_delta_copy = filtered_delta;
			typename GCType::SplitGC filtered_gc(filtered_delta_copy);
			auto contracted = filtered_gc.d_contraction();
			auto target_piece = filter_graphs_on_exact_odd_pairs<GCType>(contracted.data(), target_grade);
			if (target_piece.size() == 0) {
				continue;
			}

			solver_stage.emplace(graph, std::move(target_piece));
		}

		solver_maps.push_back(std::move(solver_stage));
	}

	return solver_maps;
}

template <typename GCType>
std::optional<GCType> try_reduce_split_obstructions_by_contraction(GCType current) {
	using GraphType = typename GCType::GraphType;
	using L = typename GCType::L;
	using SplitL = typename GCType::SplitL;

	std::vector<std::unordered_map<GraphType, SplitL>> staged_delta_maps;
	std::unordered_set<GraphType> known_graphs;

	while (true) {
		std::unordered_set<GraphType> current_support;
		current_support.reserve(current.size());
		for (const auto& be : current.data()) {
			current_support.insert(be.getValue());
		}

		auto split = current.delta();
		if (split.size() == 0) {
			return current;
		}

		const signedInt n = split.data().begin()->getValue().n_odd_pairs();
		for (const auto& be : split.data()) {
			if (be.getValue().n_odd_pairs() != n) {
				std::cout << "delta(G) is not concentrated in a single odd-pair grade: "
					  << be.getValue().n_odd_pairs() << " vs " << n << "\n";
				return std::nullopt;
			}
		}

		auto h = split.d_contraction();
		if (h.size() == 0) {
			std::cout << "d_contraction(delta(G)) is zero while delta(G) is non-zero\n";
			return std::nullopt;
		}

		auto h_down = filter_graphs_on_exact_odd_pairs<GCType>(h.data(), n - 1);
		if (h_down.size() == 0) {
			std::cout << "d_contraction(delta(G)) has no grade " << (n - 1) << " part\n";
			return std::nullopt;
		}

		std::unordered_map<GraphType, SplitL> stage_delta_map;
		for (const auto& be : h.data()) {
			const GraphType& graph = be.getValue();
			if (current_support.contains(graph)) {
				continue;
			}
			if (known_graphs.contains(graph)) {
				continue;
			}
			known_graphs.insert(graph);
			stage_delta_map.emplace(graph, GCType(graph, AssumeBasisOrderTag{}).delta().data());
		}
		staged_delta_maps.push_back(std::move(stage_delta_map));

		auto solver_maps = make_staged_contracted_solver_maps<GCType>(staged_delta_maps, n - 1, n - 1);
		auto solver = VectorSpace::wiedemann_primitive_finder<GraphType, GraphType, fieldType>(solver_maps);
		auto gamma_opt = solver.find_primitive_or_empty(h_down);
		if (!gamma_opt.has_value()) {
			std::cout << "failed to solve contracted split correction at grade " << (n - 1) << "\n";
			return std::nullopt;
		}

		GCType gamma(*gamma_opt);
		current += gamma.scalar_multiply(fieldType{-1});
	}
}

int run_split_playground_demo() {
	using WheelGC = OddGCdegZero<8>;

	WheelGC start(wheel_graph<7>());
	auto split = start.delta();

	std::cout << "start size = " << start.size() << '\n';
	std::cout << "delta size = " << split.size() << '\n';
	print_split_grading<WheelGC>(split);

	return 0;
}
