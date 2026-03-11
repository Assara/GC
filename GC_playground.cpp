#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <limits>
#include <optional>

#include "examplegraphs.hpp"

using namespace std;

constexpr bool kVerboseLocalSearch = false;

template <typename GCType>
std::vector<typename GCType::SplitGraphType> valence_4_splits(const typename GCType::GraphType& graph);

template <typename GCType>
void add_valence_4_split_boundaries(
	const GCType& cycle,
	std::unordered_map<typename GCType::SplitGraphType, typename GCType::L>& boundary_map
) {
	using SplitGraphType = typename GCType::SplitGraphType;
	using SplitGC = typename GCType::SplitGC;

	for (const auto& be : cycle.data()) {
		for (const SplitGraphType& split_graph : valence_4_splits<GCType>(be.getValue())) {
			if (boundary_map.contains(split_graph)) {
				continue;
			}
			boundary_map.emplace(
				split_graph,
				SplitGC(split_graph, AssumeBasisOrderTag{}).d_contraction().data()
			);
		}
	}
}

template <typename GCType>
std::optional<GCType> try_find_quadratic_contraction_representative_via_solver(GCType cycle) {
	using GraphType = typename GCType::GraphType;
	using SplitGraphType = typename GCType::SplitGraphType;
	using SplitGC = typename GCType::SplitGC;
	using L = typename GCType::L;

	std::unordered_map<SplitGraphType, L> boundary_map;
	SplitGC splits = cycle.delta();

	for (signedInt i = cycle.find_max_odd_pairs(); cycle.frontValence() > 4; --i) {
		for (const auto& split_be : splits.data()) {
			const SplitGraphType& split_graph = split_be.getValue();
			if (boundary_map.contains(split_graph)) {
				continue;
			}
			boundary_map.emplace(
				split_graph,
				SplitGC(split_graph, AssumeBasisOrderTag{}).d_contraction().data()
			);
		}

		auto solver =
			VectorSpace::wiedemann_primitive_finder<GraphType, SplitGraphType, fieldType>::create_filtered(
				boundary_map,
				[i](const GraphType& graph) { return graph.n_odd_pairs() >= i; }
			);

		auto primitive_opt = solver.find_primitive_or_empty(cycle.data());
		if (!primitive_opt.has_value()) {
			return std::nullopt;
		}

		SplitGC primitive(*primitive_opt);
		GCType full_correction = primitive.d_contraction();
		cycle += full_correction.scalar_multiply(fieldType{-1});
		splits = cycle.delta();
	}

	return cycle;
}

template <typename GCType>
std::optional<GCType> try_find_quadratic_contraction_representative_via_4valent_split_solver(GCType cycle) {
	using GraphType = typename GCType::GraphType;
	using SplitGraphType = typename GCType::SplitGraphType;
	using SplitGC = typename GCType::SplitGC;
	using L = typename GCType::L;

	std::unordered_map<SplitGraphType, L> boundary_map;

	for (signedInt i = cycle.find_max_odd_pairs(); cycle.frontValence() > 4; --i) {
		add_valence_4_split_boundaries(cycle, boundary_map);

		auto solver =
			VectorSpace::wiedemann_primitive_finder<GraphType, SplitGraphType, fieldType>::create_filtered(
				boundary_map,
				[i](const GraphType& graph) { return graph.n_odd_pairs() >= i; }
			);

		auto primitive_opt = solver.find_primitive_or_empty(cycle.data());
		if (!primitive_opt.has_value()) {
			return std::nullopt;
		}

		SplitGC primitive(*primitive_opt);
		GCType full_correction = primitive.d_contraction();
		cycle += full_correction.scalar_multiply(fieldType{-1});
	}

	return cycle;
}

template <typename GCType>
signedInt val3_weight_graph(const typename GCType::GraphType& graph) {
	signedInt count = 0;
	for (Int v : graph.valence_array()) {
		if (v == 3) {
			++count;
		}
	}
	return count - 4;
}

template <typename GCType>
signedInt odd_val_weight_graph(const typename GCType::GraphType& graph) {
	signedInt count = 0;
	for (Int v : graph.valence_array()) {
		if (v % 2 != 0) {
			++count;
		}
	}
	return count - 4;
}

template <typename GCType>
signedInt val3_weight_class(const GCType& gamma) {
	signedInt total = 0;
	for (const auto& be : gamma.data()) {
		total += val3_weight_graph<GCType>(be.getValue());
	}
	return total;
}

template <typename GCType>
signedInt odd_weight_class(const GCType& gamma) {
	signedInt total = 0;
	for (const auto& be : gamma.data()) {
		total += odd_val_weight_graph<GCType>(be.getValue());
	}
	return total;
}

struct ClassScore {
	signedInt val3_weight = 0;
	signedInt odd_weight = 0;
	bigInt support_size = 0;

	auto as_tuple() const {
		return std::tie(val3_weight, odd_weight, support_size);
	}

	bool operator<(const ClassScore& other) const {
		return as_tuple() < other.as_tuple();
	}
};

template <typename GCType>
ClassScore score_class(const GCType& gamma) {
	return ClassScore{
		val3_weight_class(gamma),
		odd_weight_class(gamma),
		gamma.data().size()
	};
}

template <typename GCType>
void normalize_class(GCType& gamma) {
	gamma.standardize_all();
	gamma.sort_elements();
}

template <typename GCType>
GCType single_graph_class(const typename GCType::GraphType& graph, fieldType coeff = fieldType{1}) {
	return GCType(BasisElement<typename GCType::GraphType, fieldType>(graph, coeff));
}

template <typename GCType>
struct ScoredClass {
	GCType gamma;
	ClassScore score;
};

template <typename GCType>
ScoredClass<GCType> make_scored_class(GCType gamma) {
	normalize_class(gamma);
	ClassScore score = score_class(gamma);
	return ScoredClass<GCType>{std::move(gamma), score};
}

template <typename GCType>
void push_unique_scored_class(std::vector<ScoredClass<GCType>>& classes, GCType gamma) {
	normalize_class(gamma);
	for (auto& existing : classes) {
		if (existing.gamma.data() == gamma.data()) {
			return;
		}
	}
	classes.push_back(make_scored_class<GCType>(std::move(gamma)));
}

template <typename GCType>
void sort_and_trim(
	std::vector<ScoredClass<GCType>>& classes,
	size_t max_size,
	size_t max_per_score = std::numeric_limits<size_t>::max()
) {
	std::sort(classes.begin(), classes.end(),
		[](const auto& a, const auto& b) { return a.score < b.score; });

	std::vector<ScoredClass<GCType>> trimmed;
	trimmed.reserve(std::min(max_size, classes.size()));

	for (auto& candidate : classes) {
		if (trimmed.size() >= max_size) {
			break;
		}

		size_t same_score_count = 0;
		bool same_support_for_score = false;
		for (const auto& kept : trimmed) {
			if (kept.score.val3_weight != candidate.score.val3_weight ||
				kept.score.odd_weight != candidate.score.odd_weight) {
				continue;
			}

			++same_score_count;
			if (kept.score.support_size == candidate.score.support_size) {
				same_support_for_score = true;
				break;
			}
		}

		if (same_support_for_score || same_score_count >= max_per_score) {
			continue;
		}

		trimmed.push_back(std::move(candidate));
	}

	classes = std::move(trimmed);
}

inline vector<vector<Int>> choose_k_indices(Int n, Int k) {
	vector<vector<Int>> result;
	vector<Int> current;

	auto rec = [&](auto&& self, Int start) -> void {
		if (static_cast<Int>(current.size()) == k) {
			result.push_back(current);
			return;
		}

		for (Int i = start; i <= n - (k - static_cast<Int>(current.size())); ++i) {
			current.push_back(i);
			self(self, i + 1);
			current.pop_back();
		}
	};

	if (k >= 0 && k <= n) {
		rec(rec, 0);
	}
	return result;
}

template <typename GCType>
std::vector<typename GCType::SplitGraphType> valence_4_splits(const typename GCType::GraphType& graph) {
	using SplitGraphType = typename GCType::SplitGraphType;

	std::vector<SplitGraphType> splits;
	const auto valences = graph.valence_array();

	for (Int v = 0; v < static_cast<Int>(valences.size()); ++v) {
		const auto adjacent = graph.adjacent(v);
		if (static_cast<Int>(adjacent.size()) <= 4) {
			continue;
		}

		for (const auto& moved_half_edges : choose_k_indices(static_cast<Int>(adjacent.size()), 3)) {
			splits.push_back(graph.splitGraph(v, adjacent, moved_half_edges));
		}
	}

	return splits;
}

template <typename GCType>
GCType trunk_contraction(const typename GCType::SplitGraphType& split_graph) {
	using GraphType = typename GCType::GraphType;
	using SplitGraphType = typename GCType::SplitGraphType;

	std::vector<BasisElement<GraphType, fieldType>> elems;
	elems.reserve(SplitGraphType::N_EDGES_ - 1);

	for (Int i = 0; i < SplitGraphType::N_EDGES_ - 1; ++i) {
		auto contracted = split_graph.contract_edge(i, fieldType{1});
		if (contracted.getCoefficient() != fieldType{0}) {
			elems.push_back(std::move(contracted));
		}
	}

	return GCType(std::move(elems));
}

template <typename GCType>
std::vector<GCType> split_con_candidates(const typename GCType::GraphType& graph) {
	std::vector<GCType> candidates;
	for (const auto& split_graph : valence_4_splits<GCType>(graph)) {
		candidates.push_back(trunk_contraction<GCType>(split_graph));
	}
	return candidates;
}

template <typename GCType>
const std::vector<GCType>& cached_split_con_candidates(const typename GCType::GraphType& graph) {
	using GraphType = typename GCType::GraphType;

	static std::unordered_map<GraphType, std::vector<GCType>> cache;

	auto it = cache.find(graph);
	if (it != cache.end()) {
		return it->second;
	}

	auto raw_candidates = split_con_candidates<GCType>(graph);
	std::vector<ScoredClass<GCType>> scored_candidates;
	scored_candidates.reserve(raw_candidates.size());
	for (auto& candidate : raw_candidates) {
		push_unique_scored_class(scored_candidates, std::move(candidate));
	}

	sort_and_trim(scored_candidates, scored_candidates.size(), scored_candidates.size());

	std::vector<GCType> sorted_candidates;
	sorted_candidates.reserve(scored_candidates.size());
	for (auto& scored : scored_candidates) {
		sorted_candidates.push_back(std::move(scored.gamma));
	}

	return cache.emplace(graph, std::move(sorted_candidates)).first->second;
}

template <typename GCType>
std::vector<GCType> best_split_con_candidates(const typename GCType::GraphType& graph, size_t branch_width) {
	const auto& sorted_candidates = cached_split_con_candidates<GCType>(graph);
	if (sorted_candidates.empty()) {
		return {single_graph_class<GCType>(graph)};
	}

	const size_t kept_count = std::min(branch_width, sorted_candidates.size());

	if (kVerboseLocalSearch) {
		std::cout << "cached candidates=" << sorted_candidates.size()
			  << " kept candidates=" << kept_count;
		if (!sorted_candidates.empty()) {
			const auto best = score_class(sorted_candidates.front());
			std::cout << " best=(" << best.val3_weight
				  << ", " << best.odd_weight
				  << ", " << best.support_size << ")";
		}
		std::cout << std::endl;
	}

	std::vector<GCType> best_candidates;
	best_candidates.reserve(kept_count);
	for (size_t i = 0; i < kept_count; ++i) {
		best_candidates.push_back(sorted_candidates[i]);
	}
	return best_candidates;
}

template <typename GCType>
std::vector<ScoredClass<GCType>> beam_split_con_list(const GCType& gamma, size_t branch_width, size_t beam_width) {
	std::vector<ScoredClass<GCType>> beam;
	beam.push_back(make_scored_class<GCType>(GCType{}));

	size_t term_index = 0;
	for (const auto& be : gamma.data()) {
		std::vector<GCType> term_candidates;
		if (val3_weight_graph<GCType>(be.getValue()) == 0) {
			term_candidates.push_back(single_graph_class<GCType>(be.getValue(), be.getCoefficient()));
		} else {
			term_candidates = best_split_con_candidates<GCType>(be.getValue(), branch_width);
			for (auto& candidate : term_candidates) {
				candidate.scalar_multiply(be.getCoefficient());
				normalize_class(candidate);
			}
		}

		std::vector<ScoredClass<GCType>> next_beam;
		next_beam.reserve(beam.size() * std::max<size_t>(size_t{1}, term_candidates.size()));

		for (const auto& partial : beam) {
			for (const auto& term_candidate : term_candidates) {
				GCType combined = partial.gamma;
				combined += term_candidate;
				push_unique_scored_class(next_beam, std::move(combined));
			}
		}

		sort_and_trim(next_beam, beam_width, 4);
		beam = std::move(next_beam);

		if (kVerboseLocalSearch && !beam.empty()) {
			const auto& best = beam.front().score;
			std::cout << "term " << term_index
				  << " beam_size=" << beam.size()
				  << " best=(" << best.val3_weight
				  << ", " << best.odd_weight
				  << ", " << best.support_size << ")"
				  << std::endl;
		}
		++term_index;
	}

	return beam;
}

template <typename GCType>
std::vector<GCType> make_list(
	const GCType& start,
	size_t length = 10,
	size_t branch_width = 8,
	size_t class_beam_width = 32,
	size_t level_beam_width = 32
) {
	std::vector<GCType> result;
	std::vector<ScoredClass<GCType>> current_beam;
	current_beam.push_back(make_scored_class<GCType>(start));
	result.push_back(current_beam.front().gamma);

	for (size_t i = 0; i < length; ++i) {
		std::vector<ScoredClass<GCType>> next_beam;

		for (const auto& state : current_beam) {
			auto local_beam = beam_split_con_list(state.gamma, branch_width, class_beam_width);
			for (auto& local_state : local_beam) {
				push_unique_scored_class(next_beam, std::move(local_state.gamma));
			}
		}

		sort_and_trim(next_beam, level_beam_width, 8);

		if (next_beam.empty()) {
			break;
		}

		current_beam = std::move(next_beam);
		const auto& best = current_beam.front().score;
		std::cout << "depth " << i + 1
			  << " global_beam_size=" << current_beam.size()
			  << " best=(" << best.val3_weight
			  << ", " << best.odd_weight
			  << ", " << best.support_size << ")"
			  << std::endl;

		result.push_back(current_beam.front().gamma);
		if (current_beam.front().score.val3_weight == 0) {
			break;
		}
	}

	return result;
}

int main() {
	using WheelGC = OddGCdegZero<16>;

	WheelGC start(wheel_graph<15>());
	auto levels = make_list(start, 14, 8, 16, 32);

	for (size_t i = 0; i < levels.size(); ++i) {
		auto level = levels[i];
		level.standardize_all();
		std::cout << "level " << i
			  << ": size=" << level.size()
			  << " val3_weight=" << val3_weight_class(level)
			  << " odd_weight=" << odd_weight_class(level)
			  << std::endl;
	}

	auto final_level = levels.back();
	final_level.standardize_all();
	std::cout << "final level:" << std::endl;
	final_level.print();

	auto d_final = final_level.d_contraction();
	std::cout << "d_contraction size = " << d_final.size() << std::endl;
	if (d_final.size() != 0) {
		d_final.print();
	}

	return 0;
}
