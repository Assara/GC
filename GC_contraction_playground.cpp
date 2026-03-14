#include <iostream>
#include <vector>
#include <algorithm>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
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
std::unordered_map<typename GCType::SplitGraphType, typename GCType::L> make_min_triangle_drop_valence_4_split_boundary_stage(
	const GCType& cycle,
	const std::unordered_set<typename GCType::SplitGraphType>& known_split_graphs
) {
	using GraphType = typename GCType::GraphType;
	using SplitGraphType = typename GCType::SplitGraphType;
	using SplitGC = typename GCType::SplitGC;

	std::unordered_map<SplitGraphType, typename GCType::L> stage_map;

	for (const auto& be : cycle.data()) {
		const GraphType& graph = be.getValue();
		const bigInt triangle_count = graph.count_triangles();
		for (const SplitGraphType& split_graph : valence_4_splits<GCType>(graph)) {
			if (triangle_count - split_graph.count_triangles() > 2) {
				continue;
			}
			if (known_split_graphs.contains(split_graph)) {
				continue;
			}
			stage_map.emplace(
				split_graph,
				SplitGC(split_graph, AssumeBasisOrderTag{}).d_contraction().data()
			);
		}
	}

	return stage_map;
}

template <typename GCType>
void add_min_triangle_drop_valence_4_split_boundaries(
	const GCType& cycle,
	std::unordered_map<typename GCType::SplitGraphType, typename GCType::L>& boundary_map
) {
	using GraphType = typename GCType::GraphType;
	using SplitGraphType = typename GCType::SplitGraphType;
	using SplitGC = typename GCType::SplitGC;

	for (const auto& be : cycle.data()) {
		const GraphType& graph = be.getValue();
		const bigInt triangle_count = graph.count_triangles();
		for (const SplitGraphType& split_graph : valence_4_splits<GCType>(graph)) {
			if (triangle_count - split_graph.count_triangles() > 2) {
				continue;
			}
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

template <typename GraphType, typename SplitGraphType, typename L, typename Predicate>
std::optional<VectorSpace::LinComb<SplitGraphType, fieldType>> find_staged_primitive_or_empty(
	const std::vector<std::unordered_map<SplitGraphType, L>>& staged_boundary_maps,
	const VectorSpace::LinComb<GraphType, fieldType>& target,
	const Predicate& predicate
) {
	using Solver = VectorSpace::wiedemann_primitive_finder<GraphType, SplitGraphType, fieldType>;

	if (staged_boundary_maps.empty()) {
		return std::nullopt;
	}

	auto full_solver = Solver::create_filtered(staged_boundary_maps, predicate);
	return full_solver.find_primitive_or_empty(target);
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

	signedInt stage = 0;
	for (signedInt i = cycle.find_max_odd_pairs(); cycle.frontValence() > 4; --i) {
		stage++;
		std::cout << "stage: " << stage << std::endl;
		const signedInt current_4val = stage*2;
		add_valence_4_split_boundaries(cycle, boundary_map);

		auto solver =
			VectorSpace::wiedemann_primitive_finder<GraphType, SplitGraphType, fieldType>::create_filtered(
				boundary_map,
				[i, current_4val](const GraphType& graph) {
					return graph.n_odd_pairs() >= i &&
						graph.n_4valent_vertices() < current_4val;
				}
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
std::optional<GCType> try_find_quadratic_contraction_representative_via_min_triangle_split_solver(GCType cycle) {
	using GraphType = typename GCType::GraphType;
	using SplitGraphType = typename GCType::SplitGraphType;
	using SplitGC = typename GCType::SplitGC;
	using L = typename GCType::L;

	std::vector<std::unordered_map<SplitGraphType, L>> staged_boundary_maps;
	std::unordered_set<SplitGraphType> known_split_graphs;

	for (signedInt i = cycle.find_max_odd_pairs(); cycle.frontValence() > 4; --i) {
		const signedInt current_4val = min_4valent_vertices_class(cycle);
		const signedInt max_4val = max_4valent_vertices_class(cycle);
		if (current_4val != max_4val) {
			std::cout << "warning: mixed 4-valent counts in class: min="
				  << current_4val << " max=" << max_4val << std::endl;
		}

		auto stage_map = make_min_triangle_drop_valence_4_split_boundary_stage<GCType>(cycle, known_split_graphs);
		for (const auto& [split_graph, _] : stage_map) {
			known_split_graphs.insert(split_graph);
		}
		staged_boundary_maps.push_back(std::move(stage_map));
		if (staged_boundary_maps.back().empty() && staged_boundary_maps.size() == 1) {
			return std::nullopt;
		}

		auto primitive_opt = find_staged_primitive_or_empty<GraphType, SplitGraphType, L>(
			staged_boundary_maps,
			cycle.data(),
			[i, current_4val](const GraphType& graph) {
				return graph.n_odd_pairs() >= i &&
					graph.n_4valent_vertices() <= current_4val;
			}
		);
		if (!primitive_opt.has_value()) {
			return std::nullopt;
		}

		SplitGC primitive(*primitive_opt);
		GCType full_correction = primitive.d_contraction();
		if (full_correction.size() == 0) {
			return std::nullopt;
		}
		cycle += full_correction.scalar_multiply(fieldType{-1});
	}

	return cycle;
}

template <typename GCType>
std::optional<GCType> try_find_odd_preserving_split_lift(const GCType& seed, signedInt n) {
	using GraphType = typename GCType::GraphType;
	using SplitGraphType = typename GCType::SplitGraphType;
	using SplitGC = typename GCType::SplitGC;
	using SplitL = typename GCType::SplitL;

	const signedInt target_grade = n + 1;
	SplitGC target = seed.delta();

	if (target.size() == 0) {
		std::cout << "delta(G) is zero\n";
		return GCType{};
	}

	for (const auto& be : target.data()) {
		if (be.getValue().n_odd_pairs() != target_grade) {
			std::cout << "delta(G) contains grade " << be.getValue().n_odd_pairs()
				  << ", expected only grade " << target_grade << "\n";
			return std::nullopt;
		}
	}

	std::unordered_map<GraphType, SplitL> boundary_map;
	for (const auto& split_be : target.data()) {
		const SplitGraphType& split_graph = split_be.getValue();
		for (Int i = 0; i < SplitGraphType::N_EDGES_; ++i) {
			auto contracted = split_graph.contract_edge(i, fieldType{1});
			if (contracted.getCoefficient() == fieldType{0}) {
				continue;
			}
			if (contracted.getValue().n_odd_pairs() != target_grade) {
				continue;
			}

			GraphType::std(contracted);
			if (contracted.getCoefficient() == fieldType{0}) {
				continue;
			}

			const GraphType& graph = contracted.getValue();
			if (boundary_map.contains(graph)) {
				continue;
			}

			boundary_map.emplace(
				graph,
				GCType(graph, AssumeBasisOrderTag{}).delta().data()
			);
		}
	}

	std::cout << "odd-preserving contraction domain size = " << boundary_map.size() << "\n";

	auto solver =
		VectorSpace::wiedemann_primitive_finder<SplitGraphType, GraphType, fieldType>::create_filtered(
			boundary_map,
			[target_grade](const SplitGraphType& split_graph) {
				return split_graph.n_odd_pairs() <= target_grade;
			}
		);

	auto primitive_opt = solver.find_primitive_or_empty(target.data());
	if (!primitive_opt.has_value()) {
		return std::nullopt;
	}

	GCType primitive(*primitive_opt);
	primitive.standardize_all();
	return primitive;
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

template <typename GCType>
signedInt min_4valent_vertices_class(const GCType& gamma) {
	if (gamma.data().raw_elements().empty()) {
		return 0;
	}

	signedInt result = std::numeric_limits<signedInt>::max();
	for (const auto& be : gamma.data()) {
		result = std::min(result, be.getValue().n_4valent_vertices());
	}
	return result;
}

template <typename GCType>
signedInt max_4valent_vertices_class(const GCType& gamma) {
	signedInt result = 0;
	for (const auto& be : gamma.data()) {
		result = std::max(result, be.getValue().n_4valent_vertices());
	}
	return result;
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
	using WheelGC = OddGCdegZero<18>;

	WheelGC start(wheel_graph<17>());
	auto rep = try_find_quadratic_contraction_representative_via_min_triangle_split_solver(start);
	if (!rep.has_value()) {
		std::cout << "no solution" << std::endl;
		return 1;
	}

	auto gamma = *rep;
	gamma.standardize_all();
	gamma.sort_elements();
	std::cout << "representative size = " << gamma.size() << std::endl;
	std::cout << "val3_weight = " << val3_weight_class(gamma) << std::endl;
	std::cout << "odd_weight = " << odd_weight_class(gamma) << std::endl;

	auto d_final = gamma.d_contraction();
	std::cout << "d_contraction size = " << d_final.size() << std::endl;
	if (d_final.size() != 0) {
		d_final.print();
	}

	return 0;
}
