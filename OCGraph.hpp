#pragma once

#include <algorithm>
#include <array>
#include <map>
#include <optional>
#include <vector>

#include "GraphDirections.hpp"
#include "GraphIsomorphism.hpp"
#include "GraphStandardizer.hpp"
#include "VectorSpace/LinComb.hpp"

template <typename GraphType>
class OCGraph : public VectorSpace::LinComb<GraphDirections<GraphType>, fieldType> {
	public:
		using DirectionType = GraphDirections<GraphType>;
		using DirectionComb = VectorSpace::LinComb<DirectionType, fieldType>;
		using Basis = BasisElement<OCGraph, fieldType>;
		using Element = BasisElement<DirectionType, fieldType>;
		using SplitGraphType = typename GraphType::SplitGraph;
		using SplitOCGraphType = OCGraph<SplitGraphType>;
		using SplitLinComb = VectorSpace::LinComb<SplitOCGraphType, fieldType>;
		using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;
		using Standardizer = GraphStandardizer<
			GraphType::N_VERTICES_,
			GraphType::N_EDGES_,
			GraphType::N_OUT_HAIR_,
			GraphType::N_IN_HAIR_,
			GraphType::C_,
			GraphType::D_,
			fieldType
		>;

		OCGraph() = default;

		explicit OCGraph(const GraphType& input_graph)
			: graph(input_graph) {}

		explicit OCGraph(GraphType input_graph, DirectionComb input_directions, std::vector<IsoType> input_automorphisms = {})
			: DirectionComb(std::move(input_directions)),
			  graph(std::move(input_graph)),
			  automorphisms(std::move(input_automorphisms)) {}

		signedInt compare(const OCGraph& other) const {
			return graph.compare(other.graph);
		}

		bool operator<(const OCGraph& other) const {
			return compare(other) < 0;
		}

		bool operator==(const OCGraph& other) const {
			return compare(other) == 0;
		}

		static Basis merge(Basis a, Basis b) {
			DirectionComb merged = static_cast<const DirectionComb&>(a.getValue()) * a.getCoefficient();
			merged += static_cast<const DirectionComb&>(b.getValue()) * b.getCoefficient();

			OCGraph merged_graph(
				a.getValue().graph,
				std::move(merged),
				a.getValue().automorphisms.empty() ? b.getValue().automorphisms : a.getValue().automorphisms
			);

			const fieldType outer_coeff = merged_graph.size() == 0 ? fieldType{0} : fieldType{1};
			return Basis(std::move(merged_graph), outer_coeff);
		}

		OCGraph delta_e() const {
			OCGraph result(graph, DirectionComb{}, automorphisms);

			for (const Element& elem : this->raw_elements()) {
				std::array<Int, GraphType::N_VERTICES_> true_counts{};
				for (Int i = 0; i < GraphType::SIZE; ++i) {
					if (elem.getValue()[i]) {
						++true_counts[graph.half_edges[i]];
					}
				}

				fieldType sign = fieldType{1};
				for (Int offset = 0; offset < GraphType::SIZE; ++offset) {
					const Int i = GraphType::SIZE - 1 - offset;
					if (!elem.getValue()[i]) {
						continue;
					}
					if (true_counts[graph.half_edges[i]] <= 1) {
						sign = -sign;
						continue;
					}

					DirectionType next = elem.getValue();
					next[i] = false;
					result.append_in_basis_order(next, elem.getCoefficient() * sign);
					sign = -sign;
				}
			}

			if (result.automorphisms.empty()) {
				result.sort_elements();
			} else {
				result.standardize_directions_without_vertex_check();
				result.sort_elements();
			}
			return result;
		}

		SplitLinComb delta_v() const {
			SplitLinComb result;
			const auto graph_splits = graph.unsorted_full_splits(fieldType{1});
			result.reserve(graph_splits.size());

			for (const auto& split_be : graph_splits.raw_elements()) {
				typename SplitOCGraphType::DirectionComb lifted_directions;

				for (const Element& elem : this->raw_elements()) {
					typename SplitOCGraphType::DirectionType lifted;
					lifted.fill(false);
					for (Int i = 0; i < GraphType::SIZE; ++i) {
						lifted[i] = elem.getValue()[i];
					}
					lifted[SplitGraphType::SIZE - 1] = true;
					lifted_directions.append_in_basis_order(lifted, elem.getCoefficient());
				}

				if (lifted_directions.size() == 0) {
					continue;
				}

				SplitOCGraphType split_ocgraph(split_be.getValue(), std::move(lifted_directions));
				result.append_in_basis_order(split_ocgraph, split_be.getCoefficient());
			}

			return result;
		}

		static OCGraph OCG_F0(const GraphType& input_graph) {
			OCGraph result(input_graph, DirectionComb{});
			result.standardize_graph_only();

			const Element minimal = result.minimal_one_true_direction_element();
			if (minimal.getCoefficient() != fieldType{0}) {
				result.append_in_basis_order(minimal);
			}

			result.sort_elements();
			return result;
		}

		std::optional<OCGraph> greedy_delta_e_primitive_step() const {
			OCGraph source = normalized_copy();
			if (source.size() == 0) {
				return OCGraph(graph, DirectionComb{}, automorphisms);
			}

			const DirectionType current_max = source.raw_elements().back().getValue();
			const fieldType current_max_coefficient = source.raw_elements().back().getCoefficient();

			bool found = false;
			DirectionType best_lift{};
			fieldType best_lift_coefficient = fieldType{0};
			OCGraph best_residual;

			for (const Element& elem : source.raw_elements()) {
				for (Int i = 0; i < GraphType::SIZE; ++i) {
					if (elem.getValue()[i]) {
						continue;
					}

					DirectionType lifted = elem.getValue();
					lifted[i] = true;

					OCGraph lifted_graph(graph, single_direction_lincomb(lifted, fieldType{1}), automorphisms);
					OCGraph lifted_boundary = lifted_graph.delta_e();

					const fieldType boundary_max_coefficient =
						coefficient_of_direction(static_cast<const DirectionComb&>(lifted_boundary), current_max);
					if (boundary_max_coefficient == fieldType{0}) {
						continue;
					}

					const fieldType lift_coefficient = -current_max_coefficient / boundary_max_coefficient;
					DirectionComb residual_combination = static_cast<const DirectionComb&>(source);
					residual_combination += static_cast<const DirectionComb&>(lifted_boundary) * lift_coefficient;

					OCGraph residual(graph, std::move(residual_combination), automorphisms);
					residual.standardize_and_sort();

					if (compare_descending_direction_support(
						static_cast<const DirectionComb&>(residual),
						static_cast<const DirectionComb&>(source)
					) >= 0) {
						continue;
					}

					if (!found || compare_descending_direction_support(
						static_cast<const DirectionComb&>(residual),
						static_cast<const DirectionComb&>(best_residual)
					) < 0) {
						found = true;
						best_lift = lifted;
						best_lift_coefficient = lift_coefficient;
						best_residual = std::move(residual);
					}
				}
			}

			if (!found) {
				return std::nullopt;
			}

			OCGraph primitive_step(graph, single_direction_lincomb(best_lift, best_lift_coefficient), automorphisms);
			primitive_step.standardize_and_sort();
			return primitive_step;
		}

		std::optional<OCGraph> greedy_delta_e_primitive(bigInt max_steps = 10000) const {
			(void)max_steps;
			return morse_delta_e_primitive();
		}

		std::optional<OCGraph> exact_delta_e_primitive() const {
			OCGraph source = normalized_copy();
			if (source.size() == 0) {
				return OCGraph(graph, DirectionComb{}, automorphisms);
			}

			if (all_terms_have_exactly_one_true_half_edge_at_each_vertex(static_cast<const DirectionComb&>(source))) {
				OCGraph correction = source.homotopy_standardize();
				DirectionComb residual = static_cast<const DirectionComb&>(source);
				residual += static_cast<const DirectionComb&>(correction.delta_e());
				OCGraph reduced(graph, std::move(residual), automorphisms);
				reduced.standardize_and_sort();
				if (reduced.size() == 0) {
					correction.standardize_and_sort();
					return correction;
				}
			}

			std::map<Int, DirectionComb> by_true_count;
			for (const Element& elem : source.raw_elements()) {
				by_true_count[count_true(elem.getValue())] .append_in_basis_order(elem);
			}

			DirectionComb primitive_combination;

			for (const auto& [true_count, component] : by_true_count) {
				const auto lifted_component = exact_delta_e_primitive_for_true_count(component, true_count);
				if (!lifted_component.has_value()) {
					return std::nullopt;
				}
				static_cast<DirectionComb&>(primitive_combination) += lifted_component.value();
			}

			OCGraph primitive(graph, std::move(primitive_combination), automorphisms);
			primitive.standardize_and_sort();

			DirectionComb check = static_cast<const DirectionComb&>(primitive.delta_e());
			DirectionComb target = static_cast<const DirectionComb&>(source);
			check.sort_elements();
			target.sort_elements();
			if (!(check == target)) {
				return std::nullopt;
			}

			return primitive;
		}

		std::optional<OCGraph> morse_delta_e_homotopy() const {
			OCGraph source = normalized_copy();
			OCGraph correction(graph, DirectionComb{}, automorphisms);

			const DirectionType critical_direction = minimal_one_true_direction();
			const auto critical_indices = critical_true_indices_by_vertex(critical_direction);

			for (const Element& elem : source.raw_elements()) {
				const auto lifted = morse_delta_e_homotopy_of_element(elem, critical_indices);
				if (!lifted.has_value()) {
					return std::nullopt;
				}
				static_cast<DirectionComb&>(correction) += lifted.value();
			}

			correction.standardize_and_sort();
			return correction;
		}

		std::optional<OCGraph> morse_delta_e_primitive(bigInt max_steps = 10000) const {
			OCGraph source = normalized_copy();
			if (source.size() == 0) {
				return OCGraph(graph, DirectionComb{}, automorphisms);
			}

			if (source.delta_e().size() != 0) {
				return std::nullopt;
			}

			const DirectionType critical_direction = minimal_one_true_direction();

			OCGraph primitive(graph, DirectionComb{}, automorphisms);
			OCGraph residual = source;

			for (bigInt step = 0; step < max_steps && residual.size() != 0; ++step) {
				const auto homotopy_step = residual.morse_delta_e_homotopy();
				if (!homotopy_step.has_value()) {
					return std::nullopt;
				}

				if (homotopy_step->size() == 0) {
					break;
				}

				static_cast<DirectionComb&>(primitive) +=
					static_cast<const DirectionComb&>(homotopy_step.value());

				DirectionComb next_residual = static_cast<const DirectionComb&>(residual);
				next_residual += static_cast<const DirectionComb&>(homotopy_step.value().delta_e()) * fieldType{-1};
				OCGraph reduced(graph, std::move(next_residual), automorphisms);
				reduced.standardize_and_sort();
				residual = std::move(reduced);
			}

			if (residual.size() == 0) {
				primitive.standardize_and_sort();
				return primitive;
			}

			if (all_terms_have_exactly_one_true_half_edge_at_each_vertex(
				static_cast<const DirectionComb&>(residual)
			)) {
				OCGraph correction = residual.homotopy_standardize();
				DirectionComb final_residual = static_cast<const DirectionComb&>(residual);
				final_residual += static_cast<const DirectionComb&>(correction.delta_e()) * fieldType{-1};
				OCGraph reduced(graph, std::move(final_residual), automorphisms);
				reduced.standardize_and_sort();

				if (reduced.size() == 0) {
					static_cast<DirectionComb&>(primitive) +=
						static_cast<const DirectionComb&>(correction);
					primitive.standardize_and_sort();
					return primitive;
				}

				residual = std::move(reduced);
			}

			if (residual.size() == 1
				&& residual.raw_elements().front().getValue() == critical_direction) {
				return std::nullopt;
			}

			return std::nullopt;
		}

		OCGraph homotopy_standardize() const {
			OCGraph source = normalized_copy();
			OCGraph correction(graph, DirectionComb{}, automorphisms);

			const DirectionType d0 = minimal_one_true_direction();

			for (const Element& elem : source.raw_elements()) {
				if (!has_exactly_one_true_half_edge_at_each_vertex(elem.getValue())) {
					continue;
				}

				DirectionComb term_correction = homotopy_standardize_one_true_direction(
					elem.getValue(),
					d0,
					elem.getCoefficient()
				);
				static_cast<DirectionComb&>(correction) += term_correction;
			}

			correction.normalize_in_place();
			return correction;
		}

		void standardize_and_sort() {
			standardize_graph_only();

			standardize_directions();
			this->sort_elements();
		}

		GraphType graph{};
		std::vector<IsoType> automorphisms{};

	private:
		Int count_true(const DirectionType& directions) const {
			Int result = 0;
			for (Int i = 0; i < GraphType::SIZE; ++i) {
				if (directions[i]) {
					++result;
				}
			}
			return result;
		}

		std::vector<DirectionType> enumerate_canonical_admissible_directions_with_true_count(Int true_count) const {
			std::vector<DirectionType> raw_directions;
			DirectionType current;
			current.fill(false);

			auto backtrack = [&](auto&& self, Int index, Int chosen) -> void {
				if (chosen > true_count) {
					return;
				}
				if (chosen + (GraphType::SIZE - index) < true_count) {
					return;
				}
				if (index == GraphType::SIZE) {
					if (chosen != true_count || !has_true_half_edge_at_each_vertex(current)) {
						return;
					}
					const Element canon = standardize_direction_element(Element(current, fieldType{1}));
					if (canon.getCoefficient() != fieldType{0}) {
						raw_directions.push_back(canon.getValue());
					}
					return;
				}

				current[index] = false;
				self(self, index + 1, chosen);

				current[index] = true;
				self(self, index + 1, chosen + 1);
				current[index] = false;
			};

			backtrack(backtrack, 0, 0);
			std::sort(raw_directions.begin(), raw_directions.end());
			raw_directions.erase(
				std::unique(raw_directions.begin(), raw_directions.end()),
				raw_directions.end()
			);
			return raw_directions;
		}

		std::optional<DirectionComb> exact_delta_e_primitive_for_true_count(
			const DirectionComb& target_component,
			Int true_count
		) const {
			const std::vector<DirectionType> row_basis =
				enumerate_canonical_admissible_directions_with_true_count(true_count);
			const std::vector<DirectionType> col_basis =
				enumerate_canonical_admissible_directions_with_true_count(true_count + 1);

			if (row_basis.empty()) {
				return target_component.size() == 0 ? std::optional<DirectionComb>(DirectionComb{}) : std::nullopt;
			}

			std::map<DirectionType, Int> row_index;
			for (Int i = 0; i < static_cast<Int>(row_basis.size()); ++i) {
				row_index.emplace(row_basis[i], i);
			}

			std::vector<std::vector<fieldType>> matrix(
				row_basis.size(),
				std::vector<fieldType>(col_basis.size() + 1, fieldType{0})
			);

			for (const Element& elem : target_component.raw_elements()) {
				const auto it = row_index.find(elem.getValue());
				if (it == row_index.end()) {
					return std::nullopt;
				}
				matrix[it->second][col_basis.size()] = elem.getCoefficient();
			}

			for (Int col = 0; col < static_cast<Int>(col_basis.size()); ++col) {
				OCGraph lifted_graph(graph, single_direction_lincomb(col_basis[col], fieldType{1}), automorphisms);
				const OCGraph boundary = lifted_graph.delta_e();
				for (const Element& boundary_term : boundary.raw_elements()) {
					const auto it = row_index.find(boundary_term.getValue());
					if (it != row_index.end()) {
						matrix[it->second][col] = boundary_term.getCoefficient();
					}
				}
			}

			const std::vector<fieldType> solution = solve_linear_system(matrix, static_cast<Int>(col_basis.size()));
			if (solution.empty() && !col_basis.empty()) {
				return std::nullopt;
			}

			DirectionComb primitive_component;
			for (Int col = 0; col < static_cast<Int>(col_basis.size()); ++col) {
				if (solution[col] != fieldType{0}) {
					primitive_component.append_in_basis_order(col_basis[col], solution[col]);
				}
			}
			primitive_component.sort_elements();
			return primitive_component;
		}

		std::vector<fieldType> solve_linear_system(
			std::vector<std::vector<fieldType>>& augmented,
			Int n_variables
		) const {
			const Int n_rows = static_cast<Int>(augmented.size());
			Int pivot_row = 0;
			std::vector<Int> pivot_column_for_row;
			pivot_column_for_row.reserve(std::min(n_rows, n_variables));

			for (Int col = 0; col < n_variables && pivot_row < n_rows; ++col) {
				Int selected_row = pivot_row;
				while (selected_row < n_rows && augmented[selected_row][col] == fieldType{0}) {
					++selected_row;
				}
				if (selected_row == n_rows) {
					continue;
				}

				if (selected_row != pivot_row) {
					std::swap(augmented[selected_row], augmented[pivot_row]);
				}

				const fieldType inv_pivot = fieldType{1} / augmented[pivot_row][col];
				for (Int j = col; j <= n_variables; ++j) {
					augmented[pivot_row][j] *= inv_pivot;
				}

				for (Int row = 0; row < n_rows; ++row) {
					if (row == pivot_row || augmented[row][col] == fieldType{0}) {
						continue;
					}
					const fieldType factor = augmented[row][col];
					for (Int j = col; j <= n_variables; ++j) {
						augmented[row][j] -= factor * augmented[pivot_row][j];
					}
				}

				pivot_column_for_row.push_back(col);
				++pivot_row;
			}

			for (Int row = pivot_row; row < n_rows; ++row) {
				if (augmented[row][n_variables] != fieldType{0}) {
					return {};
				}
			}

			std::vector<fieldType> solution(n_variables, fieldType{0});
			for (Int row = 0; row < static_cast<Int>(pivot_column_for_row.size()); ++row) {
				solution[pivot_column_for_row[row]] = augmented[row][n_variables];
			}
			return solution;
		}

		static signedInt compare_descending_direction_support(
			const DirectionComb& lhs,
			const DirectionComb& rhs
		) {
			const auto& left = lhs.raw_elements();
			const auto& right = rhs.raw_elements();

			bigInt i = static_cast<bigInt>(left.size());
			bigInt j = static_cast<bigInt>(right.size());

			while (i > 0 && j > 0) {
				const signedInt comparison =
					left[i - 1].getValue().compare(right[j - 1].getValue());
				if (comparison != 0) {
					return comparison;
				}
				--i;
				--j;
			}

			if (i == 0 && j == 0) {
				return 0;
			}
			return i == 0 ? -1 : 1;
		}

		fieldType coefficient_of_direction(
			const DirectionComb& directions,
			const DirectionType& target
		) const {
			for (const Element& elem : directions.raw_elements()) {
				if (elem.getValue() == target) {
					return elem.getCoefficient();
				}
			}
			return fieldType{0};
		}

		signedInt graph_isomorphism_sign(const IsoType& iso) const {
			return iso.graph_permutation_sign() * iso.vertex_permutation_sign();
		}

		signedInt direction_isomorphism_sign(const IsoType& iso, const DirectionType& directions) const {
			return graph_isomorphism_sign(iso) * iso.direction_sign(directions);
		}

		Element standardize_direction_element(const Element& input) const {
			if (automorphisms.empty()) {
				return input;
			}

			const IsoType& first_iso = automorphisms.front();
			DirectionType best_value = first_iso.permute(input.getValue());
			const signedInt first_sign = direction_isomorphism_sign(first_iso, input.getValue());
			bool contains_plus = first_sign > 0;
			bool contains_minus = first_sign < 0;

			for (bigInt i = 1; i < automorphisms.size(); ++i) {
				const IsoType& iso = automorphisms[i];
				const DirectionType candidate = iso.permute(input.getValue());
				const signedInt candidate_sign = direction_isomorphism_sign(iso, input.getValue());
				const signedInt comparison = candidate.compare(best_value);

				if (comparison < 0) {
					best_value = candidate;
					contains_plus = candidate_sign > 0;
					contains_minus = candidate_sign < 0;
					continue;
				}

				if (comparison == 0) {
					contains_plus = contains_plus || (candidate_sign > 0);
					contains_minus = contains_minus || (candidate_sign < 0);
				}
			}

			if (contains_plus && contains_minus) {
				return Element(best_value, fieldType{0});
			}

			const fieldType coeff = contains_plus ? input.getCoefficient() : -input.getCoefficient();
			return Element(best_value, coeff);
		}

		OCGraph normalized_copy() const {
			OCGraph copy = *this;
			copy.normalize_in_place();
			return copy;
		}

		void normalize_in_place() {
			if (automorphisms.empty()) {
				this->sort_elements();
			} else {
				standardize_directions();
				this->sort_elements();
			}
		}

		void standardize_graph_only() {
			if (!automorphisms.empty()) {
				return;
			}

			Standardizer standardizer;
			const std::vector<IsoType> minimizing_isomorphisms = standardizer.minimizing_isomorphisms(graph);

			if (minimizing_isomorphisms.empty()) {
				return;
			}

			const IsoType sigma0 = minimizing_isomorphisms.front();
			apply_isomorphism_to_graph_and_directions(sigma0);

			const IsoType sigma0_inverse = sigma0.inverse();
			automorphisms.reserve(minimizing_isomorphisms.size());
			for (const IsoType& sigma : minimizing_isomorphisms) {
				automorphisms.push_back(sigma0_inverse.compose(sigma));
			}
		}

		Element minimal_one_true_direction_element() const {
			DirectionType result;
			result.fill(false);
			std::array<bool, GraphType::N_VERTICES_> already_has_open{};
			std::array<Int, GraphType::N_VERTICES_> perm{};
			std::array<Int, GraphType::SIZE> current_order{};
			signedInt sign = 1;

			for (Int i = 0; i < GraphType::SIZE; ++i) {
				current_order[i] = i;
			}

			for (Int target_slot = 0; target_slot < GraphType::N_VERTICES_; ++target_slot) {
				Int current_pos = target_slot;
				for (; current_pos < GraphType::SIZE; ++current_pos) {
					const Int original_index = current_order[current_pos];
					const Int vertex = graph.half_edges[original_index];
					if (!already_has_open[vertex]) {
						already_has_open[vertex] = true;
						result[original_index] = true;
						perm[target_slot] = vertex;
						break;
					}
				}

				for (Int passed_pos = target_slot; passed_pos < current_pos; ++passed_pos) {
					const bool passed_has_odd_degree =
						((passed_pos % 2) == 0) ? ((GraphType::C_ % 2) != 0) : ((GraphType::D_ % 2) != 0);
					if (passed_has_odd_degree) {
						sign = -sign;
					}
				}

				const Int chosen_index = current_order[current_pos];
				for (Int pos = current_pos; pos > target_slot; --pos) {
					current_order[pos] = current_order[pos - 1];
				}
				current_order[target_slot] = chosen_index;
			}

			sign *= Permutation<GraphType::N_VERTICES_>(perm).sign();

			const Element base(result, fieldType{sign});

			if (automorphisms.empty()) {
				return base;
			}

			return standardize_direction_element(base);
		}

		DirectionType minimal_one_true_direction() const {
			return minimal_one_true_direction_element().getValue();
		}

		DirectionComb single_direction_lincomb(const DirectionType& directions, fieldType coefficient) const {
			DirectionComb result;
			result.append_in_basis_order(directions, coefficient);
			result.sort_elements();
			return result;
		}

		std::array<Int, GraphType::N_VERTICES_> true_counts_by_vertex(const DirectionType& directions) const {
			std::array<Int, GraphType::N_VERTICES_> counts{};
			for (Int i = 0; i < GraphType::SIZE; ++i) {
				if (directions[i]) {
					++counts[graph.half_edges[i]];
				}
			}
			return counts;
		}

		std::array<Int, GraphType::N_VERTICES_> critical_true_indices_by_vertex(const DirectionType& directions) const {
			std::array<Int, GraphType::N_VERTICES_> indices{};
			for (Int i = 0; i < GraphType::SIZE; ++i) {
				if (directions[i]) {
					indices[graph.half_edges[i]] = i;
				}
			}
			return indices;
		}

		std::optional<DirectionComb> morse_delta_e_homotopy_of_element(
			const Element& elem,
			const std::array<Int, GraphType::N_VERTICES_>& critical_indices
		) const {
			const DirectionType& directions = elem.getValue();
			const auto counts = true_counts_by_vertex(directions);

			for (Int vertex = 0; vertex < GraphType::N_VERTICES_; ++vertex) {
				const Int critical_index = critical_indices[vertex];
				const bool has_critical = directions[critical_index];
				const Int count = counts[vertex];

				if (has_critical && count == 1) {
					continue;
				}

				if (has_critical) {
					return DirectionComb{};
				}

				DirectionType lifted = directions;
				lifted[critical_index] = true;

				OCGraph lifted_graph(graph, single_direction_lincomb(lifted, fieldType{1}), automorphisms);
				lifted_graph.standardize_and_sort();
				if (lifted_graph.size() == 0) {
					return std::nullopt;
				}

				const fieldType boundary_coefficient = coefficient_of_direction(
					static_cast<const DirectionComb&>(lifted_graph.delta_e()),
					directions
				);
				if (boundary_coefficient == fieldType{0}) {
					return std::nullopt;
				}

				DirectionComb result = static_cast<const DirectionComb&>(lifted_graph);
				result *= elem.getCoefficient() / boundary_coefficient;
				return result;
			}

			return DirectionComb{};
		}

		Int unique_true_index_at_vertex(const DirectionType& directions, Int vertex) const {
			for (Int i = 0; i < GraphType::SIZE; ++i) {
				if (graph.half_edges[i] == vertex && directions[i]) {
					return i;
				}
			}
			return 0;
		}

		DirectionComb homotopy_standardize_one_true_direction(
			DirectionType current,
			const DirectionType& d0,
			fieldType current_coefficient
		) const {
			DirectionComb correction;

			while (!(current == d0)) {
				Int target_index = 0;
				for (Int i = 0; i < GraphType::SIZE; ++i) {
					if (d0[i] && !current[i]) {
						target_index = i;
						break;
					}
				}

				const Int vertex = graph.half_edges[target_index];
				const Int source_index = unique_true_index_at_vertex(current, vertex);

				DirectionType lifted = current;
				lifted[target_index] = true;

				OCGraph lifted_graph(graph, single_direction_lincomb(lifted, fieldType{1}), automorphisms);
				OCGraph lifted_boundary = lifted_graph.delta_e();

				fieldType current_boundary_coeff = fieldType{0};
				fieldType next_boundary_coeff = fieldType{0};

				DirectionType next = current;
				next[target_index] = true;
				next[source_index] = false;
				if (!automorphisms.empty()) {
					next = DirectionType::standardize(Element(next, fieldType{1}), automorphisms).getValue();
				}

				for (const Element& boundary_term : lifted_boundary.raw_elements()) {
					if (boundary_term.getValue() == current) {
						current_boundary_coeff = boundary_term.getCoefficient();
					}
					if (boundary_term.getValue() == next) {
						next_boundary_coeff = boundary_term.getCoefficient();
					}
				}

				const fieldType step_coefficient = -current_coefficient / current_boundary_coeff;
				correction.append_in_basis_order(lifted, step_coefficient);
				current_coefficient = step_coefficient * next_boundary_coeff;
				current = next;
			}

			correction.sort_elements();
			return correction;
		}

		void apply_isomorphism_to_graph_and_directions(const IsoType& iso) {
			graph = iso.permute(graph);

			for (Element& elem : this->raw_elements_nonconst()) {
				elem = Element(
					iso.permute(elem.getValue()),
					elem.getCoefficient() * fieldType{
						direction_isomorphism_sign(iso, elem.getValue())
					}
				);
			}
		}

		void standardize_directions() {
			std::vector<Element> standardized;
			standardized.reserve(this->size());

			for (const Element& elem : this->raw_elements()) {
				if (!has_true_half_edge_at_each_vertex(elem.getValue())) {
					continue;
				}

				Element canon = standardize_direction_element(elem);
				if (canon.getCoefficient() != fieldType{0} && has_true_half_edge_at_each_vertex(canon.getValue())) {
					standardized.push_back(std::move(canon));
				}
			}

			this->raw_elements_nonconst() = std::move(standardized);
		}

		void standardize_directions_without_vertex_check() {
			std::vector<Element> standardized;
			standardized.reserve(this->size());

			for (const Element& elem : this->raw_elements()) {
				Element canon = standardize_direction_element(elem);
				if (canon.getCoefficient() != fieldType{0}) {
					standardized.push_back(std::move(canon));
				}
			}

			this->raw_elements_nonconst() = std::move(standardized);
		}

		template <typename SomeGraphType>
		static bool has_true_half_edge_at_each_vertex(
			const SomeGraphType& some_graph,
			const GraphDirections<SomeGraphType>& directions
		) {
			std::array<bool, SomeGraphType::N_VERTICES_> seen{};

			for (Int i = 0; i < SomeGraphType::SIZE; ++i) {
				if (directions[i]) {
					seen[some_graph.half_edges[i]] = true;
				}
			}

			return std::all_of(seen.begin(), seen.end(), [](bool value) { return value; });
		}

		bool has_true_half_edge_at_each_vertex(const DirectionType& directions) const {
			return has_true_half_edge_at_each_vertex(graph, directions);
		}

		bool has_exactly_one_true_half_edge_at_each_vertex(const DirectionType& directions) const {
			std::array<Int, GraphType::N_VERTICES_> counts{};

			for (Int i = 0; i < GraphType::SIZE; ++i) {
				if (directions[i]) {
					++counts[graph.half_edges[i]];
				}
			}

			return std::all_of(counts.begin(), counts.end(), [](Int count) { return count == 1; });
		}

		bool all_terms_have_exactly_one_true_half_edge_at_each_vertex(const DirectionComb& directions) const {
			for (const Element& elem : directions.raw_elements()) {
				if (!has_exactly_one_true_half_edge_at_each_vertex(elem.getValue())) {
					return false;
				}
			}
			return true;
		}
};
