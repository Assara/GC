#pragma once

#include <algorithm>
#include <array>
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
};
