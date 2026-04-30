#pragma once

#include <iostream>
#include <vector>

#include "GraphDirections.hpp"
#include "GraphIsomorphism.hpp"
#include "OCGC.hpp"
#include "OCGraph.hpp"

template <typename GraphType>
bool check_graph_isomorphism_action(const char* label) {
	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	GraphIsomorphism<3, 3> iso;
	iso.vertex_permutation_data() = {1, 2, 0};
	iso.edge_permutation_data() = {2, 0, 1};
	iso.edge_flip_data() = {false, true, false};
	iso.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();

	GraphType image = iso.permute(source);

	GraphType expected;
	expected.setEdge(0, 0, 2);
	expected.setEdge(1, 1, 0);
	expected.setEdge(2, 1, 2);

	const bool ok = (image == expected);
	std::cout << label << ": graph isomorphism action -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_graph_isomorphism_composition(const char* label) {
	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	GraphIsomorphism<3, 3> a;
	a.vertex_permutation_data() = {1, 2, 0};
	a.edge_permutation_data() = {2, 0, 1};
	a.edge_flip_data() = {false, true, false};
	a.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();

	GraphIsomorphism<3, 3> b;
	b.vertex_permutation_data() = {2, 0, 1};
	b.edge_permutation_data() = {1, 2, 0};
	b.edge_flip_data() = {true, false, true};
	b.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();

	const GraphType stepwise = b.permute(a.permute(source));
	const GraphIsomorphism<3, 3> ab = a.compose(b);
	const GraphType composed = ab.permute(source);

	const bool ok = (stepwise == composed);
	std::cout << label << ": graph isomorphism composition -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_graph_isomorphism_inverse(const char* label) {
	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	GraphIsomorphism<3, 3> iso;
	iso.vertex_permutation_data() = {1, 2, 0};
	iso.edge_permutation_data() = {2, 0, 1};
	iso.edge_flip_data() = {false, true, false};
	iso.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();

	const auto inv = iso.inverse();
	const GraphType recovered = inv.permute(iso.permute(source));
	const GraphType recovered_other_side = iso.permute(inv.permute(source));

	const GraphIsomorphism<3, 3> id_left = iso.compose(inv);
	const GraphIsomorphism<3, 3> id_right = inv.compose(iso);
	const GraphType identity_left_image = id_left.permute(source);
	const GraphType identity_right_image = id_right.permute(source);

	const bool ok =
		(recovered == source) &&
		(recovered_other_side == source) &&
		(identity_left_image == source) &&
		(identity_right_image == source);

	std::cout << label << ": graph isomorphism inverse -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_graph_isomorphism_graph_sign(const char* label) {
	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	GraphIsomorphism<3, 3> iso;
	iso.edge_permutation_data() = {1, 0, 2};
	iso.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();

	const bool ok = (iso.graph_sign(source) == -1);
	std::cout << label << ": graph isomorphism graph sign -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_minimizing_isomorphisms(const char* label, const GraphType& source, bigInt expected_count) {
	GraphStandardizer<
		GraphType::N_VERTICES_,
		GraphType::N_EDGES_,
		GraphType::N_OUT_HAIR_,
		GraphType::N_IN_HAIR_,
		GraphType::C_,
		GraphType::D_,
		fieldType
	> standardizer;

	const auto minimizers = standardizer.minimizing_isomorphisms(source);
	bool ok = !minimizers.empty() && minimizers.size() == expected_count;

	if (!minimizers.empty()) {
		const GraphType target = minimizers.front().permute(source);
		for (const auto& iso : minimizers) {
			if (!(iso.permute(source) == target)) {
				ok = false;
				break;
			}
		}
	}

	const bool count_matches = standardizer.automorphism_group_size(source) == minimizers.size();
	ok &= count_matches;

	std::cout << label
		  << ": minimizing isomorphisms -> "
		  << (ok ? "ok" : "failed")
		  << " (count=" << minimizers.size() << ")\n";
	return ok;
}

template <typename GraphType>
bool check_minimizing_isomorphisms_match_standardization(const char* label, const GraphType& source, bigInt expected_count) {
	GraphStandardizer<
		GraphType::N_VERTICES_,
		GraphType::N_EDGES_,
		GraphType::N_OUT_HAIR_,
		GraphType::N_IN_HAIR_,
		GraphType::C_,
		GraphType::D_,
		fieldType
	> standardizer;

	const auto minimizers = standardizer.minimizing_isomorphisms(source);
	typename GraphType::Basis standardized_input(source, fieldType{1});
	GraphType::std(standardized_input);
	const GraphType standardized_graph = standardized_input.getValue();

	bool ok = minimizers.size() == expected_count;
	for (const auto& iso : minimizers) {
		if (!(iso.permute(source) == standardized_graph)) {
			ok = false;
			break;
		}
	}

	std::cout << label
		  << ": minimizing isomorphisms match std -> "
		  << (ok ? "ok" : "failed")
		  << " (count=" << minimizers.size() << ")\n";
	return ok;
}

template <typename GraphType>
bool check_generated_isomorphism_standardization_matches_standardize(const char* label, const GraphType& source) {
	GraphStandardizer<
		GraphType::N_VERTICES_,
		GraphType::N_EDGES_,
		GraphType::N_OUT_HAIR_,
		GraphType::N_IN_HAIR_,
		GraphType::C_,
		GraphType::D_,
		fieldType
	> standardizer;

	typename GraphType::Basis standardized_input(source, fieldType{1});
	GraphType::std(standardized_input);

	const auto minimizers = standardizer.minimizing_isomorphisms(source);
	bool ok = !minimizers.empty();

	if (ok) {
		const auto& sigma0 = minimizers.front();
		bool contains_plus = false;
		bool contains_minus = false;

		for (const auto& sigma : minimizers) {
			if (sigma.graph_permutation_sign() > 0) {
				contains_plus = true;
			}
			if (sigma.graph_permutation_sign() < 0) {
				contains_minus = true;
			}
		}

		const fieldType coeff =
			(contains_plus && contains_minus) ? fieldType{0}
			: contains_plus ? fieldType{1}
			: fieldType{-1};

		typename GraphType::Basis via_isomorphism(sigma0.permute(source), coeff);
		ok = standardized_input.total_equality(via_isomorphism);
	}

	std::cout << label
		  << ": generated isomorphism standardization -> "
		  << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_graph_directions_shape(const char* label) {
	GraphDirections<GraphType> directions;
	directions.fill(false);
	directions[0] = true;
	directions[GraphType::SIZE - 1] = true;

	const bool ok =
		(directions.size() == GraphType::SIZE) &&
		directions[0] &&
		directions[GraphType::SIZE - 1];

	std::cout << label << ": graph directions shape -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_graph_directions_isomorphism_action(const char* label) {
	GraphDirections<GraphType> source;
	source.fill(false);
	source[0] = true;
	source[GraphType::N_HAIR + 0] = true;
	source[GraphType::N_HAIR + 3] = true;

	GraphIsomorphism<3, 3> iso;
	iso.vertex_permutation_data() = {1, 2, 0};
	iso.edge_permutation_data() = {2, 0, 1};
	iso.edge_flip_data() = {false, true, false};
	iso.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();

	const auto image = iso.permute(source);

	GraphDirections<GraphType> expected;
	expected.fill(false);
	expected[GraphType::N_HAIR + 0] = true;
	expected[GraphType::N_HAIR + 4] = true;

	const bool ok = (image == expected);
	std::cout << label << ": graph directions isomorphism action -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_graph_directions_order(const char* label) {
	GraphDirections<GraphType> a;
	GraphDirections<GraphType> b;
	GraphDirections<GraphType> c;

	a.fill(false);
	b.fill(false);
	c.fill(false);

	b[GraphType::SIZE - 1] = true;
	c[0] = true;

	const bool ok =
		(a.compare(a) == 0) &&
		(a < b) &&
		(b < c) &&
		(a < c) &&
		!(b < a) &&
		!(c < b);

	std::cout << label << ": graph directions order -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_graph_directions_standardize(const char* label) {
	using DirectionType = GraphDirections<GraphType>;
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	DirectionType directions;
	directions.fill(false);
	directions[GraphType::N_HAIR + 2] = true;
	directions[GraphType::N_HAIR + 3] = true;

	IsoType identity;
	IsoType flip_middle;
	flip_middle.edge_flip_data()[1] = true;
	identity.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();
	flip_middle.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();

	const BasisElement<DirectionType, fieldType> input(directions, fieldType{1});
	const std::vector<IsoType> automorphisms{identity, flip_middle};
	const auto standardized = DirectionType::standardize(input, automorphisms);

	const bool ok = (standardized.getCoefficient() == fieldType{0});
	std::cout << label << ": graph directions standardize -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_ocgraph_standardize_and_sort(const char* label) {
	using DirectionType = GraphDirections<GraphType>;
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

	GraphType source;
	source.setEdge(0, 1, 0);
	source.setEdge(1, 2, 1);
	source.setEdge(2, 2, 0);

	DirectionType directions;
	directions.fill(false);
	directions[GraphType::N_HAIR + 0] = true;
	directions[GraphType::N_HAIR + 3] = true;

	typename OCGraph<GraphType>::DirectionComb lincomb;
	lincomb.append_in_basis_order(directions, fieldType{1});

	OCGraph<GraphType> oc(source, lincomb);
	oc.standardize_and_sort();

	Standardizer standardizer;
	const std::vector<IsoType> minimizing_isomorphisms = standardizer.minimizing_isomorphisms(source);
	const IsoType sigma0 = minimizing_isomorphisms.front();
	const GraphType expected_graph = sigma0.permute(source);

	const DirectionType moved_direction = sigma0.permute(directions);
	const fieldType moved_coeff = fieldType{sigma0.graph_sign(source) * sigma0.direction_sign(directions)};
	const BasisElement<DirectionType, fieldType> moved_basis(moved_direction, moved_coeff);

	std::vector<IsoType> expected_automorphisms;
	expected_automorphisms.reserve(minimizing_isomorphisms.size());
	const IsoType sigma0_inverse = sigma0.inverse();
	for (const IsoType& sigma : minimizing_isomorphisms) {
		expected_automorphisms.push_back(sigma0_inverse.compose(sigma));
	}

	const auto expected_direction_basis = DirectionType::standardize(moved_basis, expected_automorphisms);

	bool ok = (oc.graph == expected_graph);
	ok &= (oc.automorphisms.size() == expected_automorphisms.size());
	ok &= (oc.size() == (expected_direction_basis.getCoefficient() == fieldType{0} ? 0 : 1));

	if (expected_direction_basis.getCoefficient() != fieldType{0} && oc.size() == 1) {
		const auto& be = oc.raw_elements().front();
		ok &= (be.getValue() == expected_direction_basis.getValue());
		ok &= (be.getCoefficient() == expected_direction_basis.getCoefficient());
	}

	std::cout << label << ": ocgraph standardize and sort -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_ocgraph_filters_non_covering_directions(const char* label) {
	using DirectionType = GraphDirections<GraphType>;

	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	DirectionType directions;
	directions.fill(false);
	directions[GraphType::N_HAIR + 0] = true;
	directions[GraphType::N_HAIR + 1] = true;

	typename OCGraph<GraphType>::DirectionComb lincomb;
	lincomb.append_in_basis_order(directions, fieldType{1});

	OCGraph<GraphType> oc(source, lincomb);
	oc.standardize_and_sort();

	const bool ok = oc.size() == 0;
	std::cout << label << ": ocgraph drops non-covering directions -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_ocgraph_delta_e(const char* label) {
	using DirectionType = GraphDirections<GraphType>;
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	IsoType iso;
	iso.edge_permutation_data() = {1, 0, 2};
	iso.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();

	DirectionType directions;
	directions.fill(false);
	directions[GraphType::N_HAIR + 0] = true;
	directions[GraphType::N_HAIR + 1] = true;
	directions[GraphType::N_HAIR + 2] = true;
	directions[GraphType::N_HAIR + 3] = true;
	directions[GraphType::N_HAIR + 4] = true;
	directions[GraphType::N_HAIR + 5] = true;

	typename OCGraph<GraphType>::DirectionComb lincomb;
	lincomb.append_in_basis_order(directions, fieldType{3});

	OCGraph<GraphType> oc(source, lincomb, {iso});
	OCGraph<GraphType> delta = oc.delta_e();

	std::vector<BasisElement<DirectionType, fieldType>> expected;
	const std::vector<IsoType> automorphisms{iso};
	fieldType sign = fieldType{1};
	for (Int offset = 0; offset < GraphType::SIZE; ++offset) {
		const Int i = GraphType::SIZE - 1 - offset;
		if (!directions[i]) {
			continue;
		}

		DirectionType next = directions;
		next[i] = false;
		expected.push_back(
			BasisElement<DirectionType, fieldType>(
				iso.permute(next),
				fieldType{3} * sign * fieldType{iso.graph_sign(source) * iso.direction_sign(next)}
			)
		);
		sign = -sign;
	}

	typename OCGraph<GraphType>::DirectionComb expected_lincomb;
	for (const auto& elem : expected) {
		if (elem.getCoefficient() != fieldType{0}) {
			expected_lincomb.append_in_basis_order(elem);
		}
	}
	expected_lincomb.sort_elements();

	bool ok = (delta.graph == source);
	ok &= (delta.automorphisms.size() == 1);
	ok &= (delta.automorphisms.front().edge_permutation_data() == iso.edge_permutation_data());
	ok &= static_cast<typename OCGraph<GraphType>::DirectionComb>(delta) == expected_lincomb;

	std::cout << label << ": ocgraph delta_e -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_ocgraph_delta_e_prunes_uncovered_vertices(const char* label) {
	using DirectionType = GraphDirections<GraphType>;

	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	DirectionType directions;
	directions.fill(false);
	directions[GraphType::N_HAIR + 0] = true;
	directions[GraphType::N_HAIR + 2] = true;
	directions[GraphType::N_HAIR + 5] = true;

	typename OCGraph<GraphType>::DirectionComb lincomb;
	lincomb.append_in_basis_order(directions, fieldType{1});

	OCGraph<GraphType> oc(source, lincomb);
	OCGraph<GraphType> delta = oc.delta_e();

	const bool ok = delta.size() == 0;
	std::cout << label << ": ocgraph delta_e prunes uncovered vertices -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_ocgraph_delta_e_squares_to_zero(const char* label) {
	using DirectionType = GraphDirections<GraphType>;

	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	DirectionType directions;
	directions.fill(false);
	directions[GraphType::N_HAIR + 0] = true;
	directions[GraphType::N_HAIR + 1] = true;
	directions[GraphType::N_HAIR + 2] = true;
	directions[GraphType::N_HAIR + 4] = true;

	typename OCGraph<GraphType>::DirectionComb lincomb;
	lincomb.append_in_basis_order(directions, fieldType{1});

	OCGraph<GraphType> oc(source, lincomb);
	OCGraph<GraphType> delta2 = oc.delta_e().delta_e();
	delta2.sort_elements();

	const bool ok = delta2.size() == 0;
	std::cout << label << ": ocgraph delta_e^2 -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_ocgraph_delta_e_isomorphism_equivariance_on_data(
	const char* label,
	const GraphType& source,
	const GraphDirections<GraphType>& directions,
	const std::vector<GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>>& isomorphisms
) {
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	typename OCGraph<GraphType>::DirectionComb lincomb;
	lincomb.append_in_basis_order(directions, fieldType{2});

	OCGraph<GraphType> original(source, lincomb);
	OCGraph<GraphType> delta_original = original.delta_e();

	bool ok = true;
	for (const IsoType& iso : isomorphisms) {
		typename OCGraph<GraphType>::DirectionComb moved_lincomb;
		moved_lincomb.append_in_basis_order(
			iso.permute(directions),
			fieldType{2} * fieldType{iso.graph_sign(source) * iso.direction_sign(directions)}
		);

		OCGraph<GraphType> moved_graph(iso.permute(source), moved_lincomb);
		OCGraph<GraphType> delta_moved = moved_graph.delta_e();

		typename OCGraph<GraphType>::DirectionComb expected_lincomb;
		for (const auto& elem : delta_original.raw_elements()) {
			expected_lincomb.append_in_basis_order(
				iso.permute(elem.getValue()),
				elem.getCoefficient() * fieldType{iso.graph_sign(source) * iso.direction_sign(elem.getValue())}
			);
		}

		OCGraph<GraphType> expected(iso.permute(source), expected_lincomb);
		expected.sort_elements();
		delta_moved.sort_elements();

		ok &= (delta_moved.graph == expected.graph)
			&& (static_cast<typename OCGraph<GraphType>::DirectionComb>(delta_moved)
				== static_cast<typename OCGraph<GraphType>::DirectionComb>(expected));
	}

	std::cout << label << ": ocgraph delta_e equivariant under isomorphism -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_ocgraph_delta_e_isomorphism_equivariance(const char* label) {
	using DirectionType = GraphDirections<GraphType>;
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	DirectionType directions;
	directions.fill(false);
	directions[GraphType::N_HAIR + 0] = true;
	directions[GraphType::N_HAIR + 2] = true;
	directions[GraphType::N_HAIR + 3] = true;
	directions[GraphType::N_HAIR + 5] = true;

	typename OCGraph<GraphType>::DirectionComb lincomb;
	lincomb.append_in_basis_order(directions, fieldType{2});

	OCGraph<GraphType> original(source, lincomb);
	OCGraph<GraphType> delta_original = original.delta_e();
	std::vector<IsoType> isomorphisms;

	IsoType identity;
	isomorphisms.push_back(identity);

	IsoType edge_swap;
	edge_swap.edge_permutation_data() = {1, 0, 2};
	edge_swap.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();
	isomorphisms.push_back(edge_swap);

	IsoType edge_flip;
	edge_flip.edge_flip_data() = {true, false, false};
	edge_flip.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();
	isomorphisms.push_back(edge_flip);

	IsoType cyclic_with_flip;
	cyclic_with_flip.vertex_permutation_data() = {1, 2, 0};
	cyclic_with_flip.edge_permutation_data() = {2, 0, 1};
	cyclic_with_flip.edge_flip_data() = {false, true, false};
	cyclic_with_flip.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();
	isomorphisms.push_back(cyclic_with_flip);

	IsoType odd_vertex_odd_edge;
	odd_vertex_odd_edge.vertex_permutation_data() = {1, 0, 2};
	odd_vertex_odd_edge.edge_permutation_data() = {0, 2, 1};
	odd_vertex_odd_edge.edge_flip_data() = {true, false, true};
	odd_vertex_odd_edge.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();
	isomorphisms.push_back(odd_vertex_odd_edge);

	return check_ocgraph_delta_e_isomorphism_equivariance_on_data(label, source, directions, isomorphisms);
}

inline bool check_ocgraph_OCG_F0_same_degree(const char* label) {
	using GraphType = Graph<3, 3, 0, 0, 0, 0, fieldType>;
	using DirectionType = GraphDirections<GraphType>;
	using IsoType = GraphIsomorphism<3, 3>;

	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	IsoType iso;
	iso.vertex_permutation_data() = {1, 2, 0};
	iso.edge_permutation_data() = {2, 0, 1};
	iso.edge_flip_data() = {false, true, false};
	iso.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();

	GraphType permuted = iso.permute(source);

	OCGraph<GraphType> image = OCGraph<GraphType>::OCG_F0(source);
	OCGraph<GraphType> permuted_image = OCGraph<GraphType>::OCG_F0(permuted);

	bool ok = (image == permuted_image);
	ok &= image.size() <= 1;

	if (image.size() == 1) {
		std::array<Int, GraphType::N_VERTICES_> counts{};
		const DirectionType& directions = image.raw_elements().front().getValue();
		for (Int i = 0; i < GraphType::SIZE; ++i) {
			if (directions[i]) {
				++counts[image.graph.half_edges[i]];
			}
		}
		ok &= std::all_of(counts.begin(), counts.end(), [](Int count) { return count == 1; });
	}

	std::cout << label << ": ocgraph OCG_F0 same degree -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename SourceGraphType>
bool check_OCG_F0_exact_equivariance_on_graph(
	const char* label,
	const SourceGraphType& source,
	const std::vector<GraphIsomorphism<SourceGraphType::N_VERTICES_, SourceGraphType::N_EDGES_>>& isomorphisms
) {
	using IsoType = GraphIsomorphism<SourceGraphType::N_VERTICES_, SourceGraphType::N_EDGES_>;

	const OCGraph<SourceGraphType> reference = OCGraph<SourceGraphType>::OCG_F0(source);

	bool ok = true;
	for (const IsoType& iso : isomorphisms) {
		const SourceGraphType moved = iso.permute(source);
		const OCGraph<SourceGraphType> image = OCGraph<SourceGraphType>::OCG_F0(moved);
		ok &= (image.graph == reference.graph)
			&& (static_cast<const typename OCGraph<SourceGraphType>::DirectionComb&>(image)
				== static_cast<const typename OCGraph<SourceGraphType>::DirectionComb&>(reference));
	}

	std::cout << label << ": OCG_F0 exact equivariance -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename SourceGraphType>
bool check_OCG_F0_unique_direction_on_graph(
	const char* label,
	const SourceGraphType& source,
	const std::vector<GraphIsomorphism<SourceGraphType::N_VERTICES_, SourceGraphType::N_EDGES_>>& isomorphisms
) {
	using TargetOCGraph = OCGraph<SourceGraphType>;
	using DirectionElement = typename TargetOCGraph::Element;
	using DirectionType = typename TargetOCGraph::DirectionType;
	using IsoType = GraphIsomorphism<SourceGraphType::N_VERTICES_, SourceGraphType::N_EDGES_>;

	const TargetOCGraph reference = TargetOCGraph::OCG_F0(source);
	bool ok = reference.size() == 1;

	if (!ok) {
		std::cout << label << ": OCG_F0 unique direction data -> failed\n";
		return false;
	}

	const DirectionElement& reference_element = reference.raw_elements().front();
	const DirectionType& reference_direction = reference_element.getValue();
	const fieldType reference_coefficient = reference_element.getCoefficient();

	for (const IsoType& iso : isomorphisms) {
		const SourceGraphType moved = iso.permute(source);
		const TargetOCGraph image = TargetOCGraph::OCG_F0(moved);

		ok &= image.graph == reference.graph;
		ok &= image.size() == 1;
		if (image.size() != 1) {
			continue;
		}

		const DirectionElement& image_element = image.raw_elements().front();
		const DirectionType& image_direction = image_element.getValue();
		ok &= image_element.getCoefficient() == reference_coefficient;

		for (Int i = 0; i < DirectionType::SIZE; ++i) {
			ok &= image_direction[i] == reference_direction[i];
		}
	}

	std::cout << label << ": OCG_F0 unique direction data -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
OCGC<GraphType> reduce_ocgc_via_delta_e_homotopy(const OCGC<GraphType>& input) {
	using OCGraphType = OCGraph<GraphType>;
	using OCGCType = OCGC<GraphType>;

	std::vector<typename OCGCType::Base> reduced_terms;
	reduced_terms.reserve(input.data().raw_elements().size());

	for (const auto& be : input.data().raw_elements()) {
		const OCGraphType& y = be.getValue();
		OCGraphType correction = y.homotopy_standardize();

		typename OCGraphType::DirectionComb reduced_base =
			static_cast<typename OCGraphType::DirectionComb>(y);
		reduced_base += static_cast<const typename OCGraphType::DirectionComb&>(correction.delta_e());

		OCGraphType reduced(y.graph, std::move(reduced_base), y.automorphisms);
		reduced.standardize_and_sort();
		if (reduced.size() != 0) {
			reduced_terms.emplace_back(std::move(reduced), be.getCoefficient());
		}
	}

	return OCGCType(std::move(reduced_terms));
}

template <typename GraphType>
OCGC<GraphType> absorb_outer_coefficients_into_ocgraphs(const OCGC<GraphType>& input) {
	using OCGraphType = OCGraph<GraphType>;
	using OCGCType = OCGC<GraphType>;

	std::vector<typename OCGCType::Base> normalized_terms;
	normalized_terms.reserve(input.data().raw_elements().size());

	for (const auto& be : input.data().raw_elements()) {
		typename OCGraphType::DirectionComb scaled_inner =
			static_cast<const typename OCGraphType::DirectionComb&>(be.getValue()) * be.getCoefficient();
		OCGraphType normalized(be.getValue().graph, std::move(scaled_inner), be.getValue().automorphisms);
		normalized.standardize_and_sort();
		if (normalized.size() != 0) {
			normalized_terms.emplace_back(std::move(normalized), fieldType{1});
		}
	}

	return OCGCType(std::move(normalized_terms));
}

template <typename SourceGraphType>
bool check_OCG_F0_delta_homotopy_equivalence_on_graph(
	const char* label,
	const SourceGraphType& source
) {
	using GCType = GC<
		SourceGraphType::N_VERTICES_,
		SourceGraphType::N_EDGES_,
		SourceGraphType::N_OUT_HAIR_,
		SourceGraphType::N_IN_HAIR_,
		SourceGraphType::C_,
		SourceGraphType::D_
	>;
	using SourceOCGraphType = OCGraph<SourceGraphType>;
	using SplitGraphType = typename SourceGraphType::SplitGraph;
	using SplitOCGraphType = OCGraph<SplitGraphType>;
	using SplitOCGCType = OCGC<SplitGraphType>;

	GCType gc(source);
	const auto dgc = gc.delta();

	std::vector<typename SplitOCGCType::Base> lhs_terms;
	lhs_terms.reserve(dgc.data().raw_elements().size());
	for (const auto& be : dgc.data().raw_elements()) {
		SplitOCGraphType image = SplitOCGraphType::OCG_F0(be.getValue());
		if (image.size() != 0) {
			lhs_terms.emplace_back(std::move(image), be.getCoefficient());
		}
	}
	SplitOCGCType lhs(std::move(lhs_terms));
	lhs.standardize_and_sort();

	OCGC<SourceGraphType> f0_of_g(SourceOCGraphType::OCG_F0(source));
	auto rhs = f0_of_g.delta_v();
	rhs.standardize_and_sort();
	auto rhs_reduced = reduce_ocgc_via_delta_e_homotopy(rhs);
	rhs_reduced.standardize_and_sort();

	auto lhs_flat = absorb_outer_coefficients_into_ocgraphs(lhs);
	auto rhs_flat = absorb_outer_coefficients_into_ocgraphs(rhs_reduced);
	lhs_flat.standardize_and_sort();
	rhs_flat.standardize_and_sort();

	const bool ok = lhs_flat.data() == rhs_flat.data();
	std::cout << label << ": F0(delta(G)) ~ delta_v(F0(G)) via delta_e -> "
		  << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_ocgraph_homotopy_standardize_minimal_true_case(const char* label) {
	using DirectionType = GraphDirections<GraphType>;

	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	DirectionType directions;
	directions.fill(false);
	directions[GraphType::N_HAIR + 0] = true;
	directions[GraphType::N_HAIR + 1] = true;
	directions[GraphType::N_HAIR + 3] = true;

	typename OCGraph<GraphType>::DirectionComb lincomb;
	lincomb.append_in_basis_order(directions, fieldType{1});

	OCGraph<GraphType> y(source, lincomb);
	OCGraph<GraphType> x = y.homotopy_standardize();

	typename OCGraph<GraphType>::DirectionComb reduced_base = static_cast<typename OCGraph<GraphType>::DirectionComb>(y);
	reduced_base += static_cast<const typename OCGraph<GraphType>::DirectionComb&>(x.delta_e());
	OCGraph<GraphType> reduced(source, reduced_base);
	reduced.sort_elements();

	DirectionType expected;
	expected.fill(false);
	expected[GraphType::N_HAIR + 0] = true;
	expected[GraphType::N_HAIR + 1] = true;
	expected[GraphType::N_HAIR + 3] = true;

	bool ok = reduced.size() == 1;
	if (ok) {
		ok &= reduced.raw_elements().front().getValue() == expected;
	}

	std::cout << label << ": ocgraph homotopy standardize -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_ocgraph_compare_and_merge(const char* label) {
	using DirectionType = GraphDirections<GraphType>;
	using OCBasis = typename OCGraph<GraphType>::Basis;

	GraphType source_a;
	source_a.setEdge(0, 0, 1);
	source_a.setEdge(1, 1, 2);
	source_a.setEdge(2, 0, 2);

	GraphType source_b = source_a;

	DirectionType d0;
	d0.fill(false);
	d0[GraphType::N_HAIR + 0] = true;
	d0[GraphType::N_HAIR + 3] = true;

	DirectionType d1;
	d1.fill(false);
	d1[GraphType::N_HAIR + 1] = true;
	d1[GraphType::N_HAIR + 4] = true;

	typename OCGraph<GraphType>::DirectionComb lincomb_a;
	lincomb_a.append_in_basis_order(d0, fieldType{2});

	typename OCGraph<GraphType>::DirectionComb lincomb_b;
	lincomb_b.append_in_basis_order(d1, fieldType{5});

	OCGraph<GraphType> oc_a(source_a, lincomb_a);
	OCGraph<GraphType> oc_b(source_b, lincomb_b);

	const bool compare_ok = (oc_a.compare(oc_b) == 0);

	OCBasis merged = OCGraph<GraphType>::merge(
		OCBasis(oc_a, fieldType{3}),
		OCBasis(oc_b, fieldType{-1})
	);

	bool merge_ok = (merged.getCoefficient() == fieldType{1});
	merge_ok &= (merged.getValue().compare(oc_a) == 0);
	merge_ok &= (merged.getValue().size() == 2);

	if (merged.getValue().size() == 2) {
		const auto& elems = merged.getValue().raw_elements();
		bool found_d0 = false;
		bool found_d1 = false;
		for (const auto& elem : elems) {
			if (elem.getValue() == d0 && elem.getCoefficient() == fieldType{6}) {
				found_d0 = true;
			}
			if (elem.getValue() == d1 && elem.getCoefficient() == fieldType{-5}) {
				found_d1 = true;
			}
		}
		merge_ok &= found_d0 && found_d1;
	}

	const bool ok = compare_ok && merge_ok;
	std::cout << label << ": ocgraph compare/merge -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_lincomb_of_ocgraphs_addition(const char* label) {
	using DirectionType = GraphDirections<GraphType>;
	using OCBasis = typename OCGraph<GraphType>::Basis;
	using OCLinComb = VectorSpace::LinComb<OCGraph<GraphType>, fieldType>;

	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	DirectionType d0;
	d0.fill(false);
	d0[GraphType::N_HAIR + 0] = true;
	d0[GraphType::N_HAIR + 3] = true;

	DirectionType d1;
	d1.fill(false);
	d1[GraphType::N_HAIR + 1] = true;
	d1[GraphType::N_HAIR + 4] = true;

	typename OCGraph<GraphType>::DirectionComb lincomb_a;
	lincomb_a.append_in_basis_order(d0, fieldType{2});

	typename OCGraph<GraphType>::DirectionComb lincomb_b;
	lincomb_b.append_in_basis_order(d1, fieldType{5});

	OCLinComb a;
	a.append_in_basis_order(OCBasis(OCGraph<GraphType>(source, lincomb_a), fieldType{3}));

	OCLinComb b;
	b.append_in_basis_order(OCBasis(OCGraph<GraphType>(source, lincomb_b), fieldType{-1}));

	OCLinComb sum = a + b;

	bool ok = (sum.size() == 1);
	if (sum.size() == 1) {
		const auto& merged = sum.raw_elements().front();
		ok &= (merged.getCoefficient() == fieldType{1});
		ok &= (merged.getValue().graph == source);
		ok &= (merged.getValue().size() == 2);

		if (merged.getValue().size() == 2) {
			bool found_d0 = false;
			bool found_d1 = false;
			for (const auto& elem : merged.getValue().raw_elements()) {
				if (elem.getValue() == d0 && elem.getCoefficient() == fieldType{6}) {
					found_d0 = true;
				}
				if (elem.getValue() == d1 && elem.getCoefficient() == fieldType{-5}) {
					found_d1 = true;
				}
			}
			ok &= found_d0 && found_d1;
		}
	}

	std::cout << label << ": lincomb<ocgraph> addition -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
bool check_ocgc_wrapper(const char* label) {
	using DirectionType = GraphDirections<GraphType>;
	using OCGraphType = OCGraph<GraphType>;
	using OCGCType = OCGC<GraphType>;

	GraphType source;
	source.setEdge(0, 0, 1);
	source.setEdge(1, 1, 2);
	source.setEdge(2, 0, 2);

	DirectionType d0;
	d0.fill(false);
	d0[GraphType::N_HAIR + 0] = true;
	d0[GraphType::N_HAIR + 2] = true;
	d0[GraphType::N_HAIR + 5] = true;

	DirectionType d1;
	d1.fill(false);
	d1[GraphType::N_HAIR + 1] = true;
	d1[GraphType::N_HAIR + 3] = true;
	d1[GraphType::N_HAIR + 4] = true;

	typename OCGraphType::DirectionComb lincomb_a;
	lincomb_a.append_in_basis_order(d0, fieldType{1});

	typename OCGraphType::DirectionComb lincomb_b;
	lincomb_b.append_in_basis_order(d1, fieldType{2});

	OCGCType a(OCGraphType(source, lincomb_a));
	OCGCType b(OCGraphType(source, lincomb_b));
	OCGCType sum = a + b;

	bool ok = (sum.size() == 1);
	if (sum.size() == 1) {
		const auto& merged = sum.data().raw_elements().front();
		ok &= (merged.getCoefficient() == fieldType{1});
		ok &= (merged.getValue().size() == 2);
	}

	OCGCType delta = a.delta_e();
	OCGraphType expected_delta = a.data().raw_elements().front().getValue().delta_e();

	ok &= (delta.size() == (expected_delta.size() == 0 ? 0 : 1));
	if (delta.size() == 1 && expected_delta.size() != 0) {
		const auto& delta_term = delta.data().raw_elements().front();
		ok &= (delta_term.getCoefficient() == fieldType{1});
		ok &= (delta_term.getValue().graph == expected_delta.graph);
		ok &= (static_cast<const typename OCGraphType::DirectionComb&>(delta_term.getValue())
			== static_cast<const typename OCGraphType::DirectionComb&>(expected_delta));
	}

	std::cout << label << ": ocgc wrapper -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
struct OCGraph_test {
	static bool standardize_and_sort() {
		return check_ocgraph_standardize_and_sort<GraphType>("OCGraph_test.standardize_and_sort");
	}

	static bool filters_non_covering_directions() {
		return check_ocgraph_filters_non_covering_directions<GraphType>("OCGraph_test.filters_non_covering");
	}

	static bool delta_e() {
		return check_ocgraph_delta_e<GraphType>("OCGraph_test.delta_e");
	}

	static bool delta_e_prunes_uncovered_vertices() {
		return check_ocgraph_delta_e_prunes_uncovered_vertices<GraphType>("OCGraph_test.delta_e_prune");
	}

	static bool delta_e_squares_to_zero() {
		return check_ocgraph_delta_e_squares_to_zero<GraphType>("OCGraph_test.delta_e_squared");
	}

	static bool delta_e_isomorphism_equivariance() {
		return check_ocgraph_delta_e_isomorphism_equivariance<GraphType>("OCGraph_test.delta_e_iso");
	}

	static bool homotopy_standardize() {
		return check_ocgraph_homotopy_standardize_minimal_true_case<GraphType>("OCGraph_test.homotopy_standardize");
	}

	static bool compare_and_merge() {
		return check_ocgraph_compare_and_merge<GraphType>("OCGraph_test.compare_merge");
	}

	static bool lincomb_addition() {
		return check_lincomb_of_ocgraphs_addition<GraphType>("OCGraph_test.lincomb_add");
	}

	static bool ocgc() {
		return check_ocgc_wrapper<GraphType>("OCGraph_test.ocgc");
	}

	static bool object() {
		const bool ok = standardize_and_sort()
			&& filters_non_covering_directions()
			&& compare_and_merge()
			&& lincomb_addition()
			&& ocgc()
			&& delta_e()
			&& delta_e_prunes_uncovered_vertices()
			&& delta_e_squares_to_zero()
			&& delta_e_isomorphism_equivariance()
			&& homotopy_standardize();
		std::cout << "OCGraph_test.object: " << (ok ? "ok" : "failed") << '\n';
		return ok;
	}
};
