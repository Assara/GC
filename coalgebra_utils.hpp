#pragma once

#include <algorithm>
#include <unordered_set>
#include <type_traits>
#include <vector>

#include "graph.hpp"
#include "GC.hpp"

namespace coalgebra_utils {

namespace detail {

template <typename GraphType>
bool is_three_four_valent(const GraphType& graph) {
	for (const Int valence : graph.valence_array()) {
		if (valence != 3 && valence != 4) {
			return false;
		}
	}
	return true;
}

template <typename GraphType>
std::vector<Int> three_valent_vertices(const GraphType& graph) {
	std::vector<Int> vertices;
	for (Int v = 0; v < GraphType::N_VERTICES_; ++v) {
		if (graph.valence_array()[v] == 3) {
			vertices.push_back(v);
		}
	}
	return vertices;
}

template <typename GraphType>
Int relabel_outer_vertex_after_insertion(
	Int vertex,
	Int inserted_vertex,
	Int inserted_graph_vertex_count
) {
	if (vertex < inserted_vertex) {
		return vertex;
	}
	if (vertex > inserted_vertex) {
		return static_cast<Int>(vertex + inserted_graph_vertex_count - 1);
	}
	return inserted_vertex;
}

template <typename GraphType>
Int relabel_inner_vertex_after_insertion(Int vertex, Int inserted_vertex) {
	return static_cast<Int>(inserted_vertex + vertex);
}

template <typename GraphType>
std::vector<std::vector<Int>> edge_triangle_adjacency(const GraphType& graph) {
	std::vector<std::vector<Int>> adjacency(GraphType::N_EDGES_);

	for (Int e1 = 0; e1 < GraphType::N_EDGES_; ++e1) {
		const auto [a, b] = graph.getEdge(e1);
		for (Int e2 = e1 + 1; e2 < GraphType::N_EDGES_; ++e2) {
			const auto [c, d] = graph.getEdge(e2);

			Int shared = -1;
			Int u = -1;
			Int v = -1;

			if (a == c) {
				shared = a;
				u = b;
				v = d;
			} else if (a == d) {
				shared = a;
				u = b;
				v = c;
			} else if (b == c) {
				shared = b;
				u = a;
				v = d;
			} else if (b == d) {
				shared = b;
				u = a;
				v = c;
			}

			if (shared < 0 || u == v) {
				continue;
			}

			const Int e3 = graph.find_edge_index(u, v);
			if (e3 == GraphType::N_EDGES_) {
				continue;
			}

			adjacency[e1].push_back(e2);
			adjacency[e1].push_back(e3);
			adjacency[e2].push_back(e1);
			adjacency[e2].push_back(e3);
			adjacency[e3].push_back(e1);
			adjacency[e3].push_back(e2);
		}
	}

	for (auto& nbrs : adjacency) {
		std::sort(nbrs.begin(), nbrs.end());
		nbrs.erase(std::unique(nbrs.begin(), nbrs.end()), nbrs.end());
	}

	return adjacency;
}

} // namespace detail

template <typename OuterGraph, typename InnerGraph>
using GraInsertLowValentGraph = Graph<
	static_cast<Int>(OuterGraph::N_VERTICES_ + InnerGraph::N_VERTICES_ - 1),
	static_cast<Int>(OuterGraph::N_EDGES_ + InnerGraph::N_EDGES_),
	0,
	0,
	OuterGraph::C_,
	OuterGraph::D_,
	typename OuterGraph::Field
>;

template <typename OuterGraph, typename InnerGraph>
using GraInsert4ValentGraph = GraInsertLowValentGraph<OuterGraph, InnerGraph>;

template <typename LeftGC, typename RightGC>
using GraInsertLowValentGC = GC<
	static_cast<Int>(LeftGC::GraphType::N_VERTICES_ + RightGC::GraphType::N_VERTICES_ - 1),
	static_cast<Int>(LeftGC::GraphType::N_EDGES_ + RightGC::GraphType::N_EDGES_),
	0,
	0,
	LeftGC::GraphType::C_,
	LeftGC::GraphType::D_
>;

// Returns the connected triangle banks of a graph, i.e. the connected
// components of the graph on edges in which two edges are adjacent when they
// lie in a common triangle. Edges not contained in any triangle appear as
// singleton banks.
template <typename GraphType>
std::vector<std::vector<Int>> connected_triangle_banks(const GraphType& graph) {
	const auto adjacency = detail::edge_triangle_adjacency(graph);

	std::vector<std::vector<Int>> banks;
	std::vector<bool> visited(GraphType::N_EDGES_, false);

	for (Int start = 0; start < GraphType::N_EDGES_; ++start) {
		if (visited[start]) {
			continue;
		}

		std::vector<Int> bank;
		std::vector<Int> stack{start};
		visited[start] = true;

		while (!stack.empty()) {
			const Int e = stack.back();
			stack.pop_back();
			bank.push_back(e);

			for (const Int f : adjacency[e]) {
				if (!visited[f]) {
					visited[f] = true;
					stack.push_back(f);
				}
			}
		}

		std::sort(bank.begin(), bank.end());
		banks.push_back(std::move(bank));
	}

	std::sort(
		banks.begin(),
		banks.end(),
		[](const std::vector<Int>& lhs, const std::vector<Int>& rhs) {
			if (lhs.size() != rhs.size()) {
				return lhs.size() < rhs.size();
			}
			return lhs < rhs;
		}
	);

	return banks;
}

// Gra-style insertion specialised to the low-valence degree-zero GC_2 setting:
// insert `inner` into a 3- or 4-valent vertex `outer_vertex` of `outer`, and
// sum over all injective attachments of the incident half-edges to distinct
// 3-valent vertices of `inner`. Distinct 3-valent targets are exactly the
// attachments that keep the resulting graph in the 3/4-valent world.
template <typename OuterGraph, typename InnerGraph>
VectorSpace::LinComb<
	GraInsertLowValentGraph<OuterGraph, InnerGraph>,
	typename OuterGraph::Field
> gra_insert_low_valent_vertex_raw(
	const OuterGraph& outer,
	Int outer_vertex,
	const InnerGraph& inner,
	typename OuterGraph::Field coefficient = typename OuterGraph::Field{1}
) {
	static_assert(OuterGraph::N_HAIR == 0 && InnerGraph::N_HAIR == 0,
		"gra_insert_4valent_vertex currently only supports hairless graphs");
	static_assert(OuterGraph::C_ == 0 && OuterGraph::D_ == 1,
		"gra_insert_4valent_vertex is specialised to GC_2 conventions");
	static_assert(InnerGraph::C_ == 0 && InnerGraph::D_ == 1,
		"gra_insert_4valent_vertex is specialised to GC_2 conventions");
	static_assert(std::is_same_v<typename OuterGraph::Field, typename InnerGraph::Field>,
		"outer and inner graph must use the same coefficient field");

	using Field = typename OuterGraph::Field;
	using ResultGraph = GraInsertLowValentGraph<OuterGraph, InnerGraph>;
	using ResultLinComb = VectorSpace::LinComb<ResultGraph, Field>;

	ResultLinComb result;
	if (coefficient == Field{0}) {
		return result;
	}

	if (!detail::is_three_four_valent(outer) || !detail::is_three_four_valent(inner)) {
		return result;
	}

	const auto outer_valences = outer.valence_array();
	if (outer_vertex >= OuterGraph::N_VERTICES_) {
		return result;
	}

	const Int outer_vertex_valence = outer_valences[outer_vertex];
	if (outer_vertex_valence != 3 && outer_vertex_valence != 4) {
		return result;
	}

	auto incident_half_edges = outer.adjacent(outer_vertex);
	if (incident_half_edges.size() != static_cast<std::size_t>(outer_vertex_valence)) {
		return result;
	}
	std::sort(incident_half_edges.begin(), incident_half_edges.end());

	const auto inner_three_valent_vertices = detail::three_valent_vertices(inner);
	if (inner_three_valent_vertices.size() < static_cast<std::size_t>(outer_vertex_valence)) {
		return result;
	}

	bigInt reserve_count = 1;
	for (Int i = 0; i < outer_vertex_valence; ++i) {
		reserve_count *= (inner_three_valent_vertices.size() - static_cast<std::size_t>(i));
	}
	result.reserve(reserve_count);

	std::vector<Int> chosen_targets;
	chosen_targets.reserve(outer_vertex_valence);
	std::vector<bool> used(inner_three_valent_vertices.size(), false);

	const auto append_insertion = [&](const std::vector<Int>& chosen_vertices) {
		std::array<Int, ResultGraph::SIZE> attachment_targets{};
		std::array<bool, ResultGraph::SIZE> is_attached_half_edge{};

		for (std::size_t i = 0; i < chosen_vertices.size(); ++i) {
			const Int half_edge = incident_half_edges[i];
			attachment_targets[half_edge] = chosen_vertices[i];
			is_attached_half_edge[half_edge] = true;
		}

		ResultGraph composed;

		for (Int e = 0; e < OuterGraph::N_EDGES_; ++e) {
			const Int left_half_edge = OuterGraph::N_HAIR + 2 * e;
			const Int right_half_edge = left_half_edge + 1;
			const auto [u, v] = outer.getEdge(e);

			const Int mapped_u = is_attached_half_edge[left_half_edge]
				? attachment_targets[left_half_edge]
				: detail::relabel_outer_vertex_after_insertion<OuterGraph>(
					u,
					outer_vertex,
					InnerGraph::N_VERTICES_
				);
			const Int mapped_v = is_attached_half_edge[right_half_edge]
				? attachment_targets[right_half_edge]
				: detail::relabel_outer_vertex_after_insertion<OuterGraph>(
					v,
					outer_vertex,
					InnerGraph::N_VERTICES_
				);

			composed.setEdge(e, mapped_u, mapped_v);
		}

		for (Int e = 0; e < InnerGraph::N_EDGES_; ++e) {
			const auto [u, v] = inner.getEdge(e);
			composed.setEdge(
				static_cast<Int>(OuterGraph::N_EDGES_ + e),
				detail::relabel_inner_vertex_after_insertion<InnerGraph>(u, outer_vertex),
				detail::relabel_inner_vertex_after_insertion<InnerGraph>(v, outer_vertex)
			);
		}

		result.append_in_basis_order(composed, coefficient);
	};

	const auto backtrack = [&](auto&& self) -> void {
		if (chosen_targets.size() == static_cast<std::size_t>(outer_vertex_valence)) {
			append_insertion(chosen_targets);
			return;
		}

		for (std::size_t i = 0; i < inner_three_valent_vertices.size(); ++i) {
			if (used[i]) {
				continue;
			}

			used[i] = true;
			chosen_targets.push_back(
				detail::relabel_inner_vertex_after_insertion<InnerGraph>(
					inner_three_valent_vertices[i],
					outer_vertex
				)
			);
			self(self);
			chosen_targets.pop_back();
			used[i] = false;
		}
	};

	backtrack(backtrack);

	return result;
}

template <typename OuterGraph, typename InnerGraph>
VectorSpace::LinComb<
	GraInsertLowValentGraph<OuterGraph, InnerGraph>,
	typename OuterGraph::Field
> gra_insert_low_valent_vertex(
	const OuterGraph& outer,
	Int outer_vertex,
	const InnerGraph& inner,
	typename OuterGraph::Field coefficient = typename OuterGraph::Field{1}
) {
	using Field = typename OuterGraph::Field;
	using ResultGraph = GraInsertLowValentGraph<OuterGraph, InnerGraph>;
	using ResultLinComb = VectorSpace::LinComb<ResultGraph, Field>;

	ResultLinComb result = gra_insert_low_valent_vertex_raw(
		outer,
		outer_vertex,
		inner,
		coefficient
	);

	result.standardize_and_sort();
	return result;
}

template <typename OuterGraph, typename InnerGraph>
std::unordered_set<GraInsertLowValentGraph<OuterGraph, InnerGraph>>
gra_insertion_support(
	const OuterGraph& outer,
	const InnerGraph& inner,
	bool only_3valent_vertices = false
) {
	using ResultGraph = GraInsertLowValentGraph<OuterGraph, InnerGraph>;

	std::unordered_set<ResultGraph> support;

	const auto valences = outer.valence_array();
	for (Int vertex = 0; vertex < OuterGraph::N_VERTICES_; ++vertex) {
		if (valences[vertex] != 3 && valences[vertex] != 4) {
			continue;
		}
		if (only_3valent_vertices && valences[vertex] != 3) {
			continue;
		}

		auto raw = gra_insert_low_valent_vertex_raw(outer, vertex, inner, typename OuterGraph::Field{1});
		for (const auto& be : raw) {
			typename ResultGraph::Basis canon_basis(be.getValue(), typename OuterGraph::Field{1});
			const auto canon = ResultGraph::canonized(canon_basis);
			if (canon.getCoefficient() != typename OuterGraph::Field{0}) {
				support.insert(canon.getValue());
			}
		}
	}

	return support;
}

template <typename OuterGraph, typename InnerGraph>
VectorSpace::LinComb<
	GraInsertLowValentGraph<OuterGraph, InnerGraph>,
	typename OuterGraph::Field
> gra_pre_lie(
	const OuterGraph& outer,
	const InnerGraph& inner,
	typename OuterGraph::Field coefficient = typename OuterGraph::Field{1}
) {
	using Field = typename OuterGraph::Field;
	using ResultGraph = GraInsertLowValentGraph<OuterGraph, InnerGraph>;
	using ResultLinComb = VectorSpace::LinComb<ResultGraph, Field>;

	ResultLinComb result;
	if (coefficient == Field{0}) {
		return result;
	}

	const auto valences = outer.valence_array();
	for (Int vertex = 0; vertex < OuterGraph::N_VERTICES_; ++vertex) {
		if (valences[vertex] != 3 && valences[vertex] != 4) {
			continue;
		}
		result += gra_insert_low_valent_vertex(outer, vertex, inner, coefficient);
	}

	result.standardize_and_sort();
	return result;
}

template <typename LeftGraph, typename RightGraph>
VectorSpace::LinComb<
	GraInsertLowValentGraph<LeftGraph, RightGraph>,
	typename LeftGraph::Field
> gra_lie(
	const LeftGraph& left,
	const RightGraph& right,
	typename LeftGraph::Field coefficient = typename LeftGraph::Field{1}
) {
	using Field = typename LeftGraph::Field;

	auto result = gra_pre_lie(left, right, coefficient);
	result = result.add_scaled(gra_pre_lie(right, left, coefficient), Field{-1});
	result.standardize_all();
	result.sort_elements();
	return result;
}

template <typename OuterGraph, typename InnerGraph>
VectorSpace::LinComb<
	GraInsert4ValentGraph<OuterGraph, InnerGraph>,
	typename OuterGraph::Field
> gra_insert_4valent_vertex(
	const OuterGraph& outer,
	Int outer_vertex,
	const InnerGraph& inner,
	typename OuterGraph::Field coefficient = typename OuterGraph::Field{1}
) {
	const auto valences = outer.valence_array();
	if (outer_vertex >= OuterGraph::N_VERTICES_ || valences[outer_vertex] != 4) {
		return {};
	}
	return gra_insert_low_valent_vertex(outer, outer_vertex, inner, coefficient);
}

template <typename OuterGraph, typename InnerGraph>
VectorSpace::LinComb<
	GraInsertLowValentGraph<OuterGraph, InnerGraph>,
	typename OuterGraph::Field
> gra_insert_3valent_vertex(
	const OuterGraph& outer,
	Int outer_vertex,
	const InnerGraph& inner,
	typename OuterGraph::Field coefficient = typename OuterGraph::Field{1}
) {
	const auto valences = outer.valence_array();
	if (outer_vertex >= OuterGraph::N_VERTICES_ || valences[outer_vertex] != 3) {
		return {};
	}
	return gra_insert_low_valent_vertex(outer, outer_vertex, inner, coefficient);
}

template <typename OuterGraph, typename InnerGraph>
VectorSpace::LinComb<
	GraInsertLowValentGraph<OuterGraph, InnerGraph>,
	typename OuterGraph::Field
> gra_pre_lie_3valent_only(
	const OuterGraph& outer,
	const InnerGraph& inner,
	typename OuterGraph::Field coefficient = typename OuterGraph::Field{1}
) {
	using Field = typename OuterGraph::Field;
	using ResultGraph = GraInsertLowValentGraph<OuterGraph, InnerGraph>;
	using ResultLinComb = VectorSpace::LinComb<ResultGraph, Field>;

	ResultLinComb result;
	if (coefficient == Field{0}) {
		return result;
	}

	const auto valences = outer.valence_array();
	for (Int vertex = 0; vertex < OuterGraph::N_VERTICES_; ++vertex) {
		if (valences[vertex] != 3) {
			continue;
		}
		result += gra_insert_3valent_vertex(outer, vertex, inner, coefficient);
	}

	result.standardize_and_sort();
	return result;
}

template <typename LeftGraph, typename RightGraph>
VectorSpace::LinComb<
	GraInsertLowValentGraph<LeftGraph, RightGraph>,
	typename LeftGraph::Field
> gra_lie_3valent_only(
	const LeftGraph& left,
	const RightGraph& right,
	typename LeftGraph::Field coefficient = typename LeftGraph::Field{1}
) {
	using Field = typename LeftGraph::Field;

	auto result = gra_pre_lie_3valent_only(left, right, coefficient);
	result = result.add_scaled(gra_pre_lie_3valent_only(right, left, coefficient), Field{-1});
	result.standardize_all();
	result.sort_elements();
	return result;
}

template <typename LeftGC, typename RightGC>
GraInsertLowValentGC<LeftGC, RightGC> gra_pre_lie(
	const LeftGC& left,
	const RightGC& right,
	fieldType coefficient = fieldType{1}
) {
	using ResultGC = GraInsertLowValentGC<LeftGC, RightGC>;
	using ResultGraph = typename ResultGC::GraphType;
	using ResultLinComb = VectorSpace::LinComb<ResultGraph, fieldType>;

	ResultLinComb result;
	for (const auto& left_be : left.data()) {
		for (const auto& right_be : right.data()) {
			const fieldType pair_coefficient =
				coefficient * left_be.getCoefficient() * right_be.getCoefficient();
			result += gra_pre_lie(left_be.getValue(), right_be.getValue(), pair_coefficient);
		}
	}

	result.standardize_and_sort();
	return ResultGC(std::move(result.raw_elements_nonconst()));
}

template <typename LeftGC, typename RightGC>
GraInsertLowValentGC<LeftGC, RightGC> gra_lie(
	const LeftGC& left,
	const RightGC& right,
	fieldType coefficient = fieldType{1}
) {
	using ResultGC = GraInsertLowValentGC<LeftGC, RightGC>;
	using ResultLinComb = typename ResultGC::L;

	ResultLinComb result = gra_pre_lie(left, right, coefficient).data();
	result = result.add_scaled(gra_pre_lie(right, left, coefficient).data(), fieldType{-1});
	result.standardize_all();
	result.sort_elements();
	return ResultGC(std::move(result.raw_elements_nonconst()));
}

template <typename LeftGC, typename RightGC>
GraInsertLowValentGC<LeftGC, RightGC> gra_pre_lie_3valent_only(
	const LeftGC& left,
	const RightGC& right,
	fieldType coefficient = fieldType{1}
) {
	using ResultGC = GraInsertLowValentGC<LeftGC, RightGC>;
	using ResultGraph = typename ResultGC::GraphType;
	using ResultLinComb = VectorSpace::LinComb<ResultGraph, fieldType>;

	ResultLinComb result;
	for (const auto& left_be : left.data()) {
		for (const auto& right_be : right.data()) {
			const fieldType pair_coefficient =
				coefficient * left_be.getCoefficient() * right_be.getCoefficient();
			result += gra_pre_lie_3valent_only(left_be.getValue(), right_be.getValue(), pair_coefficient);
		}
	}

	result.standardize_and_sort();
	return ResultGC(std::move(result.raw_elements_nonconst()));
}

template <typename LeftGC, typename RightGC>
GraInsertLowValentGC<LeftGC, RightGC> gra_lie_3valent_only(
	const LeftGC& left,
	const RightGC& right,
	fieldType coefficient = fieldType{1}
) {
	using ResultGC = GraInsertLowValentGC<LeftGC, RightGC>;
	using ResultLinComb = typename ResultGC::L;

	ResultLinComb result = gra_pre_lie_3valent_only(left, right, coefficient).data();
	result = result.add_scaled(gra_pre_lie_3valent_only(right, left, coefficient).data(), fieldType{-1});
	result.standardize_all();
	result.sort_elements();
	return ResultGC(std::move(result.raw_elements_nonconst()));
}

} // namespace coalgebra_utils
