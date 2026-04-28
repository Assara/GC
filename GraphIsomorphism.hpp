#pragma once

#include <array>
#include <utility>
#include <vector>

#include "GraphDirections.hpp"
#include "permutation.hpp"

template <
Int N_VERTICES,
Int N_EDGES,
Int N_OUT_HAIR,
Int N_IN_HAIR,
signedInt c,
signedInt d,
typename FieldType
>
class Graph;

template <Int N_VERTICES, Int N_EDGES>
class GraphIsomorphism {
	public:
		std::array<Int, N_VERTICES> vertex_perm{};
		std::array<Int, N_EDGES> edge_perm{};
		std::array<bool, N_EDGES> edge_flip{};

		GraphIsomorphism() {
			for (Int i = 0; i < N_VERTICES; ++i) {
				vertex_perm[i] = i;
			}
			for (Int i = 0; i < N_EDGES; ++i) {
				edge_perm[i] = i;
				edge_flip[i] = false;
			}
		}

		template <Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename FieldType>
		using GraphType = Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>;

		template <Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename FieldType>
		GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>
		permute(const GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>& source) const {
			GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType> target;

			for (Int i = 0; i < N_OUT_HAIR + N_IN_HAIR; ++i) {
				target.half_edges[i] = vertex_perm[source.half_edges[i]];
			}

			for (Int source_edge = 0; source_edge < N_EDGES; ++source_edge) {
				auto [u, v] = source.getEdge(source_edge);
				u = vertex_perm[u];
				v = vertex_perm[v];

				if (edge_flip[source_edge]) {
					std::swap(u, v);
				}

				target.setEdge(edge_perm[source_edge], u, v);
			}

			return target;
		}

		template <Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename FieldType>
		void permute_inplace(GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>& graph) const {
			graph = permute(graph);
		}

		template <Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename FieldType>
		GraphDirections<GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>>
		permute(const GraphDirections<GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>>& source) const {
			using ThisGraphType = GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>;
			GraphDirections<ThisGraphType> target;

			for (Int i = 0; i < ThisGraphType::N_HAIR; ++i) {
				target[i] = source[i];
			}

			for (Int source_edge = 0; source_edge < N_EDGES; ++source_edge) {
				const Int source_base = ThisGraphType::N_HAIR + 2 * source_edge;
				const Int target_base = ThisGraphType::N_HAIR + 2 * edge_perm[source_edge];

				if (edge_flip[source_edge]) {
					target[target_base] = source[source_base + 1];
					target[target_base + 1] = source[source_base];
				} else {
					target[target_base] = source[source_base];
					target[target_base + 1] = source[source_base + 1];
				}
			}

			return target;
		}

		template <Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename FieldType>
		signedInt direction_sign(const GraphDirections<GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>>& source) const {
			using ThisGraphType = GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>;
			std::vector<Int> mapped_true_indices;
			mapped_true_indices.reserve(ThisGraphType::SIZE);

			for (Int i = 0; i < ThisGraphType::N_HAIR; ++i) {
				if (source[i]) {
					mapped_true_indices.push_back(i);
				}
			}

			for (Int source_edge = 0; source_edge < N_EDGES; ++source_edge) {
				const Int source_base = ThisGraphType::N_HAIR + 2 * source_edge;
				const Int target_base = ThisGraphType::N_HAIR + 2 * edge_perm[source_edge];

				if (source[source_base]) {
					mapped_true_indices.push_back(edge_flip[source_edge] ? target_base + 1 : target_base);
				}
				if (source[source_base + 1]) {
					mapped_true_indices.push_back(edge_flip[source_edge] ? target_base : target_base + 1);
				}
			}

			signedInt sign = 1;
			for (bigInt i = 0; i < mapped_true_indices.size(); ++i) {
				for (bigInt j = 0; j < i; ++j) {
					if (mapped_true_indices[j] > mapped_true_indices[i]) {
						sign = -sign;
					}
				}
			}

			return sign;
		}

		template <Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename FieldType>
		void permute_inplace(GraphDirections<GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>>& directions) const {
			directions = permute(directions);
		}

		GraphIsomorphism compose(const GraphIsomorphism& after) const {
			GraphIsomorphism result;

			for (Int v = 0; v < N_VERTICES; ++v) {
				result.vertex_perm[v] = after.vertex_perm[vertex_perm[v]];
			}

			for (Int e = 0; e < N_EDGES; ++e) {
				const Int intermediate_edge = edge_perm[e];
				result.edge_perm[e] = after.edge_perm[intermediate_edge];
				result.edge_flip[e] = edge_flip[e] ^ after.edge_flip[intermediate_edge];
			}

			return result;
		}

		GraphIsomorphism inverse() const {
			GraphIsomorphism result;

			for (Int v = 0; v < N_VERTICES; ++v) {
				result.vertex_perm[vertex_perm[v]] = v;
			}

			for (Int e = 0; e < N_EDGES; ++e) {
				const Int target_edge = edge_perm[e];
				result.edge_perm[target_edge] = e;
				result.edge_flip[target_edge] = edge_flip[e];
			}

			return result;
		}

		Permutation<N_VERTICES> vertex_permutation() const {
			return Permutation<N_VERTICES>(vertex_perm);
		}
};
