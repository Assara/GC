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
		GraphIsomorphism() {
			for (Int i = 0; i < N_VERTICES; ++i) {
				vertex_perm_[i] = i;
			}
			for (Int i = 0; i < N_EDGES; ++i) {
				edge_perm_[i] = i;
				edge_flip_[i] = false;
			}
		}

		template <Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename FieldType>
		using GraphType = Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>;

		template <Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename FieldType>
		GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>
		permute(const GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>& source) const {
			GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType> target;

			for (Int i = 0; i < N_OUT_HAIR + N_IN_HAIR; ++i) {
				target.half_edges[i] = vertex_perm_[source.half_edges[i]];
			}

			for (Int source_edge = 0; source_edge < N_EDGES; ++source_edge) {
				auto [u, v] = source.getEdge(source_edge);
				u = vertex_perm_[u];
				v = vertex_perm_[v];

				if (edge_flip_[source_edge]) {
					std::swap(u, v);
				}

				target.setEdge(edge_perm_[source_edge], u, v);
			}

			return target;
		}

		template <Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename FieldType>
		void permute_inplace(GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>& graph) const {
			graph = permute(graph);
		}

		template <Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename FieldType>
		signedInt graph_sign(const GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>&) const {
			return graph_perm_sign_;
		}

		signedInt edge_permutation_sign() const {
			return edge_perm_sign_;
		}

		signedInt vertex_permutation_sign() const {
			return vertex_perm_sign_;
		}

		signedInt graph_permutation_sign() const {
			return graph_perm_sign_;
		}

		const std::array<Int, N_VERTICES>& vertex_permutation_data() const {
			return vertex_perm_;
		}

		std::array<Int, N_VERTICES>& vertex_permutation_data() {
			return vertex_perm_;
		}

		const std::array<Int, N_EDGES>& edge_permutation_data() const {
			return edge_perm_;
		}

		std::array<Int, N_EDGES>& edge_permutation_data() {
			return edge_perm_;
		}

		const std::array<bool, N_EDGES>& edge_flip_data() const {
			return edge_flip_;
		}

		std::array<bool, N_EDGES>& edge_flip_data() {
			return edge_flip_;
		}

		template <Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename FieldType>
		void compute_signs() {
			vertex_perm_sign_ = Permutation<N_VERTICES>(vertex_perm_).sign();
			edge_perm_sign_ = Permutation<N_EDGES>(edge_perm_).sign();

			signedInt edge_flip_sign = 1;
			for (Int edge_index = 0; edge_index < N_EDGES; ++edge_index) {
				if (edge_flip_[edge_index]) {
					edge_flip_sign = -edge_flip_sign;
				}
			}

			graph_perm_sign_ = 1;
			if constexpr (GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>::SWAP_VERTICES_SIGN == -1) {
				graph_perm_sign_ *= vertex_perm_sign_;
			}
			if constexpr (GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>::SWAP_EDGE_SIGN == -1) {
				graph_perm_sign_ *= edge_perm_sign_;
			}
			if constexpr (GraphType<N_OUT_HAIR, N_IN_HAIR, c, d, FieldType>::FLIP_EDGE_SIGN == -1) {
				graph_perm_sign_ *= edge_flip_sign;
			}
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
				const Int target_base = ThisGraphType::N_HAIR + 2 * edge_perm_[source_edge];

				if (edge_flip_[source_edge]) {
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
				const Int target_base = ThisGraphType::N_HAIR + 2 * edge_perm_[source_edge];

				if (source[source_base]) {
					mapped_true_indices.push_back(edge_flip_[source_edge] ? target_base + 1 : target_base);
				}
				if (source[source_base + 1]) {
					mapped_true_indices.push_back(edge_flip_[source_edge] ? target_base : target_base + 1);
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
				result.vertex_perm_[v] = after.vertex_perm_[vertex_perm_[v]];
			}

			for (Int e = 0; e < N_EDGES; ++e) {
				const Int intermediate_edge = edge_perm_[e];
				result.edge_perm_[e] = after.edge_perm_[intermediate_edge];
				result.edge_flip_[e] = edge_flip_[e] ^ after.edge_flip_[intermediate_edge];
			}

			result.vertex_perm_sign_ = vertex_perm_sign_ * after.vertex_perm_sign_;
			result.edge_perm_sign_ = edge_perm_sign_ * after.edge_perm_sign_;
			result.graph_perm_sign_ = graph_perm_sign_ * after.graph_perm_sign_;

			return result;
		}

		GraphIsomorphism inverse() const {
			GraphIsomorphism result;

			for (Int v = 0; v < N_VERTICES; ++v) {
				result.vertex_perm_[vertex_perm_[v]] = v;
			}

			for (Int e = 0; e < N_EDGES; ++e) {
				const Int target_edge = edge_perm_[e];
				result.edge_perm_[target_edge] = e;
				result.edge_flip_[target_edge] = edge_flip_[e];
			}

			result.vertex_perm_sign_ = vertex_perm_sign_;
			result.edge_perm_sign_ = edge_perm_sign_;
			result.graph_perm_sign_ = graph_perm_sign_;

			return result;
		}

		Permutation<N_VERTICES> vertex_permutation() const {
			return Permutation<N_VERTICES>(vertex_perm_);
		}

	private:
		std::array<Int, N_VERTICES> vertex_perm_{};
		std::array<Int, N_EDGES> edge_perm_{};
		std::array<bool, N_EDGES> edge_flip_{};
		signedInt vertex_perm_sign_ = 1;
		signedInt edge_perm_sign_ = 1;
		signedInt graph_perm_sign_ = 1;
};
