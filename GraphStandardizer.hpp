#include <vector>
#include <utility>
#include "VectorSpace/BasisElement.hpp"
#include "GraphIsomorphism.hpp"

template <
Int N_VERTICES,
	Int N_EDGES,
	Int N_OUT_HAIR,
	Int N_IN_HAIR,
	signedInt c,
	signedInt d,
	typename fieldType
	>
	class Graph;

	template <
	Int N_VERTICES,
	Int N_EDGES,
	Int N_OUT_HAIR,
	Int N_IN_HAIR,
	signedInt c,
	signedInt d,
	typename fieldType
	>
	class GraphStandardizer {
		public:
			using GraphType = Graph<N_VERTICES,N_EDGES,N_OUT_HAIR, N_IN_HAIR,c,d, fieldType>;
			using IsomorphismType = GraphIsomorphism<N_VERTICES, N_EDGES>;

			class CanonBuilder {
				public:
					GraphType G;
					Int n_assignedVertices;
					signedInt sign;

					CanonBuilder(const GraphType& initialGraph, Int n, signedInt s)
						: G(initialGraph), n_assignedVertices(n) {
							n_assignedVertices = n;

							sign = s*G.directAndSortEdges();
						}

					signedInt sort_edges() {
						sign *= G.directAndSortEdges();
						return sign;
					}

					Int vertex_value(Int v) const {
						if(v < n_assignedVertices) {
							return v;
						}
						return N_VERTICES;
					}

					array<array<Int, N_VERTICES+1>,N_VERTICES> score_vertices() {
						array<array<Int, N_VERTICES+1>,N_VERTICES> result{};
						for (Int i = G.N_HAIR; i < G.SIZE; i+=2) {
							Int a = G.half_edges[i];
							Int b = G.half_edges[i+1];
							++result[a][vertex_value(b)];
							++result[b][vertex_value(a)];
						}

						return result;
					}

					GraphType vertex_values_graph() const {
						GraphType fakeGraph = G;

						for (Int i = 0; i < GraphType::SIZE; ++i) {
							fakeGraph.half_edges[i] = vertex_value(G.half_edge(i));
						}

						return fakeGraph;
					}

					int compare(const CanonBuilder& other) {
						return combutils::compareHalfEdges(vertex_values_graph().half_edges, other.vertex_values_graph().half_edges);
					} 

					//we are assigning j <- n_assignedVertices
					CanonBuilder with_assigned_next(Int j) {
						GraphType copied_graph = G;
						signedInt s = copied_graph.swapVertices(j, n_assignedVertices);
						return CanonBuilder(copied_graph, n_assignedVertices +1, sign*s);
					}
			};

			class IsoCanonBuilder {
				public:
					GraphType G;
					Int n_assignedVertices;
					signedInt sign;
					IsomorphismType iso;

					IsoCanonBuilder(const GraphType& initialGraph, Int n, signedInt s, IsomorphismType input_iso = IsomorphismType{})
						: G(initialGraph), n_assignedVertices(n), sign(s), iso(input_iso) {
							n_assignedVertices = n;
							sign = s * direct_and_sort_edges_with_isomorphism(G, iso);
						}

					Int vertex_value(Int v) const {
						if(v < n_assignedVertices) {
							return v;
						}
						return N_VERTICES;
					}

					GraphType vertex_values_graph() const {
						GraphType fakeGraph = G;

						for (Int i = 0; i < GraphType::SIZE; ++i) {
							fakeGraph.half_edges[i] = vertex_value(G.half_edge(i));
						}

						return fakeGraph;
					}

					int compare(const IsoCanonBuilder& other) {
						return combutils::compareHalfEdges(vertex_values_graph().half_edges, other.vertex_values_graph().half_edges);
					}

					IsoCanonBuilder with_assigned_next(Int j) {
						GraphType copied_graph = G;
						IsomorphismType copied_iso = iso;
						signedInt s = copied_graph.swapVertices(j, n_assignedVertices);
						copied_iso = copied_iso.compose(vertex_swap_isomorphism(j, n_assignedVertices));
						return IsoCanonBuilder(copied_graph, n_assignedVertices + 1, sign * s, copied_iso);
					}

					static IsomorphismType vertex_swap_isomorphism(Int v, Int w) {
						IsomorphismType iso;
						std::swap(iso.vertex_perm[v], iso.vertex_perm[w]);
						return iso;
					}

					static IsomorphismType edge_flip_isomorphism(Int e) {
						IsomorphismType iso;
						iso.edge_flip[e] = true;
						return iso;
					}

					static IsomorphismType edge_swap_isomorphism(Int i, Int j) {
						IsomorphismType iso;
						std::swap(iso.edge_perm[i], iso.edge_perm[j]);
						return iso;
					}

					static signedInt direct_and_sort_edges_with_isomorphism(GraphType& graph, IsomorphismType& iso) {
						signedInt overall_sign = 1;

						for (Int edge_index = 0; edge_index < N_EDGES; ++edge_index) {
							const Int base = GraphType::N_HAIR + 2 * edge_index;
							if (graph.half_edges[base] > graph.half_edges[base + 1]) {
								overall_sign *= graph.flipEdge(edge_index);
								iso = iso.compose(edge_flip_isomorphism(edge_index));
							}
						}

						for (Int i = 1; i < N_EDGES; ++i) {
							Int j = i;
							while (j > 0 && graph.compareEdge(j - 1, j) > 0) {
								overall_sign *= graph.swapEdges(j - 1, j);
								iso = iso.compose(edge_swap_isomorphism(j - 1, j));
								--j;
							}
						}

						return overall_sign;
					}
			};

			BasisElement<GraphType, fieldType> standardize(BasisElement<GraphType, fieldType>& input) {
				return standardize(input.getValue(), input.getCoefficient());
			}    

			bigInt automorphism_group_size(const GraphType& input_graph) const {
				return minimizing_isomorphisms(input_graph).size();
			}

			std::vector<IsomorphismType> minimizing_isomorphisms(const GraphType& input_graph) const {
				GraphType graph = input_graph;
				auto attempts = canonical_attempts_with_isomorphisms(graph);
				const auto& minimizers = attempts[N_VERTICES % 2];

				std::vector<IsomorphismType> result;
				result.reserve(minimizers.size());
				for (const auto& attempt : minimizers) {
					result.push_back(attempt.iso);
				}
				return result;
			}


			BasisElement<GraphType, fieldType> standardize(GraphType& graph, fieldType k) const {
				if (GraphType::SWAP_EDGE_SIGN == -1 ) {
					k *= graph.directAndSortEdges();
					if (graph.has_double_edge()) {
						return BasisElement<GraphType, fieldType>(graph, 0); 
					}
				}

				CanonBuilder G = assignHair(graph);
				vector<CanonBuilder> attempts[2];

				// TODO reserve appropriate space for attempts

				attempts[G.n_assignedVertices%2].push_back(G);

				for (Int n = G.n_assignedVertices; n < N_VERTICES; n++) {                    

					attempts[(n+1)%2].clear();
					for (CanonBuilder attempt : attempts[n%2]) {
						for (Int l = attempt.n_assignedVertices; l < N_VERTICES; ++l) {
							CanonBuilder next = attempt.with_assigned_next(l);
							if (attempts[(n+1)%2].empty()) {
								attempts[(n+1)%2].push_back(next);
								continue;
							}
							int comparison = next.compare(attempts[(n+1)%2].back());
							if (comparison  > 0 ) {
								continue;
							}
							if (comparison  < 0 ) {
								attempts[(n+1)%2].clear();
							}
							attempts[(n+1)%2].push_back(next);
						}

					}
				}

				//cout << "Aut size = " << attempts[N_VERTICES%2].size() << endl;

				bool containsPlus = false;
				bool containsMinus = false;

				for (bigInt i = 0; i < attempts[N_VERTICES%2].size(); i++) {
					if (attempts[N_VERTICES%2][i].sign > 0) {
						containsPlus = true;
					} else {
						containsMinus = true;
					};

					if (containsPlus && containsMinus) {

						return BasisElement<GraphType, fieldType>(attempts[N_VERTICES%2][0].G, 0);
					}
				}


				BasisElement<GraphType, fieldType> final_form;

				if (containsPlus) {
					final_form = BasisElement<GraphType, fieldType>(attempts[N_VERTICES%2][0].G, k);
				} else {
					final_form= BasisElement<GraphType, fieldType>(attempts[N_VERTICES%2][0].G, -k);
				}

				return final_form;
			}

			CanonBuilder assignHair(GraphType G) const {
				Int n_assigned = 0;
				signedInt sign = 1;
				for (Int i = 0; i < G.N_HAIR ; ++i) {
					if (G.half_edges[i] > n_assigned) {
						sign *= G.swapVertices(G.half_edges[i], n_assigned);
						n_assigned ++;
					} 
				}
				return CanonBuilder(G, n_assigned, sign);
			}

			IsoCanonBuilder assignHair_with_isomorphisms(GraphType G) const {
				Int n_assigned = 0;
				signedInt sign = 1;
				IsomorphismType iso;
				for (Int i = 0; i < G.N_HAIR ; ++i) {
					const Int hair_vertex = G.half_edges[i];
					if (hair_vertex > n_assigned) {
						sign *= G.swapVertices(hair_vertex, n_assigned);
						iso = iso.compose(IsoCanonBuilder::vertex_swap_isomorphism(hair_vertex, n_assigned));
						n_assigned ++;
					}
				}
				return IsoCanonBuilder(G, n_assigned, sign, iso);
			}

			array<vector<IsoCanonBuilder>, 2> canonical_attempts_with_isomorphisms(GraphType& graph) const {
				IsoCanonBuilder G = assignHair_with_isomorphisms(graph);
				array<vector<IsoCanonBuilder>, 2> attempts;

				attempts[G.n_assignedVertices % 2].push_back(G);

				for (Int n = G.n_assignedVertices; n < N_VERTICES; n++) {
					attempts[(n + 1) % 2].clear();
					for (IsoCanonBuilder attempt : attempts[n % 2]) {
						for (Int l = attempt.n_assignedVertices; l < N_VERTICES; ++l) {
							IsoCanonBuilder next = attempt.with_assigned_next(l);
							if (attempts[(n + 1) % 2].empty()) {
								attempts[(n + 1) % 2].push_back(next);
								continue;
							}
							int comparison = next.compare(attempts[(n + 1) % 2].back());
							if (comparison > 0) {
								continue;
							}
							if (comparison < 0) {
								attempts[(n + 1) % 2].clear();
							}
							attempts[(n + 1) % 2].push_back(next);
						}
					}
				}

				return attempts;
			}
	};
