#include "graph.hpp"
#include <utility>
#include "VectorSpace/BasisElement.hpp"


template <
    Int N_VERTICES,
    Int N_EDGES,
    Int N_OUT_HAIR,
    Int N_IN_HAIR,
    signedInt c,
    signedInt d
>


class GraphStandadizer {
    public:
    using GraphType = Graph<N_VERTICES,N_EDGES,N_OUT_HAIR, N_IN_HAIR,c,d>;

    class CanonBuilder {
        public:
        GraphType G;
        Int n_assignedVertices;
        signedInt sign;

        CanonBuilder(const GraphType& initialGraph, Int n, signedInt s)
            : G(initialGraph), n_assignedVertices(n) {
            n_assignedVertices = n;
            sign = s;

            sign = s*G.directAndSortEdges();
        }

        signedInt sort_edges() {
            sign *= G.directAndSortEdges();
            return sign;
        }

        Int vertex_value(Int v) {
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
                ++ result[a][vertex_value(b)];
                ++result[b][vertex_value(a)];
            }

            return result;
        }

        GraphType vertex_values_graph() {
            GraphType fakeGraph = G;

            for (Int i = 0; i < GraphType::SIZE; ++i) {
                fakeGraph.half_edges[i] = vertex_value(G.half_edges[i]);
            }

            return fakeGraph;
        }

        int compare(CanonBuilder other) {
            return combutils::compareHalfEdges(vertex_values_graph().half_edges, other.vertex_values_graph().half_edges);
        } 

        //we are assigning j <- n_assignedVertices
        CanonBuilder with_assigned_next(Int j) {
            GraphType copied_graph = G;
            signedInt s = copied_graph.swapVertices(j, n_assignedVertices);
            return CanonBuilder(copied_graph, n_assignedVertices +1, sign*s);
        }
    };

    BasisElement<GraphType> standardize(BasisElement<GraphType>& input) {
        return standardize(input.getValue(), input.getCoefficient());
    }    


    BasisElement<GraphType> standardize(GraphType graph, fieldType k) {
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
        
        if(GraphType::SWAP_EDGE_SIGN == -1 && attempts[N_VERTICES%2][0].G.has_double_edge()) {
                return BasisElement<GraphType>(attempts[N_VERTICES%2][0].G, 0); 
        }


        bool containsPlus = false;
        bool containsMinus = false;

        for (bigInt i = 0; i < attempts[N_VERTICES%2].size(); i++) {
                if (attempts[N_VERTICES%2][i].sign > 0) {
                        containsPlus = true;
                } else {
                        containsMinus = true;
                };

                if (containsPlus && containsMinus) {
                        return BasisElement<GraphType>(attempts[N_VERTICES%2][0].G, 0);
                }
        }
        if (containsPlus) {
            return BasisElement<GraphType>(attempts[N_VERTICES%2][0].G, k);
        }
        return BasisElement<GraphType>(attempts[N_VERTICES%2][0].G, -k);
    }
    
    CanonBuilder assignHair(GraphType G) {
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
};