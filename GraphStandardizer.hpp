#include "graph.hpp"
#include <utility>



template <
    Int N_VERTICES,
    Int N_EDGES,
    Int N_OUT_HAIR,
    Int N_IN_HAIR,
    signedInt c,
    signedInt d
>


class GraphStandadizer {
    using Graph = Graph<N_VERTICES,N_EDGES,N_OUT_HAIR,c,d>;

    pair<Graph, fieldType> standardize(Graph graph, fieldType k) {
        
    }

    
    return signedInt compare(CanonBuilder graph1, CanonBuilder graph2) {
        //return combutils::compareHalfEdges(graph1.vertex_values.half_edges, graph2.vertex_values.half_edges )
    }







    class CanonBuilder {
        Graph G;
        Int n_assignedVertices;

        signedInt sign;

       
         CanonBuilder(const Graph& initialGraph, Int n, signedInt s)
            : G(initialGraph), n_assignedVertices(n) {
            n_assignedVertices = n;
            sign = s*G.directAndSortEdges();
        }


        
        Int vertex_value(Int v) {
            if(v < assignedVertices) {
                return v;
            }
            return N_VERTICES;
        }

        Array<Int, N_VERTICES> value() {
            Array<Int, G.SIZE> values;
            for (Int i = 0; i< G.SIZE ; ++i) {
                values[i] = vertex_value(G.half_edges[i]);
            }
            return values;
        }

        //we are assigning j <- n_assignedVertices
        GraphStandadizer with_assigned_next(Int j) {
            Graph nextG = G.with_assigned_next(j, n_assignedVertices);
        
        }


    };


    pair<Graph, fieldType> standardize(CanonBuilder graph_info) {

        pair<Int,Int> edge = graph_info.G.getEdge()

    }



};