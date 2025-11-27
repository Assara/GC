
#include "GC.hpp"
#include "VectorSpace/BasisElement.hpp"



template <Int N_VERTICES, Int N_EDGES, Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d, typename fieldType>
void Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>::
std(BasisElement<ThisGraph, fieldType>& b) {
    GraphStandadizer<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType> s;
    b = s.standardize(b);
}

  