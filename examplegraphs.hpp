#include "graph.cpp"

using namespace std;

namespace graphaliases {

    // Alias for a graph with odd degree zero condition.
    template <Int N>
    using OddGraphdegZero = Graph<N, 2 * N - 2, 0, 0, 0, 1>;

    // Alias for the GraphStandadizer corresponding to oddGraphdegZero.
    template <Int N>
    using OddGraphdegZeroStandadizer = GraphStandadizer<N, 2 * N - 2, 0, 0, 0, 1>;

}

using namespace graphaliases;

template<Int N> 
OddGraphdegZero<N+1> wheel_graph() {
    array<Int, 4*N> arr{};

    OddGraphdegZero<N+1>  W(arr);

    //add spokes
    for (Int i = 0; i<N; ++i) {
        W.setEdge(i, 0, i+1);
    }
    //add rim
    for (Int i = 0; i <N-1; ++i) {
        W.setEdge(N+i, i+1, i+2); 

    }
    W.setEdge(2*N-1, 1, N); 

    return W;

}