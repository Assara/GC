#include "GraphStandardizer.hpp"


template


template <Int N_VERTICES, Int N_EDGES, Int N_OUT_HAIR, Int N_IN_HAIR, signedInt c, signedInt d>

class GC {

    using ThisGC =  GC<N_VERTICES,N_EDGES,N_OUT_HAIR, N_IN_HAIR,c,d>;



    using GraphType = Graph<N_VERTICES,N_EDGES,N_OUT_HAIR, N_IN_HAIR,c,d>;

    using SplitGC = GC<N_VERTICES + 1, N_EDGES+1, N_OUT_HAIR, N_IN_HAIR,c,d>;

    private:
    LinComb<GraphType> vec;


    public:
    
    SplitGC delta() {
        for (auto it = vec.begin(); it < vec.end(); it++) {



        }

    }


    static SplitGC delta(&BasisElement<GraphType> G) {
        dad


    } 


};
