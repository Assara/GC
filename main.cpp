#include <iostream>
//#include "graph.hpp"
#include "examplegraphs.hpp"


#include "VectorSpace/LinComb.hpp"            // your differential graded vector space
#include "VectorSpace/BasisElement.hpp"         // BasisElement
#include "VectorSpace/ValidBasisElement.hpp" // concept enforcement

using namespace std;


const Int wheelSize = 11;

void tryFindQuadraticRepresentativeForWheel11() {
    OddGraphdegZero<wheelSize + 1> W = wheel_graph<wheelSize>();
    BasisElement<OddGraphdegZero<wheelSize + 1>, fieldType> res = BasisElement<OddGraphdegZero<wheelSize + 1>, fieldType>(W);

    res.getValue().print();
    cout << res.getCoefficient() << endl;

    VectorSpace::LinComb<OddGraphdegZero<wheelSize + 1>, fieldType> dglin = VectorSpace::LinComb<OddGraphdegZero<wheelSize + 1>, fieldType>(res);

    dglin.standardize_all();


    GC wheel = GC(W);

    cout << "wheel = ";
    wheel.print();
    int i = 0;
    while (wheel.frontValence() > 4 && i < 4) {
        cout << "attempt" <<  i << endl;
        wheel = wheel.reduce2();

        i++;
        wheel.printFront();        
        cout << "wheel.size() = " << wheel.size() << endl; 
        cout << "grade front =" << wheel.data().raw_elements().front().getValue().custom_filter() << endl;
        cout << "grade back =" << wheel.data().raw_elements().back().getValue().custom_filter() << endl;
    }

    cout << "final representation: " << endl;
    wheel.print();

}

template<Int N>
void tryFindFullWheelClassesByWaterfall() {
    GC loop(loop_graph<4*N+1>());

    auto maybe_W3 = push_down_the_waterfall(loop);

    if (!maybe_W3) {
        cout << "Could not find primitive!! " << endl;

    }

    maybe_W3 ->  add_edge_differential().print();
}


template <
    Int N_VERTICES, Int N_EDGES
> 
 std::optional<GC<N_VERTICES-1, N_EDGES, 0, 0, 0, 1>> push_down_the_waterfall(GC<N_VERTICES, N_EDGES, 0, 0, 0, 1> gamma) {
    auto with_edge_diff = gamma.add_edge_differential();    


    cout << "with edge diff: " ;
    with_edge_diff.print();

    return with_edge_diff
        .try_find_split_primitive();
}


void tryFindFullWheel5ClassByWaterfall() {
    GC loop(loop_graph<9>());

    loop.print();

    auto step1 = push_down_the_waterfall(loop);

    if (!step1) {
        cout << "Could not find primitive for step 1!! " << endl;
    }

    step1 -> print();

    auto test = step1 -> add_edge_differential();


    cout << "TESTLOG ------------------" << endl;

    auto a = step1 -> add_edge_differential();
    auto b = step1 -> delta();


    cout << "original:" << endl;

    step1 -> print();



    cout << "with add_edge_differential " << endl;
    a.print();
    cout << "add edge then split:" << endl;
    
    a.delta().print();



}





int main() {
    tryFindFullWheel5ClassByWaterfall();
    return 0;
}



