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

    cout << "using push_down_the_waterfall: " ;
    auto with_edge_diff = gamma.add_edge_differential();    
    

    auto primitive = with_edge_diff
        .try_find_split_primitive();

    cout << "finished edge diff:" << endl;
    return primitive;
}


void tryFindFullWheel5ClassByWaterfall() {
    GC loop(loop_graph<9>());

    cout << "loop: ";
    loop.print();

    auto step1 = push_down_the_waterfall(loop);

    if (!step1) {
        cout << "Could not find primitive for step 1!! " << endl;
    }

    cout << "split then add: " << endl;
    step1 -> delta().add_edge_differential().print();


    cout << "add then split: " << endl;
    step1 -> add_edge_differential ().delta().print();


        cout << "delta ^2: " << endl;
    step1 -> delta().delta().print();

    /*
    auto step2 = push_down_the_waterfall(*step1);

    if (!step2) {
        cout << "Could not find primitive for step 2!! " << endl;
    }

    auto step3 = push_down_the_waterfall(*step2);

    if (!step3) {
        cout << "Could not find primitive for step 3!! " << endl;
    }

    auto full_W5_class = step3 -> add_edge_differential();

    cout << "full W5 class" << endl;
    full_W5_class.print();

    */
}

/*

void tryFindFullWheel7ClassByWaterfall() {
    GC loop(loop_graph<13>());

    cout << "loop: ";
    loop.print();

    auto step1 = push_down_the_waterfall(loop);

    if (!step1) {
        cout << "Could not find primitive for step 1!! " << endl;
    }

    auto step2 = push_down_the_waterfall(*step1);

    if (!step2) {
        cout << "Could not find primitive for step 2!! " << endl;
    }

    auto step3 = push_down_the_waterfall(*step2);

    if (!step3) {
        cout << "Could not find primitive for step 3!! " << endl;
    }

    
    auto step4 = push_down_the_waterfall(*step3);

    if (!step4) {
        cout << "Could not find primitive for step 4!! " << endl;
    }

    auto step5 = push_down_the_waterfall(*step4);

    if (!step5) {
        cout << "Could not find primitive for step 5!! " << endl;
    }


    auto full_W7_class = step5 -> add_edge_differential();

    cout << "full W7 class" << endl;
    full_W7_class.print();
}
*/


int main() {
    tryFindFullWheel5ClassByWaterfall();
    return 0;
}



