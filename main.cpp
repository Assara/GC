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

    GC loop_with_chord = loop.add_edge_differential();

    loop_with_chord.print();


    auto first_primitive = loop_with_chord.try_find_split_primitive();

    

    auto maybe_W3 = first_primitive-> add_edge_differential();

    cout << "-------------------------------------" << endl;

    maybe_W3.print();

}



int main() {
    tryFindFullWheelClassesByWaterfall<1>();
    return 0;
}



