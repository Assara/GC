#include <iostream>
//#include "graph.hpp"
#include "examplegraphs.hpp"


#include "VectorSpace/LinComb.hpp"            // your differential graded vector space
#include "VectorSpace/BasisElement.hpp"         // BasisElement
#include "VectorSpace/ValidForDifferential.hpp" // concept enforcement


using namespace std;


const Int wheelSize = 11;

int main() {
    cout << "hello " <<endl;

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

    return 0;
}



