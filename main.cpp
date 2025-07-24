#include <iostream>
//#include "graph.hpp"
#include "examplegraphs.hpp"


#include "VectorSpace/LinComb.hpp"            // your differential graded vector space
#include "VectorSpace/BasisElement.hpp"         // BasisElement
#include "VectorSpace/ValidForDifferential.hpp" // concept enforcement


using namespace std;

int main() {
    cout << "hello " <<endl;

    OddGraphdegZero<6> W = wheel_graph<5>();
   

    BasisElement<OddGraphdegZero<6>> res = BasisElement(W);


    res.getValue().print();
    cout << res.getCoefficient() << endl;

    VectorSpace::LinComb<OddGraphdegZero<6>> dglin = VectorSpace::LinComb<OddGraphdegZero<6>>(res);

    dglin.standardize_all();


     GC wheel = GC(W);

    GC dWheel = wheel.delta();

    dWheel.print();



    GC ddWheel = dWheel.delta();


    cout << "delta squared : "<< endl;
    ddWheel.print();
    return 0;

}




