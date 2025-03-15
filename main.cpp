#include <iostream>
//#include "graph.hpp"
#include "examplegraphs.hpp"


#include "VectorSpace.hpp"
using namespace std;

int main() {
    cout << "hello " <<endl;

    OddGraphdegZero<17> W = wheel_graph<16>();
    OddGraphdegZeroStandadizer<17> G;

    auto res = G.standardize(W, 1.0);


    res.first.print();

    cout << res.second << endl;
    cout << "bajs" << endl;

}



