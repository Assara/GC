#include <iostream>
//#include "graph.hpp"
#include "examplegraphs.hpp"

using namespace std;

int main() {
    cout << "hello " <<endl;

    OddGraphdegZero<99> W = wheel_graph<98>();
    OddGraphdegZeroStandadizer<99> G;
    auto res = G.standardize(W, 1.0);

    res.first.print();

    cout << res.second << endl;

}
