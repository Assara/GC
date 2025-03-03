#include <iostream>
#include "graph.hpp"

using namespace std;


int main() {
    cout << "hello " <<endl;

    array<Int, 9> nums = {0,0,0,0,1,1,1,1,1};

    Graph<2,0,9,0,0,0> g(nums);

    auto splits = g.split_vertex_differential();

    
    cout << "got " << splits.size() << " splits" <<endl;
}