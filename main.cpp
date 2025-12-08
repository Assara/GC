#include <iostream>
//#include "graph.hpp"
#include "examplegraphs.hpp"
#include <fstream>
#include "VectorSpace/LinComb.hpp"


using namespace std;


const Int wheelSize = 11;

void tryFindQuadraticRepresentativeForWheel11() {
    OddGraphdegZero<wheelSize + 1> W = wheel_graph<wheelSize>();
    BasisElement<OddGraphdegZero<wheelSize + 1>, fieldType> res = BasisElement<OddGraphdegZero<wheelSize + 1>, fieldType>(W);

    res.getValue().print();
    cout << res.getCoefficient() << endl;

    VectorSpace::LinComb<OddGraphdegZero<wheelSize + 1>, fieldType> dglin(res);

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
    cout << "completed edge differential. size = " << with_edge_diff.size() << endl;
      
    return with_edge_diff
        .try_find_split_primitive();
}


void tryFindFullWheel5ClassByWaterfall() {
    GC loop(loop_graph<9>());

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


    cout << "W5 class:" << endl;

    step3 -> add_edge_differential().print();

}


void tryFindFullWheel7ClassByWaterfall() {
    GC loop(loop_graph<13>());

    cout << "loop: ";
    loop.print();

    auto step1 = push_down_the_waterfall(loop);

    if (!step1.has_value()) {
        cout << "Could not find primitive for step 1!! " << endl;
    }

    
    auto step2 = push_down_the_waterfall(*step1);

    if (!step2.has_value()) {
        cout << "Could not find primitive for step 2!! " << endl;
    }

    
    auto step3 = push_down_the_waterfall(*step2);

    if (!step3) {
        cout << "Could not find primitive for step 3!! " << endl;
    }
    
    auto step4 = push_down_the_waterfall(*step3);

    if (!step4) {
        cout << "Could not find primitive for step 4!! " << endl;
        step3 -> print();
        
        return;
    }
    
    auto step5 = push_down_the_waterfall(*step4);

    if (!step5) {
        cout << "Could not find primitive for step 5!! " << endl;
    }

	OddGCdegZero<8> W7_class =  step5 -> add_edge_differential();

	MetaGraph metaGraph(W7_class.map_split_differential());
	
	
	unordered_set<OddGraphdegZero<8>> important_graphs = metaGraph.component_containing(wheel_graph<7>().canonical_represesentation());
	

	cout<< "found " << important_graphs.size() << " graphs in the component of W7" << endl;


	auto W7_class_filtered = W7_class.filtered(important_graphs);
	

	
	
	auto d_W7_class = W7_class_filtered.delta();
	
	if (d_W7_class.size() != 0) {
			cout << "error! final class not a cocycle" << endl;
			
			std::ofstream out("false_W7.txt");
   			W7_class_filtered.print(out);
   			return;
	}
	
	cout << "found W7 class of size: " << W7_class_filtered.size() << endl;
		
	std::ofstream out("W7.txt");
    
    W7_class_filtered.print(out);
}


void tryFindFullWheel9ClassByWaterfall() {
    GC loop(loop_graph<17>());

    cout << "loop: ";
    loop.print();

	cout << "are we getting step 1?" << endl;

    auto step1 = push_down_the_waterfall(loop);

    if (!step1.has_value()) {
        cout << "Could not find primitive for step 1!! " << endl;
    }
	cout << "step 2" << endl;
    
    auto step2 = push_down_the_waterfall(*step1);

    if (!step2.has_value()) {
        cout << "Could not find primitive for step 2!! " << endl;
    }
	cout << "step 3" << endl;
    
    auto step3 = push_down_the_waterfall(*step2);

    if (!step3.has_value()) {
        cout << "Could not find primitive for step 3!! " << endl;
    }
    
    cout << "trying step 4. step3.size() = " << step3 -> size() << endl;
    
  
    auto step4 = push_down_the_waterfall(*step3);

    if (!step4.has_value()) {
        cout << "Could not find primitive for step 4!! " << endl;
      
        return;
    }
    
    
    
    cout << "step 5" << endl;
    auto step5 = push_down_the_waterfall(*step4);

    if (!step5.has_value()) {
        cout << "Could not find primitive for step 5!! " << endl;
    }
    
    
    
    cout << "step 6" << endl;
    auto step6 = push_down_the_waterfall(*step5);

    if (!step6.has_value()) {
        cout << "Could not find primitive for step 6!! " << endl;
    }
    
    cout << "step 7" << endl;
    
    auto step7 = push_down_the_waterfall(*step6);

    if (!step7.has_value()) {
        cout << "Could not find primitive for step 7!! " << endl;
    }

	OddGCdegZero<10> W9_class =  step7 -> add_edge_differential();

	MetaGraph metaGraph(W9_class.map_split_differential());
	
	
	unordered_set<OddGraphdegZero<10>> important_graphs = metaGraph.component_containing(wheel_graph<9>().canonical_represesentation());
	

	cout<< "found " << important_graphs.size() << " graphs in the component of W9" << endl;


	auto W9_class_filtered = W9_class.filtered(important_graphs);
	

	auto d_W9_class = W9_class_filtered.delta();
	
	if (d_W9_class.size() != 0) {
			cout << "error! final class not a cocycle" << endl;
			
			std::ofstream out("false_W7.txt");
   			d_W9_class.print(out);
   			return;
	}
	
	cout << "found W9 class of size: " << W9_class_filtered.size() << endl;
		
	std::ofstream out("W9.txt");
    
    W9_class_filtered.print(out);
}


int main() {
    tryFindFullWheel9ClassByWaterfall();
    return 0;
}



