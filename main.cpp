#include <iostream>
//#include "graph.hpp"
#include "examplegraphs.hpp"
#include <fstream>
#include "VectorSpace/LinComb.hpp"

using namespace std;


constexpr Int wheelSize = 11;

OddGCdegZero<wheelSize + 1> tryFindQuadraticRepresentativeForWheel() {
	OddGraphdegZero<wheelSize + 1> W = wheel_graph<wheelSize>();
	OddGCdegZero<wheelSize + 1> wheel_class(W);
	auto representative_opt =
		OddGCdegZero<wheelSize + 1>::try_find_quadratic_cont_representative(wheel_class);

	if (!representative_opt.has_value()) {
		cout << "FAILURE: could not find quadratic contraction representative" << std::endl;
		return OddGCdegZero<wheelSize + 1>{};
	}

	cout << "quadratic representative size = " << representative_opt->size() << std::endl;
	return *representative_opt;

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
		.try_find_split_primitive_graded();
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



	auto W5 = step3 -> add_edge_differential();


	cout << "W5 class:" << endl;
	W5.print();


	cout << "dW5 (should be 0)" << endl;


	W5.delta().print();


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
		//step3 -> print();

		return;
	}

	auto step5 = push_down_the_waterfall(*step4);

	if (!step5) {
		cout << "Could not find primitive for step 5!! " << endl;
	}

	OddGCdegZero<8> W7_class =  step5 -> add_edge_differential();

	MetaGraph metaGraph(W7_class.map_split_differential());


	cout << "W7 size before filtering: " << W7_class.size() << endl;

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

void some_truss_graph_exploration () {
	const Int size = 17;

	OddGraphdegZero<size + 1> T = truss_graph<size>();

	T.print();

	cout << "-------------------" << endl;
	GC truss = GC(T);

	auto d_truss = truss.d_contraction();


	vector< OddGraphdegZero<size + 1>> graphs_to_exclude;


	auto canon_T = T.canonical_represesentation();

	graphs_to_exclude.push_back(T.canonical_represesentation());

	auto other_primitive_optional = d_truss.try_find_cont_primitive(graphs_to_exclude); 

	auto other_primitive = *other_primitive_optional;


	bool found_truss_graph = false;
	for (const auto& be: other_primitive.data()) {
		if (be.getValue() == canon_T) {
			cout << "found truss in other primitive. Something is wrong" << endl;
			found_truss_graph= true;
		}
	}

	if (!found_truss_graph) {
		cout << "truss graph is not present in other primitive. All good" << endl;
	}

	truss.scalar_multiply(fieldType{-1}); 
	truss += other_primitive;

	cout << "truss_class size = " << truss.size() << std::endl;

	auto d_truss_class = truss.d_contraction();
	cout << "d_truss_class size = " << d_truss_class.size() << std::endl;


	truss.try_find_cont_primitive();
}

int main() {
	OddGraphdegZero<wheelSize + 1> wheel = wheel_graph<wheelSize>();
	auto representative_opt =
		OddGCdegZero<wheelSize + 1>::try_find_quadratic_cont_representative(
			OddGCdegZero<wheelSize + 1>(wheel)
		);

	if (!representative_opt.has_value()) {
		cout << "FAILURE: could not find quadratic contraction representative for wheel graph" << endl;
		return 1;
	}

	const string output_filename =
		"quadratic_cont_rep_" + std::to_string(wheelSize) + ".txt";
	ofstream out(output_filename);
	if (!out) {
		cerr << "FAILURE: could not open output file " << output_filename << endl;
		return 1;
	}

	representative_opt->print(out);
	cout << "SUCCESS: quadratic representative size = " << representative_opt->size()
		 << ", written to " << output_filename << endl;
	return 0;
}

