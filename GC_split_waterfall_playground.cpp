#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>

#include "examplegraphs.hpp"
#include "MetaGraph.hpp"

using namespace std;

template <
	Int N_VERTICES,
	Int N_EDGES
>
std::optional<GC<N_VERTICES - 1, N_EDGES, 0, 0, 0, 1>>
push_down_the_waterfall(GC<N_VERTICES, N_EDGES, 0, 0, 0, 1> gamma) {
	cout << "using push_down_the_waterfall: ";
	auto with_edge_diff = gamma.add_edge_differential();
	return with_edge_diff.try_find_split_primitive();
}

template <
	Int N_VERTICES,
	Int N_EDGES
>
std::optional<GC<N_VERTICES - 1, N_EDGES, 0, 0, 0, 1>>
push_down_with_generated_split_contract_support(
	GC<N_VERTICES, N_EDGES, 0, 0, 0, 1> gamma,
	int generated_rounds
) {
	cout << "using generated split-contract support: ";
	auto with_edge_diff = gamma.add_edge_differential();
	return with_edge_diff.try_find_split_primitive_generated(generated_rounds);
}

template <
	Int N_VERTICES,
	Int N_EDGES
>
std::optional<GC<N_VERTICES - 1, N_EDGES, 0, 0, 0, 1>>
push_down_by_selected_method(
	GC<N_VERTICES, N_EDGES, 0, 0, 0, 1> gamma,
	int generated_rounds
) {
	if (generated_rounds >= 0) {
		return push_down_with_generated_split_contract_support(gamma, generated_rounds);
	}
	return push_down_the_waterfall(gamma);
}

template <typename GCType>
void print_missing_step(const char* label) {
	cout << "Could not find primitive for " << label << "!!" << endl;
}

template <Int WheelN>
std::optional<OddGCdegZero<WheelN + 1>>
filter_to_wheel_component(OddGCdegZero<WheelN + 1> wheel_class) {
	MetaGraph meta_graph(wheel_class.map_split_differential());
	auto important_graphs =
		meta_graph.component_containing(wheel_graph<WheelN>().canonical_represesentation());

	cout << "found " << important_graphs.size()
	     << " graphs in the component of W" << +WheelN << endl;

	auto filtered = wheel_class.filtered(important_graphs);
	auto d_filtered = filtered.delta();
	if (d_filtered.size() != 0) {
		cout << "error! final class not a cocycle" << endl;
		return std::nullopt;
	}

	return filtered;
}

std::optional<OddGCdegZero<4>> tryFindFullWheel3Class(int generated_rounds) {
	GC loop(loop_graph<5>());

	cout << "loop: ";
	loop.print();

	auto step1 = push_down_by_selected_method(loop, generated_rounds);
	if (!step1) {
		print_missing_step<decltype(loop)>("step 1");
		return std::nullopt;
	}

	cout << "W3 class:" << endl;
	auto W3_class = step1->add_edge_differential();
	W3_class.print();
	return W3_class;
}

std::optional<OddGCdegZero<6>> tryFindFullWheel5Class(int generated_rounds) {
	GC loop(loop_graph<9>());

	cout << "loop: ";
	loop.print();

	auto step1 = push_down_by_selected_method(loop, generated_rounds);
	if (!step1) {
		print_missing_step<decltype(loop)>("step 1");
		return std::nullopt;
	}

	auto step2 = push_down_by_selected_method(*step1, generated_rounds);
	if (!step2) {
		print_missing_step<decltype(*step1)>("step 2");
		return std::nullopt;
	}

	auto step3 = push_down_by_selected_method(*step2, generated_rounds);
	if (!step3) {
		print_missing_step<decltype(*step2)>("step 3");
		return std::nullopt;
	}

	cout << "W5 class:" << endl;
	auto W5_class = step3->add_edge_differential();
	W5_class.print();
	return W5_class;
}

std::optional<OddGCdegZero<8>> tryFindFullWheel7Class(int generated_rounds) {
	GC loop(loop_graph<13>());

	cout << "loop: ";
	loop.print();

	auto step1 = push_down_by_selected_method(loop, generated_rounds);
	if (!step1) {
		print_missing_step<decltype(loop)>("step 1");
		return std::nullopt;
	}

	auto step2 = push_down_by_selected_method(*step1, generated_rounds);
	if (!step2) {
		print_missing_step<decltype(*step1)>("step 2");
		return std::nullopt;
	}

	auto step3 = push_down_by_selected_method(*step2, generated_rounds);
	if (!step3) {
		print_missing_step<decltype(*step2)>("step 3");
		return std::nullopt;
	}

	auto step4 = push_down_by_selected_method(*step3, generated_rounds);
	if (!step4) {
		print_missing_step<decltype(*step3)>("step 4");
		step3->print();
		return std::nullopt;
	}

	auto step5 = push_down_by_selected_method(*step4, generated_rounds);
	if (!step5) {
		print_missing_step<decltype(*step4)>("step 5");
		return std::nullopt;
	}

	OddGCdegZero<8> W7_class = step5->add_edge_differential();
	return filter_to_wheel_component<7>(W7_class);
}

std::optional<OddGCdegZero<10>> tryFindFullWheel9Class(int generated_rounds) {
	GC loop(loop_graph<17>());

	cout << "loop: ";
	loop.print();

	cout << "are we getting step 1?" << endl;
	auto step1 = push_down_by_selected_method(loop, generated_rounds);
	if (!step1) {
		print_missing_step<decltype(loop)>("step 1");
		return std::nullopt;
	}

	cout << "step 2" << endl;
	auto step2 = push_down_by_selected_method(*step1, generated_rounds);
	if (!step2) {
		print_missing_step<decltype(*step1)>("step 2");
		return std::nullopt;
	}

	cout << "step 3" << endl;
	auto step3 = push_down_by_selected_method(*step2, generated_rounds);
	if (!step3) {
		print_missing_step<decltype(*step2)>("step 3");
		return std::nullopt;
	}

	cout << "trying step 4. step3.size() = " << step3->size() << endl;
	auto step4 = push_down_by_selected_method(*step3, generated_rounds);
	if (!step4) {
		print_missing_step<decltype(*step3)>("step 4");
		return std::nullopt;
	}

	cout << "step 5" << endl;
	auto step5 = push_down_by_selected_method(*step4, generated_rounds);
	if (!step5) {
		print_missing_step<decltype(*step4)>("step 5");
		return std::nullopt;
	}

	cout << "step 6" << endl;
	auto step6 = push_down_by_selected_method(*step5, generated_rounds);
	if (!step6) {
		print_missing_step<decltype(*step5)>("step 6");
		return std::nullopt;
	}

	cout << "step 7" << endl;
	auto step7 = push_down_by_selected_method(*step6, generated_rounds);
	if (!step7) {
		print_missing_step<decltype(*step6)>("step 7");
		return std::nullopt;
	}

	OddGCdegZero<10> W9_class = step7->add_edge_differential();
	return filter_to_wheel_component<9>(W9_class);
}

template <Int WheelN>
int write_split_waterfall_class(
	const string& output_path,
	std::optional<OddGCdegZero<WheelN + 1>> wheel_class
) {
	if (!wheel_class) {
		return EXIT_FAILURE;
	}

	cout << "found W" << +WheelN << " class of size: " << wheel_class->size() << endl;

	ofstream out(output_path);
	if (!out) {
		cerr << "failed to open " << output_path << endl;
		return EXIT_FAILURE;
	}
	wheel_class->print(out);
	return EXIT_SUCCESS;
}

int main(int argc, char** argv) {
	const int wheel_n = argc >= 2 ? std::stoi(argv[1]) : 7;
	const string output_path = argc >= 3 ? argv[2] : "W" + std::to_string(wheel_n) + ".txt";
	const int generated_rounds = argc >= 4 ? std::stoi(argv[3]) : -1;

	switch (wheel_n) {
	case 3:
		return write_split_waterfall_class<3>(output_path, tryFindFullWheel3Class(generated_rounds));
	case 5:
		return write_split_waterfall_class<5>(output_path, tryFindFullWheel5Class(generated_rounds));
	case 7:
		return write_split_waterfall_class<7>(output_path, tryFindFullWheel7Class(generated_rounds));
	case 9:
		return write_split_waterfall_class<9>(output_path, tryFindFullWheel9Class(generated_rounds));
	default:
		cerr << "usage: " << argv[0] << " [3|5|7|9] [output_path] [generated_rounds]" << endl;
		return EXIT_FAILURE;
	}
}
