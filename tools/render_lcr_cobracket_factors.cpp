#include <filesystem>

#include "GraphVisualizer.hpp"
#include "graph.hpp"

namespace {

using LeftGraph = Graph<6, 10, 0, 0, 0, 1, fieldType>;
using RightGraph = Graph<8, 14, 0, 0, 0, 1, fieldType>;

LeftGraph make_left_graph() {
	LeftGraph graph;
	graph.setEdge(0, 0, 1);
	graph.setEdge(1, 0, 2);
	graph.setEdge(2, 0, 3);
	graph.setEdge(3, 0, 4);
	graph.setEdge(4, 1, 2);
	graph.setEdge(5, 1, 3);
	graph.setEdge(6, 1, 5);
	graph.setEdge(7, 2, 4);
	graph.setEdge(8, 3, 5);
	graph.setEdge(9, 4, 5);
	return graph;
}

RightGraph make_right_graph_1() {
	RightGraph graph;
	graph.setEdge(0, 0, 1);
	graph.setEdge(1, 0, 2);
	graph.setEdge(2, 0, 3);
	graph.setEdge(3, 0, 4);
	graph.setEdge(4, 1, 2);
	graph.setEdge(5, 1, 3);
	graph.setEdge(6, 1, 5);
	graph.setEdge(7, 2, 4);
	graph.setEdge(8, 2, 6);
	graph.setEdge(9, 3, 5);
	graph.setEdge(10, 3, 7);
	graph.setEdge(11, 4, 6);
	graph.setEdge(12, 5, 7);
	graph.setEdge(13, 6, 7);
	return graph;
}

RightGraph make_right_graph_2() {
	RightGraph graph;
	graph.setEdge(0, 0, 1);
	graph.setEdge(1, 0, 2);
	graph.setEdge(2, 0, 3);
	graph.setEdge(3, 0, 4);
	graph.setEdge(4, 1, 2);
	graph.setEdge(5, 1, 3);
	graph.setEdge(6, 1, 5);
	graph.setEdge(7, 2, 4);
	graph.setEdge(8, 2, 6);
	graph.setEdge(9, 3, 5);
	graph.setEdge(10, 4, 7);
	graph.setEdge(11, 5, 6);
	graph.setEdge(12, 5, 7);
	graph.setEdge(13, 6, 7);
	return graph;
}

} // namespace

int main() {
	namespace fs = std::filesystem;
	const fs::path out_dir = "docs/notice/figures";
	fs::create_directories(out_dir);

	GraphVisualizer::Options options;
	options.show_edge_labels = false;
	options.show_valence_labels = true;
	options.color_by_valence = true;

	options.graph_name = "cobracket_left";
	GraphVisualizer::write_dot(out_dir / "cobracket_left.dot", make_left_graph(), options);

	options.graph_name = "cobracket_right_1";
	GraphVisualizer::write_dot(out_dir / "cobracket_right_1.dot", make_right_graph_1(), options);

	options.graph_name = "cobracket_right_2";
	GraphVisualizer::write_dot(out_dir / "cobracket_right_2.dot", make_right_graph_2(), options);

	return 0;
}
