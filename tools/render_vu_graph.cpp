#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "examplegraphs.hpp"

namespace {

struct Point {
	double x;
	double y;
};

std::vector<bool> parse_sequence(const std::string& text) {
	std::vector<bool> sequence;
	sequence.reserve(text.size());
	for (char ch : text) {
		switch (ch) {
			case '0':
			case 'L':
			case 'l':
				sequence.push_back(false);
				break;
			case '1':
			case 'R':
			case 'r':
				sequence.push_back(true);
				break;
			case ',':
			case ' ':
			case '_':
			case '-':
				break;
			default:
				throw std::runtime_error("sequence may only use 0/1 or L/R");
		}
	}
	return sequence;
}

template <typename GraphType>
void write_tikz_picture(
	const std::filesystem::path& output_path,
	const GraphType& graph,
	const std::vector<Point>& positions
) {
	constexpr int n_vertices = std::tuple_size_v<decltype(graph.valence_array())>;
	constexpr int n_edges = GraphType::N_EDGES_;

	std::ofstream out(output_path);
	out << std::fixed << std::setprecision(3);
	out << "\\begin{tikzpicture}[x=0.8cm,y=0.8cm,baseline=(current bounding box.center)]\n";
	out << "  \\tikzset{graphvertex/.style={circle,fill=black,draw=black,inner sep=0pt,minimum size=2.2pt}}\n";
	out << "  \\tikzset{graphedge/.style={draw=black,line width=0.4pt}}\n";
	out << "  \\tikzset{edgelabel/.style={font=\\scriptsize,inner sep=1pt,fill=white,text=black}}\n";

	for (int v = 0; v < n_vertices; ++v) {
		out << "  \\node[graphvertex] (v" << v << ") at (" << positions[v].x << "," << positions[v].y << ") {};\n";
	}

	for (int e = 0; e < n_edges; ++e) {
		auto [u, v] = graph.getEdge(static_cast<Int>(e));
		const double mx = 0.5 * (positions[static_cast<int>(u)].x + positions[static_cast<int>(v)].x);
		const double my = 0.5 * (positions[static_cast<int>(u)].y + positions[static_cast<int>(v)].y);
		const double dx = positions[static_cast<int>(v)].x - positions[static_cast<int>(u)].x;
		const double dy = positions[static_cast<int>(v)].y - positions[static_cast<int>(u)].y;
		double nx = -dy;
		double ny = dx;
		const double norm = std::sqrt(nx * nx + ny * ny);
		if (norm > 0.0) {
			nx /= norm;
			ny /= norm;
		}
		const bool curved_closing_edge = std::abs(dx) > 5.0;
		if (curved_closing_edge) {
			const double xu = positions[static_cast<int>(u)].x;
			const double yu = positions[static_cast<int>(u)].y;
			const double xv = positions[static_cast<int>(v)].x;
			const double yv = positions[static_cast<int>(v)].y;
			const double ctrl_y = std::min(yu, yv) - 2.0;
			const double c1x = xu + 1.8;
			const double c2x = xv - 1.8;
			out << "  \\draw[graphedge] (v" << static_cast<int>(u) << ") .. controls (" << c1x << "," << ctrl_y << ") and (" << c2x << "," << ctrl_y << ") .. (v" << static_cast<int>(v) << ");\n";
			out << "  \\node[edgelabel] at (" << mx << "," << (ctrl_y + 0.25) << ") {" << e << "};\n";
		} else {
			const double ox = 0.16 * nx;
			const double oy = 0.16 * ny;
			out << "  \\draw[graphedge] (v" << static_cast<int>(u) << ") -- (v" << static_cast<int>(v) << ");\n";
			out << "  \\node[edgelabel] at (" << (mx + ox) << "," << (my + oy) << ") {" << e << "};\n";
		}
	}

	out << "\\end{tikzpicture}\n";
}

std::vector<Point> build_positions(const std::vector<bool>& sequence, bool include_extra_center, int total_vertices) {
	std::vector<Point> positions(total_vertices, {0.0, 0.0});

	const double left_x = -3.0;
	const double center_x = 0.0;
	const double right_x = 3.0;
	const double top_y = 2.2;
	const double base_y = 0.0;
	const double row_step = 2.2;
	const double pi = std::acos(-1.0);

	positions[0] = {center_x, top_y};
	positions[1] = {left_x, base_y};
	positions[2] = {center_x, base_y};
	positions[3] = {right_x, base_y};

	int left = 1;
	int center = 2;
	int right = 3;
	int next = 4;

	for (std::size_t i = 0; i < sequence.size(); ++i) {
		const double row_y = base_y - row_step * static_cast<double>(i + 1);

		const int new_center = next++;
		positions[new_center] = {center_x, row_y};
		center = new_center;

		const int new_side = next++;
		if (sequence[i]) {
			positions[new_side] = {right_x, row_y};
			right = new_side;
		} else {
			positions[new_side] = {left_x, row_y};
			left = new_side;
		}
	}

	if (include_extra_center) {
		const double row_y = base_y - row_step * static_cast<double>(sequence.size() + 1);
		const int new_center = next++;
		positions[new_center] = {center_x, row_y};
		center = new_center;
	}

	(void)center;

	const int remaining = total_vertices - next;
	if (remaining > 0) {
		const Point a = positions[left];
		const Point b = positions[right];
		const double mx = 0.5 * (a.x + b.x);
		const double my = 0.5 * (a.y + b.y);
		const double dx = b.x - a.x;
		const double dy = b.y - a.y;
		const double chord = std::sqrt(dx * dx + dy * dy);
		if (chord == 0.0) {
			throw std::runtime_error("arc endpoints coincide in placement algorithm");
		}
		const double ux = dx / chord;
		const double uy = dy / chord;
		double nx = uy;
		double ny = -ux;
		if (ny > 0.0) {
			nx = -nx;
			ny = -ny;
		}
		const double radius = 0.5 * chord;

		const double arc_y_shift = include_extra_center ? -row_step : 0.0;
		for (int j = 0; j < remaining; ++j) {
			const double theta = pi - (pi * static_cast<double>(j + 1)) / static_cast<double>(remaining + 1);
			const double x = mx + radius * std::cos(theta) * ux + radius * std::sin(theta) * nx;
			const double y = my + radius * std::cos(theta) * uy + radius * std::sin(theta) * ny + arc_y_shift;
			positions[next + j] = {x, y};
		}
	}

	return positions;
}

template <Int N>
int dispatch_render(
	const std::string& family,
	const std::vector<bool>& sequence,
	const std::filesystem::path& output_path,
	const std::string& graph_name
) {
	if (family == "V" || family == "v") {
		auto graph = V_graph<N>(const_cast<std::vector<bool>&>(sequence));
		constexpr int n_vertices = std::tuple_size_v<decltype(graph.valence_array())>;
		auto tikz_path = output_path;
		tikz_path.replace_extension(".tex");
		write_tikz_picture(tikz_path, graph, build_positions(sequence, false, n_vertices));
		return 0;
	}
	if (family == "U" || family == "u") {
		auto graph = U_graph<N>(const_cast<std::vector<bool>&>(sequence));
		constexpr int n_vertices = std::tuple_size_v<decltype(graph.valence_array())>;
		auto tikz_path = output_path;
		tikz_path.replace_extension(".tex");
		write_tikz_picture(tikz_path, graph, build_positions(sequence, true, n_vertices));
		return 0;
	}
	throw std::runtime_error("family must be V or U");
}

int render_for_n(
	int n,
	const std::string& family,
	const std::vector<bool>& sequence,
	const std::filesystem::path& output_path,
	const std::string& graph_name
) {
	switch (n) {
		case 5: return dispatch_render<5>(family, sequence, output_path, graph_name);
		case 7: return dispatch_render<7>(family, sequence, output_path, graph_name);
		case 9: return dispatch_render<9>(family, sequence, output_path, graph_name);
		case 11: return dispatch_render<11>(family, sequence, output_path, graph_name);
		case 13: return dispatch_render<13>(family, sequence, output_path, graph_name);
		case 15: return dispatch_render<15>(family, sequence, output_path, graph_name);
		case 17: return dispatch_render<17>(family, sequence, output_path, graph_name);
		case 19: return dispatch_render<19>(family, sequence, output_path, graph_name);
		case 21: return dispatch_render<21>(family, sequence, output_path, graph_name);
		case 23: return dispatch_render<23>(family, sequence, output_path, graph_name);
		case 25: return dispatch_render<25>(family, sequence, output_path, graph_name);
		case 27: return dispatch_render<27>(family, sequence, output_path, graph_name);
		case 29: return dispatch_render<29>(family, sequence, output_path, graph_name);
		default:
			throw std::runtime_error("unsupported N for renderer");
	}
}

} // namespace

int main(int argc, char** argv) {
	if (argc < 5 || argc > 6) {
		std::cerr << "usage: render_vu_graph <V|U> <odd N> <sequence> <output-base> [graph_name]\n";
		return 1;
	}

	try {
		const std::string family = argv[1];
		const int n = std::stoi(argv[2]);
		const auto sequence = parse_sequence(argv[3]);
		const std::filesystem::path output_path = argv[4];
		const std::string graph_name = argc == 6 ? argv[5] : (family + std::to_string(n));

		std::filesystem::create_directories(output_path.parent_path());
		return render_for_n(n, family, sequence, output_path, graph_name);
	} catch (const std::exception& ex) {
		std::cerr << "render_vu_graph: " << ex.what() << '\n';
		return 1;
	}
}
