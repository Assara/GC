#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "examplegraphs.hpp"

namespace {

template <typename GraphType>
void write_graph_line(std::ostream& out, const GraphType& graph) {
	for (Int edge = 0; edge < GraphType::N_EDGES_; ++edge) {
		const auto [u, v] = graph.getEdge(edge);
		out << "(" << +u << "," << +v << ")";
		if (edge + 1 < GraphType::N_EDGES_) {
			out << ", ";
		}
	}
	out << '\n';
}

std::uint64_t parse_number_of_graphs(const std::string& line) {
	const std::string prefix = "number_of_graphs: ";
	if (!line.starts_with(prefix)) {
		throw std::runtime_error("expected number_of_graphs header");
	}
	return std::stoull(line.substr(prefix.size()));
}

template <typename GraphType>
GraphType parse_graph_line(const std::string& line) {
	GraphType graph;
	std::size_t pos = 0;

	for (Int edge = 0; edge < GraphType::N_EDGES_; ++edge) {
		pos = line.find('(', pos);
		if (pos == std::string::npos) {
			throw std::runtime_error("missing '(' while parsing graph line");
		}
		++pos;

		std::size_t next_pos = 0;
		const int u = std::stoi(line.substr(pos), &next_pos);
		pos += next_pos;

		if (pos >= line.size() || line[pos] != ',') {
			throw std::runtime_error("missing ',' while parsing graph line");
		}
		++pos;

		const int v = std::stoi(line.substr(pos), &next_pos);
		pos += next_pos;

		if (pos >= line.size() || line[pos] != ')') {
			throw std::runtime_error("missing ')' while parsing graph line");
		}
		++pos;

		graph.setEdge(edge, static_cast<Int>(u), static_cast<Int>(v));
	}

	return graph;
}

template <typename GraphType, typename Rng>
GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_> random_isomorphism(Rng& rng) {
	GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_> iso;

	auto& vertex_perm = iso.vertex_permutation_data();
	for (Int v = 0; v < GraphType::N_VERTICES_; ++v) {
		vertex_perm[v] = v;
	}
	std::shuffle(vertex_perm.begin(), vertex_perm.end(), rng);

	auto& edge_perm = iso.edge_permutation_data();
	for (Int e = 0; e < GraphType::N_EDGES_; ++e) {
		edge_perm[e] = e;
	}
	std::shuffle(edge_perm.begin(), edge_perm.end(), rng);

	std::bernoulli_distribution flip_dist(0.5);
	auto& edge_flip = iso.edge_flip_data();
	for (Int e = 0; e < GraphType::N_EDGES_; ++e) {
		edge_flip[e] = flip_dist(rng);
	}

	return iso;
}

template <typename GraphType>
int run_case(
	const std::filesystem::path& input_path,
	const std::filesystem::path& output_path,
	std::uint64_t seed
) {
	std::ifstream in(input_path);
	if (!in) {
		std::cerr << "failed to open input: " << input_path << '\n';
		return 1;
	}

	if (output_path.has_parent_path()) {
		std::filesystem::create_directories(output_path.parent_path());
	}
	std::ofstream out(output_path, std::ios::trunc);
	if (!out) {
		std::cerr << "failed to open output: " << output_path << '\n';
		return 1;
	}

	std::string graph_size_line;
	std::string field_type_line;
	std::string count_line;
	if (!std::getline(in, graph_size_line)
		|| !std::getline(in, field_type_line)
		|| !std::getline(in, count_line)) {
		std::cerr << "input header is incomplete\n";
		return 1;
	}

	const std::uint64_t total_graphs = parse_number_of_graphs(count_line);
	std::mt19937_64 rng(seed);

	out << "graph_size: (" << +GraphType::N_VERTICES_ << "," << +GraphType::N_EDGES_ << ")\n";
	out << "field_type: " << TypeName<fieldType>::name() << "\n";
	out << "number_of_graphs: " << total_graphs << "\n";

	std::string line;
	std::uint64_t written = 0;
	while (std::getline(in, line)) {
		if (line.empty()) {
			continue;
		}
		GraphType graph = parse_graph_line<GraphType>(line);
		const auto iso = random_isomorphism<GraphType>(rng);
		write_graph_line(out, iso.permute(graph));
		++written;
	}

	if (written != total_graphs) {
		std::cerr << "warning: header count = " << total_graphs
		          << ", written = " << written << '\n';
	}

	std::cout << "input = " << input_path << '\n';
	std::cout << "output = " << output_path << '\n';
	std::cout << "seed = " << seed << '\n';
	std::cout << "graphs written = " << written << '\n';
	return 0;
}

void print_usage(const char* argv0) {
	std::cerr << "usage: " << argv0 << " <wheel> <input> <output> [seed]\n";
}

} // namespace

int main(int argc, char** argv) {
	if (argc < 4 || argc > 5) {
		print_usage(argv[0]);
		return 1;
	}

	const int wheel = std::stoi(argv[1]);
	const std::filesystem::path input_path = argv[2];
	const std::filesystem::path output_path = argv[3];
	const std::uint64_t seed = argc >= 5 ? std::stoull(argv[4]) : 123456789ULL;

	switch (wheel) {
	case 7:
		return run_case<OddGraphdegZero<8>>(input_path, output_path, seed);
	case 9:
		return run_case<OddGraphdegZero<10>>(input_path, output_path, seed);
	case 11:
		return run_case<OddGraphdegZero<12>>(input_path, output_path, seed);
	default:
		std::cerr << "unsupported wheel: " << wheel << '\n';
		return 1;
	}
}
