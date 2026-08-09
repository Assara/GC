#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include "coalgebra_utils.hpp"
#include "examplegraphs.hpp"

namespace {

template <typename GCType>
int run_for_file(const std::filesystem::path& input_path, const std::filesystem::path& output_path) {
	using GraphType = typename GCType::GraphType;

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
		|| !std::getline(in, field_type_line)) {
		std::cerr << "input header is incomplete\n";
		return 1;
	}
	while (std::getline(in, count_line) && count_line.rfind("number_of_graphs:", 0) != 0) {}
	if (count_line.rfind("number_of_graphs:", 0) != 0) {
		std::cerr << "input header is incomplete\n";
		return 1;
	}

	std::regex edge_pattern(R"(\((\d+),(\d+)\))");
	std::vector<typename GCType::Base> kept_terms;

	std::string line;
	std::size_t total_terms = 0;
	std::size_t filtered_terms = 0;
	while (std::getline(in, line)) {
		if (line.empty()) {
			continue;
		}

		const auto semicolon = line.find(';');
		if (semicolon == std::string::npos) {
			std::cerr << "malformed line: " << line << '\n';
			return 1;
		}

		fieldType coeff;
		{
			std::istringstream coeff_stream(line.substr(0, semicolon));
			if (!(coeff_stream >> coeff)) {
				std::cerr << "failed to parse coefficient in line: " << line << '\n';
				return 1;
			}
		}

		GraphType graph;
		std::size_t edge_index = 0;
		const std::string edges_text = line.substr(semicolon + 1);
		for (std::sregex_iterator it(edges_text.begin(), edges_text.end(), edge_pattern), end;
		     it != end;
		     ++it) {
			if (edge_index >= GraphType::N_EDGES_) {
				std::cerr << "too many edges in line: " << line << '\n';
				return 1;
			}
			graph.half_edges[2 * edge_index] = static_cast<Int>(std::stoi((*it)[1].str()));
			graph.half_edges[2 * edge_index + 1] = static_cast<Int>(std::stoi((*it)[2].str()));
			++edge_index;
		}

		if (edge_index != GraphType::N_EDGES_) {
			std::cerr << "expected " << +GraphType::N_EDGES_ << " edges, got " << edge_index << '\n';
			return 1;
		}

		++total_terms;
		if (!coalgebra_utils::detail::is_three_four_valent(graph)) {
			++filtered_terms;
			continue;
		}
		kept_terms.emplace_back(graph, coeff);
	}

	out << "graph_size: (" << +GraphType::N_VERTICES_ << "," << +GraphType::N_EDGES_ << ")\n";
	out << "field_type: " << fieldType::name() << "\n";
	out << "number_of_graphs: " << kept_terms.size() << '\n';
	for (const auto& be : kept_terms) {
		out << be.getCoefficient() << "; ";
		for (Int edge = 0; edge < GraphType::N_EDGES_; ++edge) {
			const auto [u, v] = be.getValue().getEdge(edge);
			out << "(" << +u << "," << +v << ")";
			if (edge + 1 < GraphType::N_EDGES_) {
				out << ", ";
			}
		}
		out << '\n';
	}

	std::cout << "input = " << input_path << '\n';
	std::cout << "output = " << output_path << '\n';
	std::cout << "terms read = " << total_terms << '\n';
	std::cout << "terms kept = " << kept_terms.size() << '\n';
	std::cout << "terms filtered = " << filtered_terms << '\n';
	return 0;
}

void print_usage(const char* argv0) {
	std::cerr << "usage: " << argv0 << " <wheel> <input> <output>\n";
}

} // namespace

int main(int argc, char** argv) {
	if (argc != 4) {
		print_usage(argv[0]);
		return 1;
	}

	const int wheel = std::stoi(argv[1]);
	const std::filesystem::path input_path = argv[2];
	const std::filesystem::path output_path = argv[3];

	switch (wheel) {
	case 7:
		return run_for_file<OddGCdegZero<8>>(input_path, output_path);
	case 9:
		return run_for_file<OddGCdegZero<10>>(input_path, output_path);
	case 11:
		return run_for_file<OddGCdegZero<12>>(input_path, output_path);
	default:
		std::cerr << "unsupported wheel: " << wheel << '\n';
		return 1;
	}
}
