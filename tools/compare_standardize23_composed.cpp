#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "examplegraphs.hpp"

namespace {

enum class WorkloadKind {
	SplitContract,
	MaximalV
};

WorkloadKind parse_workload_kind(const std::string& value) {
	if (value == "split-contract" || value == "split_contract" || value == "split") {
		return WorkloadKind::SplitContract;
	}
	if (value == "maximal-v" || value == "maximal_v" || value == "vmax" || value == "v") {
		return WorkloadKind::MaximalV;
	}
	throw std::invalid_argument("unknown workload kind: " + value);
}

std::vector<std::vector<bool>> left_right_sequences(int length) {
	std::vector<std::vector<bool>> result;
	if (length < 0) {
		return result;
	}

	const int count = 1 << length;
	result.reserve(count);

	for (int mask = 0; mask < count; ++mask) {
		std::vector<bool> seq;
		seq.reserve(length);
		for (int bit = 0; bit < length; ++bit) {
			seq.push_back(((mask >> bit) & 1) != 0);
		}
		result.push_back(std::move(seq));
	}

	return result;
}

template <Int N>
std::vector<OddGraphdegZero<N + 1>> make_maximal_v_workload_graphs() {
	std::vector<OddGraphdegZero<N + 1>> graphs;
	const int maximal_v_sequence_length = (N - 3) / 2;
	const int sequence_count = 1 << maximal_v_sequence_length;
	graphs.reserve(sequence_count);

	for (auto sequence : left_right_sequences(maximal_v_sequence_length)) {
		graphs.push_back(V_graph<N>(sequence));
	}

	std::sort(graphs.begin(), graphs.end());
	return graphs;
}

template <Int N>
std::string split_contract_cache_path(int rounds) {
	return "/tmp/gc_standardizer_split_contract_W" + std::to_string(N)
		+ "_rounds" + std::to_string(rounds)
		+ "_int" + std::to_string(sizeof(Int))
		+ ".bin";
}

template <typename GraphType>
bool load_graph_workload_cache(const std::string& path, std::vector<GraphType>& graphs) {
	std::ifstream in(path, std::ios::binary);
	if (!in) {
		return false;
	}

	std::uint64_t count = 0;
	in.read(reinterpret_cast<char*>(&count), sizeof(count));
	if (!in) {
		return false;
	}

	graphs.resize(static_cast<std::size_t>(count));
	for (auto& graph : graphs) {
		in.read(
			reinterpret_cast<char*>(graph.half_edges.data()),
			static_cast<std::streamsize>(graph.half_edges.size() * sizeof(Int))
		);
		if (!in) {
			graphs.clear();
			return false;
		}
	}

	return true;
}

template <Int N>
std::vector<OddGraphdegZero<N + 1>> load_split_contract_workload_graphs(int rounds) {
	using WorkGraph = OddGraphdegZero<N + 1>;

	std::vector<WorkGraph> graphs;
	const std::string cache_path = split_contract_cache_path<N>(rounds);
	if (!load_graph_workload_cache(cache_path, graphs)) {
		throw std::runtime_error("missing workload cache: " + cache_path);
	}

	return graphs;
}

template <typename GraphType>
bool equal_basis_elements(
	const BasisElement<GraphType, fieldType>& lhs,
	const BasisElement<GraphType, fieldType>& rhs
) {
	return lhs.getValue().compare(rhs.getValue()) == 0 &&
	       lhs.getCoefficient() == rhs.getCoefficient();
}

template <typename GraphType>
void compare_compositions(
	const std::vector<GraphType>& graphs,
	const char* label
) {
	using Element = BasisElement<GraphType, fieldType>;
	using Standardizer = GraphStandardizer<
		GraphType::N_VERTICES_,
		GraphType::N_EDGES_,
		GraphType::N_OUT_HAIR_,
		GraphType::N_IN_HAIR_,
		GraphType::C_,
		GraphType::D_,
		fieldType
	>;

	Standardizer standardizer;
	std::size_t mismatch_32_vs_2 = 0;
	std::size_t mismatch_23_vs_3 = 0;

	for (const auto& graph : graphs) {
		Element input2(graph, fieldType{1});
		Element input3(graph, fieldType{1});

		Element result2 = standardizer.standardize2(input2);
		Element result3 = standardizer.standardize3(input3);
		Element result32_input = result3;
		Element result23_input = result2;
		const Element result32 = standardizer.standardize2(result32_input);
		const Element result23 = standardizer.standardize3(result23_input);

		if (!equal_basis_elements(result32, result2)) {
			++mismatch_32_vs_2;
		}
		if (!equal_basis_elements(result23, result3)) {
			++mismatch_23_vs_3;
		}
	}

	std::cout << label << '\n';
	std::cout << "graphs = " << graphs.size() << '\n';
	std::cout << "standardize2(standardize3(g)) != standardize2(g): "
	          << mismatch_32_vs_2 << '\n';
	std::cout << "standardize3(standardize2(g)) != standardize3(g): "
	          << mismatch_23_vs_3 << '\n';
}

template <Int N>
int run_case(int rounds, WorkloadKind workload_kind) {
	std::vector<OddGraphdegZero<N + 1>> graphs =
		(workload_kind == WorkloadKind::MaximalV)
			? make_maximal_v_workload_graphs<N>()
			: load_split_contract_workload_graphs<N>(rounds);

	const std::string label =
		(workload_kind == WorkloadKind::MaximalV)
			? ("W" + std::to_string(N) + " maximal V")
			: ("W" + std::to_string(N) + " split-contract rounds=" + std::to_string(rounds));

	compare_compositions(graphs, label.c_str());
	return EXIT_SUCCESS;
}

void print_usage(const char* argv0) {
	std::cerr << "usage: " << argv0 << " [wheel_N] [split_contract_rounds] [workload]\n";
	std::cerr << "  workload: split-contract (default) or vmax\n";
	std::cerr << "  wheel_N must be one of 11,25,27\n";
}

} // namespace

int main(int argc, char** argv) {
	const int wheel_n = argc >= 2 ? std::stoi(argv[1]) : 25;
	const int rounds = argc >= 3 ? std::stoi(argv[2]) : 2;
	WorkloadKind workload_kind = WorkloadKind::SplitContract;
	try {
		workload_kind = argc >= 4 ? parse_workload_kind(argv[3]) : WorkloadKind::SplitContract;
	} catch (const std::invalid_argument& error) {
		std::cerr << error.what() << '\n';
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	switch (wheel_n) {
	case 11:
		return run_case<11>(rounds, workload_kind);
	case 25:
		return run_case<25>(rounds, workload_kind);
	case 27:
		return run_case<27>(rounds, workload_kind);
	default:
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}
}
