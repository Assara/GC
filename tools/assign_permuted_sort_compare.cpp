#include <chrono>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "examplegraphs.hpp"

namespace {

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
	graphs.reserve(static_cast<std::size_t>(1) << maximal_v_sequence_length);

	for (auto sequence : left_right_sequences(maximal_v_sequence_length)) {
		graphs.push_back(V_graph<N>(sequence));
	}

	std::sort(graphs.begin(), graphs.end());
	return graphs;
}

template <Int N>
std::vector<OddGraphdegZero<N + 1>> make_split_contract_workload_graphs(int rounds) {
	using WorkGraph = OddGraphdegZero<N + 1>;

	std::vector<WorkGraph> graphs;
	graphs.push_back(wheel_graph<N>());

	for (int round = 0; round < rounds; ++round) {
		std::unordered_set<WorkGraph> next_graphs;

		for (const auto& graph : graphs) {
			const auto splits = graph.unsorted_splits(fieldType{1});
			for (const auto& split_be : splits.raw_elements()) {
				if (split_be.getCoefficient() == fieldType{}) {
					continue;
				}

				for (Int edge = 0; edge < WorkGraph::SplitGraph::N_EDGES_; ++edge) {
					auto contracted = split_be.getValue().contract_edge(
						edge,
						split_be.getCoefficient()
					);
					auto& contracted_graph = contracted.getValue();
					const signedInt sort_sign = contracted_graph.directAndSortEdges();
					if (sort_sign != 0 && !contracted_graph.has_double_edge()) {
						next_graphs.emplace(contracted.getValue());
					}
				}
			}
		}

		graphs.assign(next_graphs.begin(), next_graphs.end());
		std::sort(graphs.begin(), graphs.end());
	}

	return graphs;
}

template <typename GraphType>
std::vector<Permutation<GraphType::N_VERTICES_>> benchmark_permutations() {
	using PermType = Permutation<GraphType::N_VERTICES_>;
	std::vector<PermType> perms(3);

	for (Int v = 0; v < GraphType::N_VERTICES_; ++v) {
		perms[0].p[v] = v;
	}

	perms[1].p[0] = 0;
	for (Int v = 1; v < GraphType::N_VERTICES_; ++v) {
		perms[1].p[v] = (v + 1 == GraphType::N_VERTICES_) ? 1 : static_cast<Int>(v + 1);
	}

	perms[2].p[0] = 0;
	perms[2].p[1] = 1;
	for (Int v = 2; v < GraphType::N_VERTICES_; ++v) {
		perms[2].p[v] = static_cast<Int>(GraphType::N_VERTICES_ + 1 - v);
	}

	return perms;
}

template <typename Fn>
double time_seconds(Fn&& fn) {
	const auto start = std::chrono::steady_clock::now();
	fn();
	const auto stop = std::chrono::steady_clock::now();
	return std::chrono::duration<double>(stop - start).count();
}

template <typename GraphType>
std::uint64_t graph_checksum(const GraphType& graph) {
	std::uint64_t checksum = 0;
	for (Int i = 0; i < GraphType::SIZE; ++i) {
		checksum = checksum * 1315423911ULL
			+ static_cast<std::uint64_t>(graph.half_edges[i] + 17);
	}
	return checksum;
}

template <typename GraphType, typename AssignFn>
double benchmark_assign(
	const std::vector<GraphType>& graphs,
	const std::vector<Permutation<GraphType::N_VERTICES_>>& perms,
	int iterations,
	std::uint64_t& checksum,
	AssignFn assign
) {
	checksum = 0;
	return time_seconds([&]() {
		for (int iteration = 0; iteration < iterations; ++iteration) {
			for (const GraphType& source : graphs) {
				for (const auto& perm : perms) {
					GraphType output;
					const signedInt sign = assign(output, source, perm);
					checksum += static_cast<std::uint64_t>(sign + 2);
					checksum ^= graph_checksum(output);
				}
			}
		}
	});
}

template <Int N>
int run_case(int rounds, int iterations, const std::string& workload_kind) {
	using GraphType = OddGraphdegZero<N + 1>;

	std::vector<GraphType> graphs;
	if (workload_kind == "vmax") {
		graphs = make_maximal_v_workload_graphs<N>();
	} else if (workload_kind == "split-contract") {
		graphs = make_split_contract_workload_graphs<N>(rounds);
	} else {
		std::cerr << "unknown workload kind = " << workload_kind << '\n';
		return EXIT_FAILURE;
	}

	const auto perms = benchmark_permutations<GraphType>();
	const std::uint64_t calls =
		static_cast<std::uint64_t>(iterations)
		* static_cast<std::uint64_t>(graphs.size())
		* static_cast<std::uint64_t>(perms.size());

	std::uint64_t old_path_checksum = 0;
	const double old_path_seconds = benchmark_assign(
		graphs,
		perms,
		iterations,
		old_path_checksum,
		[](GraphType& output, const GraphType& source, const auto& perm) {
			output = source;
			return output.permuteVertices(perm)
				* output.directEdges()
				* output.sortEdgesQuick();
		}
	);

	std::uint64_t fused_checksum = 0;
	const double fused_seconds = benchmark_assign(
		graphs,
		perms,
		iterations,
		fused_checksum,
		[](GraphType& output, const GraphType& source, const auto& perm) {
			return output.assignPermutedDirectedSortedEdges(source, perm);
		}
	);

	std::cout << "wheel = W" << +N << '\n';
	std::cout << "workload kind = " << workload_kind << '\n';
	std::cout << "rounds = " << rounds << '\n';
	std::cout << "graphs = " << graphs.size() << '\n';
	std::cout << "permutations = " << perms.size() << '\n';
	std::cout << "iterations = " << iterations << '\n';
	std::cout << "calls = " << calls << '\n';
	std::cout << "old path checksum = " << old_path_checksum << '\n';
	std::cout << "fused checksum = " << fused_checksum << '\n';
	std::cout << "checksums match = " << (old_path_checksum == fused_checksum ? "yes" : "no") << '\n';
	std::cout << "permute + direct + sortEdgesQuick = " << old_path_seconds << " s\n";
	std::cout << "assignPermutedDirectedSortedEdges = " << fused_seconds << " s\n";
	if (fused_seconds > 0) {
		std::cout << "old / fused = " << (old_path_seconds / fused_seconds) << '\n';
	}

	return old_path_checksum == fused_checksum ? EXIT_SUCCESS : EXIT_FAILURE;
}

void print_usage(const char* argv0) {
	std::cerr << "usage: " << argv0 << " wheel_N rounds iterations workload_kind\n";
	std::cerr << "  workload_kind: vmax or split-contract\n";
}

} // namespace

int main(int argc, char** argv) {
	if (argc != 5) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	const int wheel_n = std::stoi(argv[1]);
	const int rounds = std::stoi(argv[2]);
	const int iterations = std::stoi(argv[3]);
	const std::string workload_kind = argv[4];

	if (rounds < 0 || iterations <= 0) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	switch (wheel_n) {
		case 7: return run_case<7>(rounds, iterations, workload_kind);
		case 9: return run_case<9>(rounds, iterations, workload_kind);
		case 11: return run_case<11>(rounds, iterations, workload_kind);
		case 13: return run_case<13>(rounds, iterations, workload_kind);
		case 15: return run_case<15>(rounds, iterations, workload_kind);
		case 17: return run_case<17>(rounds, iterations, workload_kind);
		case 21: return run_case<21>(rounds, iterations, workload_kind);
		case 25: return run_case<25>(rounds, iterations, workload_kind);
		case 27: return run_case<27>(rounds, iterations, workload_kind);
		default:
			print_usage(argv[0]);
			return EXIT_FAILURE;
	}
}
