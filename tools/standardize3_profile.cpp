#include <chrono>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#include "examplegraphs.hpp"

namespace {

std::vector<std::vector<bool>> left_right_sequences_starting_left(int length) {
	std::vector<std::vector<bool>> result;
	if (length <= 0) {
		return result;
	}

	const int tail_length = length - 1;
	const int count = 1 << tail_length;
	result.reserve(count);

	for (int mask = 0; mask < count; ++mask) {
		std::vector<bool> seq;
		seq.reserve(length);
		seq.push_back(false);
		for (int bit = 0; bit < tail_length; ++bit) {
			seq.push_back(((mask >> bit) & 1) != 0);
		}
		result.push_back(std::move(seq));
	}

	return result;
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
		+ "_v2"
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

template <typename GraphType>
void save_graph_workload_cache(const std::string& path, const std::vector<GraphType>& graphs) {
	std::ofstream out(path, std::ios::binary);
	if (!out) {
		return;
	}

	const std::uint64_t count = graphs.size();
	out.write(reinterpret_cast<const char*>(&count), sizeof(count));
	for (const auto& graph : graphs) {
		out.write(
			reinterpret_cast<const char*>(graph.half_edges.data()),
			static_cast<std::streamsize>(graph.half_edges.size() * sizeof(Int))
		);
	}
}

template <Int N>
std::vector<OddGraphdegZero<N + 1>> load_or_make_split_contract_workload_graphs(int rounds) {
	using WorkGraph = OddGraphdegZero<N + 1>;

	std::vector<WorkGraph> graphs;
	const std::string cache_path = split_contract_cache_path<N>(rounds);
	if (load_graph_workload_cache(cache_path, graphs)) {
		std::cout << "loaded workload cache = " << cache_path << '\n';
		return graphs;
	}

	graphs = make_split_contract_workload_graphs<N>(rounds);
	save_graph_workload_cache(cache_path, graphs);
	std::cout << "saved workload cache = " << cache_path << '\n';
	return graphs;
}

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

template <typename Fn>
double time_seconds(Fn&& fn) {
	const auto start = std::chrono::steady_clock::now();
	fn();
	const auto stop = std::chrono::steady_clock::now();
	return std::chrono::duration<double>(stop - start).count();
}

template <typename GraphType>
double benchmark_standardize3(const std::vector<GraphType>& graphs, int iterations) {
	using Element = BasisElement<GraphType, fieldType>;
	double total = 0.0;
	std::size_t checksum = 0;

	for (int i = 0; i < iterations; ++i) {
		total += time_seconds([&]() {
			for (const auto& graph : graphs) {
				Element be(graph, fieldType{1});
				auto canon = GraphType::canonized(be);
				checksum += canon.getValue().half_edges[0];
			}
		});
	}

	std::cout << "checksum = " << checksum << '\n';
	return total / static_cast<double>(iterations);
}

template <Int N>
int run_standardize3_profile(
	int expansion_rounds,
	int repeat,
	int iterations,
	WorkloadKind workload_kind
) {
	const auto base_graphs = workload_kind == WorkloadKind::MaximalV
		? make_maximal_v_workload_graphs<N>()
		: load_or_make_split_contract_workload_graphs<N>(expansion_rounds);
	std::vector<OddGraphdegZero<N + 1>> graphs;
	graphs.reserve(base_graphs.size() * repeat);
	for (int r = 0; r < repeat; ++r) {
		graphs.insert(graphs.end(), base_graphs.begin(), base_graphs.end());
	}

	std::cout << "standardize3 profile\n";
	std::cout << "wheel = W" << +N << '\n';
	if (workload_kind == WorkloadKind::MaximalV) {
		std::cout << "workload kind = maximal V_N sequences\n";
		std::cout << "maximal sequence length = " << ((N - 3) / 2) << '\n';
		std::cout << "base maximal V graphs = " << base_graphs.size() << '\n';
	} else {
		std::cout << "workload kind = split-contract\n";
		std::cout << "split-contract rounds = " << expansion_rounds << '\n';
		std::cout << "base split-contract graphs = " << base_graphs.size() << '\n';
	}
	std::cout << "repeat = " << repeat << '\n';
	std::cout << "workload graphs = " << graphs.size() << '\n';
	std::cout << "iterations = " << iterations << '\n';

#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
	gc_standardizer_sort_profile::reset();
#endif
	const double avg = benchmark_standardize3(graphs, iterations);
	std::cout << "standardize3 average = " << avg << " s\n";

#if defined(GC_PROFILE_STANDARDIZER_SORT) || defined(GC_PROFILE_STANDARDIZER_LABELING_ONLY)
	const double labeling_search_total_seconds =
		static_cast<double>(gc_standardizer_sort_profile::labeling_search_nanoseconds.load(std::memory_order_relaxed)) / 1'000'000'000.0;
	const double labeling_search_avg_seconds = labeling_search_total_seconds / static_cast<double>(iterations);
	const auto attempts_created = gc_standardizer_sort_profile::attempts_created.load(std::memory_order_relaxed);
	const auto true_final_attempts = gc_standardizer_sort_profile::true_final_attempts.load(std::memory_order_relaxed);
	const auto false_final_attempts = gc_standardizer_sort_profile::false_final_attempts.load(std::memory_order_relaxed);
	const auto total_final_attempts = true_final_attempts + false_final_attempts;
	std::cout << "labeling search average = " << labeling_search_avg_seconds << " s\n";
	std::cout << "attempts created = " << attempts_created << '\n';
	std::cout << "true final attempts = " << true_final_attempts << '\n';
	std::cout << "false final attempts = " << false_final_attempts << '\n';
	if (total_final_attempts > 0) {
		std::cout << "true/false final ratio = "
		          << (static_cast<double>(true_final_attempts) / static_cast<double>(false_final_attempts == 0 ? 1 : false_final_attempts))
		          << '\n';
		std::cout << "true final share = "
		          << (100.0 * static_cast<double>(true_final_attempts) / static_cast<double>(total_final_attempts))
		          << "%\n";
	}
#endif

#if defined(GC_PROFILE_STANDARDIZER_SORT)
	const double init_colors_total_seconds =
		static_cast<double>(gc_standardizer_sort_profile::init_colors_nanoseconds.load(std::memory_order_relaxed)) / 1'000'000'000.0;
	const double init_colors_avg_seconds =
		init_colors_total_seconds / static_cast<double>(iterations);
	const double update_colors_total_seconds =
		static_cast<double>(gc_standardizer_sort_profile::vertex_color_update_nanoseconds.load(std::memory_order_relaxed)) / 1'000'000'000.0;
	const double update_colors_avg_seconds =
		update_colors_total_seconds / static_cast<double>(iterations);
	const double update_groups_total_seconds =
		static_cast<double>(gc_standardizer_sort_profile::vertex_bucket_init_nanoseconds.load(std::memory_order_relaxed)) / 1'000'000'000.0;
	const double update_groups_avg_seconds =
		update_groups_total_seconds / static_cast<double>(iterations);
	const double push_next_attempts_total_seconds =
		static_cast<double>(gc_standardizer_sort_profile::push_next_attempts_nanoseconds.load(std::memory_order_relaxed)) / 1'000'000'000.0;
	const double push_next_attempts_avg_seconds =
		push_next_attempts_total_seconds / static_cast<double>(iterations);
	const double final_graph_build_total_seconds =
		static_cast<double>(gc_standardizer_sort_profile::final_graph_build_nanoseconds.load(std::memory_order_relaxed)) / 1'000'000'000.0;
	const double final_graph_build_avg_seconds =
		final_graph_build_total_seconds / static_cast<double>(iterations);
	const double final_compare_total_seconds =
		static_cast<double>(gc_standardizer_sort_profile::final_compare_nanoseconds.load(std::memory_order_relaxed)) / 1'000'000'000.0;
	const double final_compare_avg_seconds =
		final_compare_total_seconds / static_cast<double>(iterations);
	const double labeling_only_avg_seconds =
		init_colors_avg_seconds + update_colors_avg_seconds + update_groups_avg_seconds + push_next_attempts_avg_seconds;
	const double sign_and_filter_avg_seconds = avg - labeling_only_avg_seconds;
	const double tracked_sign_and_filter_avg_seconds =
		final_graph_build_avg_seconds +
		final_compare_avg_seconds;
	const double other_sign_and_filter_avg_seconds =
		sign_and_filter_avg_seconds - tracked_sign_and_filter_avg_seconds;
	std::cout << "init_colors average = " << init_colors_avg_seconds << " s\n";
	std::cout << "init_colors share = "
	          << (100.0 * init_colors_avg_seconds / avg) << "%\n";
	std::cout << "update_colors share = "
	          << (100.0 * update_colors_avg_seconds / avg) << "%\n";
	std::cout << "update_groups share = "
	          << (100.0 * update_groups_avg_seconds / avg) << "%\n";
	std::cout << "push_next_attempts average = " << push_next_attempts_avg_seconds << " s\n";
	std::cout << "push_next_attempts share = "
	          << (100.0 * push_next_attempts_avg_seconds / avg) << "%\n";
	std::cout << "final_graph_build average = " << final_graph_build_avg_seconds << " s\n";
	std::cout << "final_graph_build share = "
	          << (100.0 * final_graph_build_avg_seconds / avg) << "%\n";
	std::cout << "final_compare average = " << final_compare_avg_seconds << " s\n";
	std::cout << "final_compare share = "
	          << (100.0 * final_compare_avg_seconds / avg) << "%\n";
	std::cout << "sign/filter share = "
	          << (100.0 * sign_and_filter_avg_seconds / avg) << "%\n";
	if (sign_and_filter_avg_seconds > 0.0) {
		std::cout << "  final_graph_build / sign-filter = "
		          << (100.0 * final_graph_build_avg_seconds / sign_and_filter_avg_seconds) << "%\n";
		std::cout << "  final_compare / sign-filter = "
		          << (100.0 * final_compare_avg_seconds / sign_and_filter_avg_seconds) << "%\n";
		std::cout << "  tracked / sign-filter = "
		          << (100.0 * tracked_sign_and_filter_avg_seconds / sign_and_filter_avg_seconds) << "%\n";
		std::cout << "  other / sign-filter = "
		          << (100.0 * other_sign_and_filter_avg_seconds / sign_and_filter_avg_seconds) << "%\n";
	}
	std::cout << "other sign/filter average = " << other_sign_and_filter_avg_seconds << " s\n";
	std::cout << "other sign/filter share = "
	          << (100.0 * other_sign_and_filter_avg_seconds / avg) << "%\n";
#endif

	return EXIT_SUCCESS;
}

void print_usage(const char* argv0) {
	std::cerr << "usage: " << argv0 << " [wheel_N] [split_contract_rounds] [repeat] [iterations] [workload]\n";
	std::cerr << "  workload: split-contract (default) or vmax\n";
	std::cerr << "  wheel_N must be one of 7,9,11,13,15,17,21,25,27,29,31,33,35\n";
}

} // namespace

int main(int argc, char** argv) {
	const int wheel_n = argc >= 2 ? std::stoi(argv[1]) : 11;
	const int expansion_rounds = argc >= 3 ? std::stoi(argv[2]) : 2;
	const int repeat = argc >= 4 ? std::stoi(argv[3]) : 1;
	const int iterations = argc >= 5 ? std::stoi(argv[4]) : 3;
	WorkloadKind workload_kind = WorkloadKind::SplitContract;
	try {
		workload_kind = argc >= 6 ? parse_workload_kind(argv[5]) : WorkloadKind::SplitContract;
	} catch (const std::invalid_argument& error) {
		std::cerr << error.what() << '\n';
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	if (expansion_rounds < 0 || repeat <= 0 || iterations <= 0) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	switch (wheel_n) {
	case 7:
		return run_standardize3_profile<7>(expansion_rounds, repeat, iterations, workload_kind);
	case 9:
		return run_standardize3_profile<9>(expansion_rounds, repeat, iterations, workload_kind);
	case 11:
		return run_standardize3_profile<11>(expansion_rounds, repeat, iterations, workload_kind);
	case 13:
		return run_standardize3_profile<13>(expansion_rounds, repeat, iterations, workload_kind);
	case 15:
		return run_standardize3_profile<15>(expansion_rounds, repeat, iterations, workload_kind);
	case 17:
		return run_standardize3_profile<17>(expansion_rounds, repeat, iterations, workload_kind);
	case 21:
		return run_standardize3_profile<21>(expansion_rounds, repeat, iterations, workload_kind);
	case 25:
		return run_standardize3_profile<25>(expansion_rounds, repeat, iterations, workload_kind);
	case 27:
		return run_standardize3_profile<27>(expansion_rounds, repeat, iterations, workload_kind);
	case 29:
		return run_standardize3_profile<29>(expansion_rounds, repeat, iterations, workload_kind);
	case 31:
		return run_standardize3_profile<31>(expansion_rounds, repeat, iterations, workload_kind);
	case 33:
		return run_standardize3_profile<33>(expansion_rounds, repeat, iterations, workload_kind);
	case 35:
		return run_standardize3_profile<35>(expansion_rounds, repeat, iterations, workload_kind);
	default:
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}
}
