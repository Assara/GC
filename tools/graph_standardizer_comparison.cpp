#include <chrono>
#include <array>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

#if defined(GC_PERF_WITH_NAUTY)
#if __has_include(<nauty/nauty.h>)
#include <nauty/nauty.h>
#elif __has_include(<nauty.h>)
#include <nauty.h>
#else
#error "GC_PERF_WITH_NAUTY was set, but nauty headers were not found"
#endif
#endif

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
std::vector<typename OddGraphdegZero<N + 1>::SplitGraph> make_wheel_homotopy_graphs() {
	using WorkGraph = typename OddGraphdegZero<N + 1>::SplitGraph;

	std::vector<WorkGraph> graphs;
	std::vector<bool> empty;
	graphs.push_back(U_graph<N>(empty));

	for (int length = 1; length <= (N - 5) / 2; ++length) {
		for (auto seq : left_right_sequences_starting_left(length)) {
			graphs.push_back(U_graph<N>(seq));
		}
	}

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
double benchmark_own_standardizer(const std::vector<GraphType>& graphs, int iterations) {
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

	std::cout << "own checksum = " << checksum << '\n';
	return total / static_cast<double>(iterations);
}

#if defined(GC_PERF_WITH_NAUTY)
std::vector<std::array<Int, 255>>* current_automorphisms = nullptr;

void collect_automorphism_generator(int, int* perm, int*, int, int, int n) {
	if (current_automorphisms == nullptr) {
		return;
	}

	std::array<Int, 255> stored{};
	for (int i = 0; i < n; ++i) {
		stored[i] = static_cast<Int>(perm[i]);
	}
	current_automorphisms->push_back(stored);
}

template <typename GraphType>
struct NautyLabeling {
	std::array<Int, GraphType::N_VERTICES_> lab_new_to_old{};
	std::vector<std::array<Int, 255>> automorphism_generators;
};

template <typename GraphType>
struct PreparedNautyGraph {
	static constexpr int n = GraphType::N_VERTICES_;
	static constexpr int m = SETWORDSNEEDED(n);
	std::vector<graph> g = std::vector<graph>(static_cast<std::size_t>(m * n), 0);
};

template <typename GraphType>
NautyLabeling<GraphType> nauty_labeling(const PreparedNautyGraph<GraphType>& prepared_graph, bool collect_automorphisms);

template <typename GraphType>
PreparedNautyGraph<GraphType> prepare_nauty_graph(const GraphType& input_graph) {
	PreparedNautyGraph<GraphType> prepared;
	EMPTYGRAPH(prepared.g.data(), PreparedNautyGraph<GraphType>::m, PreparedNautyGraph<GraphType>::n);

	for (Int e = 0; e < GraphType::N_EDGES_; ++e) {
		auto [u, v] = input_graph.getEdge(e);
		if (u != v) {
			ADDONEEDGE(
				prepared.g.data(),
				static_cast<int>(u),
				static_cast<int>(v),
				PreparedNautyGraph<GraphType>::m
			);
		}
	}

	return prepared;
}

template <typename GraphType>
NautyLabeling<GraphType> nauty_labeling(const GraphType& input_graph, bool collect_automorphisms) {
	return nauty_labeling(prepare_nauty_graph(input_graph), collect_automorphisms);
}

template <typename GraphType>
NautyLabeling<GraphType> nauty_labeling(const PreparedNautyGraph<GraphType>& prepared_graph, bool collect_automorphisms) {
	const int n = PreparedNautyGraph<GraphType>::n;
	const int m = PreparedNautyGraph<GraphType>::m;

	DYNALLSTAT(graph, cg, cg_sz);
	DYNALLOC2(graph, cg, cg_sz, m, n, "malloc");

	std::vector<int> lab(n);
	std::vector<int> ptn(n);
	std::vector<int> orbits(n);
	DEFAULTOPTIONS_GRAPH(options);
	statsblk stats;
	options.getcanon = TRUE;

	NautyLabeling<GraphType> result;
	if (collect_automorphisms) {
		current_automorphisms = &result.automorphism_generators;
		options.userautomproc = collect_automorphism_generator;
	}

	densenauty(
		const_cast<graph*>(prepared_graph.g.data()),
		lab.data(),
		ptn.data(),
		orbits.data(),
		&options,
		&stats,
		m,
		n,
		cg
	);

	if (collect_automorphisms) {
		current_automorphisms = nullptr;
	}

	for (Int i = 0; i < GraphType::N_VERTICES_; ++i) {
		result.lab_new_to_old[i] = static_cast<Int>(lab[i]);
	}

	DYNFREE(cg, cg_sz);
	return result;
}

template <typename GraphType>
std::array<Int, GraphType::N_VERTICES_> inverse_vertex_permutation(
	const std::array<Int, GraphType::N_VERTICES_>& perm
) {
	std::array<Int, GraphType::N_VERTICES_> inverse{};
	for (Int i = 0; i < GraphType::N_VERTICES_; ++i) {
		inverse[perm[i]] = i;
	}
	return inverse;
}

template <typename GraphType>
BasisElement<GraphType, fieldType> vertex_permute_and_normalize(
	const GraphType& input_graph,
	const std::array<Int, GraphType::N_VERTICES_>& old_to_new
) {
	GraphType graph = input_graph;
	signedInt sign = graph.permuteVertices(Permutation<GraphType::N_VERTICES_>(old_to_new));
	sign *= graph.directAndSortEdges();

	if constexpr (GraphType::SWAP_EDGE_SIGN == -1) {
		if (graph.has_double_edge()) {
			return BasisElement<GraphType, fieldType>(graph, fieldType{0});
		}
	}

	return BasisElement<GraphType, fieldType>(graph, fieldType{sign});
}

template <typename GraphType>
signedInt field_sign(fieldType value) {
	if (value == fieldType{0}) {
		return 0;
	}
	if (value == fieldType{1}) {
		return 1;
	}
	if (value == fieldType{-1}) {
		return -1;
	}
	return 2;
}

template <typename GraphType>
std::optional<signedInt> automorphism_generator_sign(
	const GraphType& input_graph,
	const std::array<Int, 255>& generator
) {
	std::array<Int, GraphType::N_VERTICES_> old_to_new{};
	for (Int i = 0; i < GraphType::N_VERTICES_; ++i) {
		old_to_new[i] = generator[i];
	}

	std::array<Int, GraphType::N_VERTICES_> identity{};
	for (Int i = 0; i < GraphType::N_VERTICES_; ++i) {
		identity[i] = i;
	}
	const auto base = vertex_permute_and_normalize(input_graph, identity);

	auto moved = vertex_permute_and_normalize(input_graph, old_to_new);
	if (moved.getValue() == base.getValue()) {
		return field_sign<GraphType>(moved.getCoefficient() * base.getCoefficient());
	}

	auto inverse = inverse_vertex_permutation<GraphType>(old_to_new);
	moved = vertex_permute_and_normalize(input_graph, inverse);
	if (moved.getValue() == base.getValue()) {
		return field_sign<GraphType>(moved.getCoefficient() * base.getCoefficient());
	}

	return std::nullopt;
}

template <typename GraphType>
BasisElement<GraphType, fieldType> nauty_standardize_for_sign_test(const GraphType& input_graph) {
	const auto labeling = nauty_labeling(input_graph, true);
	const auto old_to_new = inverse_vertex_permutation<GraphType>(labeling.lab_new_to_old);
	auto nauty_form = vertex_permute_and_normalize(input_graph, old_to_new);

	for (const auto& generator : labeling.automorphism_generators) {
		const auto sign = automorphism_generator_sign(input_graph, generator);
		if (sign.has_value() && *sign < 0) {
			nauty_form = BasisElement<GraphType, fieldType>(nauty_form.getValue(), fieldType{0});
			break;
		}
	}

	return nauty_form;
}

template <typename GraphType>
void check_nauty_sign_correctness(const std::vector<GraphType>& graphs) {
	using Element = BasisElement<GraphType, fieldType>;

	std::size_t canonical_mismatches = 0;
	std::size_t sign_mismatches = 0;
	std::size_t zero_mismatches = 0;
	std::size_t unrecognized_signs = 0;

	std::size_t graph_index = 0;
	for (const auto& graph : graphs) {
		Element own_input(graph, fieldType{1});
		const auto own = GraphType::canonized(own_input);
		auto nauty = nauty_standardize_for_sign_test(graph);

		const signedInt nauty_to_own_sign = [&]() -> signedInt {
			if (nauty.getCoefficient() == fieldType{0}) {
				return 0;
			}
			auto graph_to_standardize = nauty.getValue();
			auto adjusted = GraphStandardizer<
				GraphType::N_VERTICES_,
				GraphType::N_EDGES_,
				GraphType::N_OUT_HAIR_,
				GraphType::N_IN_HAIR_,
				GraphType::C_,
				GraphType::D_,
				fieldType
			>{}.lexicographical_standardize(graph_to_standardize, nauty.getCoefficient());
			if (adjusted.getValue() != own.getValue()) {
				++canonical_mismatches;
			}
			return field_sign<GraphType>(adjusted.getCoefficient());
		}();

		const signedInt own_sign = field_sign<GraphType>(own.getCoefficient());
		if (own_sign == 2 || nauty_to_own_sign == 2) {
			++unrecognized_signs;
			continue;
		}
		if ((own_sign == 0) != (nauty_to_own_sign == 0)) {
			++zero_mismatches;
		}
		if (own_sign != nauty_to_own_sign) {
			if (sign_mismatches == 0) {
				std::cout << "first nauty sign mismatch index = " << graph_index
				          << ", own = " << own_sign
				          << ", nauty = " << nauty_to_own_sign << '\n';
			}
			++sign_mismatches;
		}
		++graph_index;
	}

	std::cout << "nauty sign correctness checked = " << graphs.size() << '\n';
	std::cout << "nauty sign mismatches = " << sign_mismatches << '\n';
	std::cout << "nauty zero mismatches = " << zero_mismatches << '\n';
	std::cout << "nauty canonical mismatches after adjustment = " << canonical_mismatches << '\n';
	std::cout << "nauty unrecognized signs = " << unrecognized_signs << '\n';
}

template <typename GraphType>
void nauty_canonicalize_once(const PreparedNautyGraph<GraphType>& prepared_graph) {
	(void)nauty_labeling(prepared_graph, false);
}

template <typename GraphType>
double benchmark_nauty(const std::vector<GraphType>& graphs, int iterations) {
	std::vector<PreparedNautyGraph<GraphType>> prepared_graphs;
	prepared_graphs.reserve(graphs.size());
	for (const auto& graph : graphs) {
		prepared_graphs.push_back(prepare_nauty_graph(graph));
	}

	double total = 0.0;
	for (int i = 0; i < iterations; ++i) {
		total += time_seconds([&]() {
			for (const auto& prepared_graph : prepared_graphs) {
				nauty_canonicalize_once<GraphType>(prepared_graph);
			}
		});
	}
	return total / static_cast<double>(iterations);
}

template <typename GraphType>
double benchmark_nauty_standardizer(const std::vector<GraphType>& graphs, int iterations) {
	double total = 0.0;
	std::size_t checksum = 0;

	for (int i = 0; i < iterations; ++i) {
		total += time_seconds([&]() {
			for (const auto& graph : graphs) {
				auto canon = nauty_standardize_for_sign_test(graph);
				checksum += canon.getValue().half_edges[0];
			}
		});
	}

	std::cout << "nauty full checksum = " << checksum << '\n';
	return total / static_cast<double>(iterations);
}
#endif

template <Int N>
int run_graph_standardizer_comparison(
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

	std::cout << "graph standardizer comparison\n";
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
	const double own_avg = benchmark_own_standardizer(graphs, iterations);
	std::cout << "own standardizer average = " << own_avg << " s\n";
#if defined(GC_PROFILE_STANDARDIZER_SORT)
	const double create_final_attempts_avg_seconds =
		static_cast<double>(gc_standardizer_sort_profile::create_final_attempts_nanoseconds.load(std::memory_order_relaxed))
		/ (1'000'000'000.0 * static_cast<double>(iterations));
	const double sort_and_filter_avg_seconds =
		static_cast<double>(gc_standardizer_sort_profile::sort_and_filter_nanoseconds.load(std::memory_order_relaxed))
		/ (1'000'000'000.0 * static_cast<double>(iterations));
	std::cout << "create final attempts share = "
	          << (100.0 * create_final_attempts_avg_seconds / own_avg) << "%\n";
	std::cout << "sort/filter share = "
	          << (100.0 * sort_and_filter_avg_seconds / own_avg) << "%\n";
#endif

#if defined(GC_PERF_WITH_NAUTY)
	const double nauty_avg = benchmark_nauty(graphs, iterations);
	std::cout << "nauty labeling average = " << nauty_avg << " s\n";
#else
	std::cout << "nauty average = disabled; rebuild this tool with GC_PERF_WITH_NAUTY\n";
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
		return run_graph_standardizer_comparison<7>(expansion_rounds, repeat, iterations, workload_kind);
	case 9:
		return run_graph_standardizer_comparison<9>(expansion_rounds, repeat, iterations, workload_kind);
	case 11:
		return run_graph_standardizer_comparison<11>(expansion_rounds, repeat, iterations, workload_kind);
	case 13:
		return run_graph_standardizer_comparison<13>(expansion_rounds, repeat, iterations, workload_kind);
	case 15:
		return run_graph_standardizer_comparison<15>(expansion_rounds, repeat, iterations, workload_kind);
	case 17:
		return run_graph_standardizer_comparison<17>(expansion_rounds, repeat, iterations, workload_kind);
	case 21:
		return run_graph_standardizer_comparison<21>(expansion_rounds, repeat, iterations, workload_kind);
	case 25:
		return run_graph_standardizer_comparison<25>(expansion_rounds, repeat, iterations, workload_kind);
	case 27:
		return run_graph_standardizer_comparison<27>(expansion_rounds, repeat, iterations, workload_kind);
	case 29:
		return run_graph_standardizer_comparison<29>(expansion_rounds, repeat, iterations, workload_kind);
	case 31:
		return run_graph_standardizer_comparison<31>(expansion_rounds, repeat, iterations, workload_kind);
	case 33:
		return run_graph_standardizer_comparison<33>(expansion_rounds, repeat, iterations, workload_kind);
	case 35:
		return run_graph_standardizer_comparison<35>(expansion_rounds, repeat, iterations, workload_kind);
	default:
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}
}
