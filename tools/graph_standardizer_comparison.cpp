#include <chrono>
#include <array>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <optional>
#include <string>
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
NautyLabeling<GraphType> nauty_labeling(const GraphType& input_graph, bool collect_automorphisms) {
	const int n = GraphType::N_VERTICES_;
	const int m = SETWORDSNEEDED(n);

	DYNALLSTAT(graph, g, g_sz);
	DYNALLSTAT(graph, cg, cg_sz);
	DYNALLOC2(graph, g, g_sz, m, n, "malloc");
	DYNALLOC2(graph, cg, cg_sz, m, n, "malloc");
	EMPTYGRAPH(g, m, n);

	for (Int e = 0; e < GraphType::N_EDGES_; ++e) {
		auto [u, v] = input_graph.getEdge(e);
		if (u != v) {
			ADDONEEDGE(g, static_cast<int>(u), static_cast<int>(v), m);
		}
	}

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
		g,
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

	DYNFREE(g, g_sz);
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
			>{}.standardize(graph_to_standardize, nauty.getCoefficient());
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
void nauty_canonicalize_once(const GraphType& input_graph) {
	(void)nauty_labeling(input_graph, false);
}

template <typename GraphType>
double benchmark_nauty(const std::vector<GraphType>& graphs, int iterations) {
	double total = 0.0;
	for (int i = 0; i < iterations; ++i) {
		total += time_seconds([&]() {
			for (const auto& graph : graphs) {
				nauty_canonicalize_once(graph);
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
int run_graph_standardizer_comparison(int repeat, int iterations) {
	auto base_graphs = make_wheel_homotopy_graphs<N>();
	std::vector<typename OddGraphdegZero<N + 1>::SplitGraph> graphs;
	graphs.reserve(base_graphs.size() * repeat);
	for (int r = 0; r < repeat; ++r) {
		graphs.insert(graphs.end(), base_graphs.begin(), base_graphs.end());
	}

	std::cout << "graph standardizer comparison\n";
	std::cout << "wheel = W" << +N << '\n';
	std::cout << "base homotopy graphs = " << base_graphs.size() << '\n';
	std::cout << "repeat = " << repeat << '\n';
	std::cout << "workload graphs = " << graphs.size() << '\n';
	std::cout << "iterations = " << iterations << '\n';

#if defined(GC_PROFILE_STANDARDIZER_SORT)
	gc_standardizer_sort_profile::reset();
#endif
	const double own_avg = benchmark_own_standardizer(graphs, iterations);
	std::cout << "own standardizer average = " << own_avg << " s\n";
#if defined(GC_PROFILE_STANDARDIZER_SORT)
	const double sort_total_seconds =
		static_cast<double>(gc_standardizer_sort_profile::insertion_sort_nanoseconds.load(std::memory_order_relaxed)) / 1'000'000'000.0;
	const double sort_avg_seconds = sort_total_seconds / static_cast<double>(iterations);
	const auto sort_calls = gc_standardizer_sort_profile::insertion_sort_calls.load(std::memory_order_relaxed);
	const auto sort_swaps = gc_standardizer_sort_profile::insertion_sort_swaps.load(std::memory_order_relaxed);
	std::cout << "own edge insertion sort average = " << sort_avg_seconds << " s\n";
	std::cout << "own edge insertion sort share = " << (100.0 * sort_avg_seconds / own_avg) << "%\n";
	std::cout << "own edge insertion sort calls = " << sort_calls << '\n';
	std::cout << "own edge insertion sort swaps = " << sort_swaps << '\n';

	const double color_update_total_seconds =
		static_cast<double>(gc_standardizer_sort_profile::vertex_color_update_nanoseconds.load(std::memory_order_relaxed)) / 1'000'000'000.0;
	const double color_update_avg_seconds = color_update_total_seconds / static_cast<double>(iterations);
	const auto color_update_calls = gc_standardizer_sort_profile::vertex_color_update_calls.load(std::memory_order_relaxed);
	std::cout << "own vertex color update average = " << color_update_avg_seconds << " s\n";
	std::cout << "own vertex color update share = " << (100.0 * color_update_avg_seconds / own_avg) << "%\n";
	std::cout << "own vertex color update calls = " << color_update_calls << '\n';

	const double bucket_init_total_seconds =
		static_cast<double>(gc_standardizer_sort_profile::vertex_bucket_init_nanoseconds.load(std::memory_order_relaxed)) / 1'000'000'000.0;
	const double bucket_init_avg_seconds = bucket_init_total_seconds / static_cast<double>(iterations);
	const auto bucket_init_calls = gc_standardizer_sort_profile::vertex_bucket_init_calls.load(std::memory_order_relaxed);
	std::cout << "own vertex bucket init average = " << bucket_init_avg_seconds << " s\n";
	std::cout << "own vertex bucket init share = " << (100.0 * bucket_init_avg_seconds / own_avg) << "%\n";
	std::cout << "own vertex bucket init calls = " << bucket_init_calls << '\n';
#endif

#if defined(GC_PERF_WITH_NAUTY)
	check_nauty_sign_correctness(graphs);
	const double nauty_avg = benchmark_nauty(graphs, iterations);
	std::cout << "nauty labeling average = " << nauty_avg << " s\n";
	const double nauty_full_avg = benchmark_nauty_standardizer(graphs, iterations);
	std::cout << "nauty full standardizer average = " << nauty_full_avg << " s\n";
#else
	std::cout << "nauty average = disabled; rebuild this tool with GC_PERF_WITH_NAUTY\n";
#endif

	return EXIT_SUCCESS;
}

void print_usage(const char* argv0) {
	std::cerr << "usage: " << argv0 << " [wheel_N] [repeat] [iterations]\n";
	std::cerr << "  wheel_N must be one of 9,13,15,17,21,25,29,33\n";
}

} // namespace

int main(int argc, char** argv) {
	const int wheel_n = argc >= 2 ? std::stoi(argv[1]) : 25;
	const int repeat = argc >= 3 ? std::stoi(argv[2]) : 1;
	const int iterations = argc >= 4 ? std::stoi(argv[3]) : 3;

	if (repeat <= 0 || iterations <= 0) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	switch (wheel_n) {
	case 9:
		return run_graph_standardizer_comparison<9>(repeat, iterations);
	case 13:
		return run_graph_standardizer_comparison<13>(repeat, iterations);
	case 15:
		return run_graph_standardizer_comparison<15>(repeat, iterations);
	case 17:
		return run_graph_standardizer_comparison<17>(repeat, iterations);
	case 21:
		return run_graph_standardizer_comparison<21>(repeat, iterations);
	case 25:
		return run_graph_standardizer_comparison<25>(repeat, iterations);
	case 29:
		return run_graph_standardizer_comparison<29>(repeat, iterations);
	case 33:
		return run_graph_standardizer_comparison<33>(repeat, iterations);
	default:
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}
}
