#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "examplegraphs.hpp"

namespace {

int max_threads() {
#ifdef _OPENMP
	return omp_get_max_threads();
#else
	return 1;
#endif
}

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
std::vector<BasisElement<typename OddGraphdegZero<N + 1>::SplitGraph, fieldType>>
make_wheel_homotopy_terms() {
	using WorkGraph = typename OddGraphdegZero<N + 1>::SplitGraph;
	using Element = BasisElement<WorkGraph, fieldType>;

	std::vector<Element> terms;

	std::vector<bool> empty;
	terms.emplace_back(U_graph<N>(empty), fieldType{1});

	for (int length = 1; length <= (N - 5) / 2; ++length) {
		for (auto seq : left_right_sequences_starting_left(length)) {
			terms.emplace_back(U_graph<N>(seq), fieldType{2});
		}
	}

	return terms;
}

template <Int N>
int run_standardize_perf(int repeat, int iterations) {
	using WorkGraph = typename OddGraphdegZero<N + 1>::SplitGraph;
	using Element = BasisElement<WorkGraph, fieldType>;

	auto homotopy_terms = make_wheel_homotopy_terms<N>();
	if (repeat == 0) {
		constexpr std::size_t target_terms = 10000;
		repeat = static_cast<int>(
			(target_terms + homotopy_terms.size() - 1) / homotopy_terms.size()
		);
	}

	std::vector<Element> workload;
	workload.reserve(homotopy_terms.size() * repeat);
	for (int r = 0; r < repeat; ++r) {
		workload.insert(
			workload.end(),
			homotopy_terms.begin(),
			homotopy_terms.end()
		);
	}

	std::cout << "standardize_all perf: wheel homotopy U_N\n";
	std::cout << "wheel = W" << +N << '\n';
	std::cout << "homotopy terms = " << homotopy_terms.size() << '\n';
	std::cout << "repeat = " << repeat << '\n';
	std::cout << "workload terms = " << workload.size() << '\n';
	std::cout << "iterations = " << iterations << '\n';
	std::cout << "omp max threads = " << max_threads() << '\n';

	std::vector<double> seconds;
	seconds.reserve(iterations);
	bigInt output_size = 0;

	for (int i = 0; i < iterations; ++i) {
		VectorSpace::LinComb<WorkGraph, fieldType> lc(
			std::vector<Element>(workload),
			AssumeBasisOrderTag{}
		);

		const auto start = std::chrono::steady_clock::now();
		lc.standardize_all();
		const auto stop = std::chrono::steady_clock::now();

		const double elapsed =
			std::chrono::duration<double>(stop - start).count();
		seconds.push_back(elapsed);
		output_size = lc.size();

		std::cout << "iteration " << (i + 1)
		          << ": " << elapsed << " s"
		          << ", output terms = " << output_size << '\n';
	}

	const double total = std::accumulate(seconds.begin(), seconds.end(), 0.0);
	std::cout << "average = " << (total / seconds.size()) << " s\n";
	return EXIT_SUCCESS;
}

void print_usage(const char* argv0) {
	std::cerr << "usage: " << argv0 << " [wheel_N] [repeat] [iterations]\n";
	std::cerr << "  wheel_N must be one of 5,7,9,11,13,15,17,19,21,23,25,27,29,31,33\n";
	std::cerr << "  repeat = 0 means repeat enough homotopies to reach at least 10000 terms\n";
}

} // namespace

int main(int argc, char** argv) {
	const int wheel_n = argc >= 2 ? std::stoi(argv[1]) : 33;
	const int repeat = argc >= 3 ? std::stoi(argv[2]) : 1;
	const int iterations = argc >= 4 ? std::stoi(argv[3]) : 3;

	if (repeat < 0 || iterations <= 0) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	switch (wheel_n) {
	case 5:
		return run_standardize_perf<5>(repeat, iterations);
	case 7:
		return run_standardize_perf<7>(repeat, iterations);
	case 9:
		return run_standardize_perf<9>(repeat, iterations);
	case 11:
		return run_standardize_perf<11>(repeat, iterations);
	case 13:
		return run_standardize_perf<13>(repeat, iterations);
	case 15:
		return run_standardize_perf<15>(repeat, iterations);
	case 17:
		return run_standardize_perf<17>(repeat, iterations);
	case 19:
		return run_standardize_perf<19>(repeat, iterations);
	case 21:
		return run_standardize_perf<21>(repeat, iterations);
	case 23:
		return run_standardize_perf<23>(repeat, iterations);
	case 25:
		return run_standardize_perf<25>(repeat, iterations);
	case 27:
		return run_standardize_perf<27>(repeat, iterations);
	case 29:
		return run_standardize_perf<29>(repeat, iterations);
	case 31:
		return run_standardize_perf<31>(repeat, iterations);
	case 33:
		return run_standardize_perf<33>(repeat, iterations);
	default:
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}
}
