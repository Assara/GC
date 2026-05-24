#include <chrono>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "examplegraphs.hpp"

namespace {

template <typename Fn>
double time_seconds(Fn&& fn) {
	const auto start = std::chrono::steady_clock::now();
	fn();
	const auto stop = std::chrono::steady_clock::now();
	return std::chrono::duration<double>(stop - start).count();
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
std::vector<typename OddGraphdegZero<N + 1>::SplitGraph> make_u_graphs() {
	using GraphType = typename OddGraphdegZero<N + 1>::SplitGraph;

	std::vector<GraphType> graphs;
	std::vector<bool> empty;
	graphs.push_back(U_graph<N>(empty));

	for (int length = 1; length <= (N - 5) / 2; ++length) {
		for (auto seq : left_right_sequences_starting_left(length)) {
			graphs.push_back(U_graph<N>(seq));
		}
	}

	return graphs;
}

template <typename GraphType>
std::size_t checksum_graph_basis(const BasisElement<GraphType, fieldType>& element) {
	std::size_t checksum = 0;
	for (const Int half_edge : element.getValue().half_edges) {
		checksum = checksum * 131 + static_cast<std::size_t>(half_edge);
	}
	checksum = checksum * 131 + static_cast<std::size_t>(element.getCoefficient().value());
	return checksum;
}

template <typename GraphType>
int run_graphs_perf(const std::vector<GraphType>& sources, int iterations, const char* label) {
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

	std::vector<double> standardize_seconds;
	std::vector<double> standardize2_seconds;
	standardize_seconds.reserve(iterations);
	standardize2_seconds.reserve(iterations);

	std::size_t standardize_checksum = 0;
	std::size_t standardize2_checksum = 0;

	for (int i = 0; i < iterations; ++i) {
		standardize_seconds.push_back(time_seconds([&]() {
			for (const GraphType& source : sources) {
				GraphType graph = source;
				const Element result = standardizer.standardize(graph, fieldType{1});
				standardize_checksum ^= checksum_graph_basis(result) + static_cast<std::size_t>(i);
			}
		}));

		standardize2_seconds.push_back(time_seconds([&]() {
			for (const GraphType& source : sources) {
				Element input(source, fieldType{1});
				const Element result = standardizer.standardize2(input);
				standardize2_checksum ^= checksum_graph_basis(result) + static_cast<std::size_t>(i);
			}
		}));
	}

	const double standardize_total = std::accumulate(
		standardize_seconds.begin(),
		standardize_seconds.end(),
		0.0
	);
	const double standardize2_total = std::accumulate(
		standardize2_seconds.begin(),
		standardize2_seconds.end(),
		0.0
	);
	const double standardize_avg = standardize_total / static_cast<double>(iterations);
	const double standardize2_avg = standardize2_total / static_cast<double>(iterations);

	std::cout << label << '\n';
	std::cout << "graphs = " << sources.size() << '\n';
	std::cout << "iterations = " << iterations << '\n';
	std::cout << "standardize total = " << standardize_total << " s\n";
	std::cout << "standardize avg = " << standardize_avg << " s\n";
	std::cout << "standardize2 total = " << standardize2_total << " s\n";
	std::cout << "standardize2 avg = " << standardize2_avg << " s\n";
	std::cout << "speedup standardize/standardize2 = "
	          << (standardize_avg / standardize2_avg) << "x\n";
	std::cout << "standardize checksum = " << standardize_checksum << '\n';
	std::cout << "standardize2 checksum = " << standardize2_checksum << '\n';

	return EXIT_SUCCESS;
}

template <Int N>
int run_wheel_perf(int iterations) {
	using GraphType = OddGraphdegZero<N + 1>;
	const std::vector<GraphType> sources{wheel_graph<N>()};
	const std::string label = "wheel = W" + std::to_string(N);
	return run_graphs_perf(sources, iterations, label.c_str());
}

template <Int N>
int run_u_perf(int iterations) {
	const auto sources = make_u_graphs<N>();
	const std::string label = "U_" + std::to_string(N) + " sequence sum";
	return run_graphs_perf(sources, iterations, label.c_str());
}

} // namespace

int main(int argc, char** argv) {
	const std::string mode = argc >= 2 ? argv[1] : "u15";
	const int iterations = argc >= 3 ? std::stoi(argv[2]) : 10;
	if (iterations <= 0) {
		std::cerr << "usage: " << argv[0] << " [u15|w15|w17|w19] [iterations]\n";
		return EXIT_FAILURE;
	}
	if (mode == "u15") {
		return run_u_perf<15>(iterations);
	}
	if (mode == "w15") {
		return run_wheel_perf<15>(iterations);
	}
	if (mode == "w17") {
		return run_wheel_perf<17>(iterations);
	}
	if (mode == "w19") {
		return run_wheel_perf<19>(iterations);
	}
	std::cerr << "mode must be one of u15, w15, w17, w19\n";
	return EXIT_FAILURE;
}
