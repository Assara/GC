#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>

#include "GraphGeneration/FinalCanonicalization.hpp"
#include "GraphGeneration/MappedFinalGraphFile.hpp"
#include "GraphGeneration/MappedGraphFile.hpp"
#include "graph.hpp"

#ifndef GC_TRICONNECTED_ORACLE_LOOP
#define GC_TRICONNECTED_ORACLE_LOOP 8
#endif

namespace {

constexpr int LOOP_NUMBER = GC_TRICONNECTED_ORACLE_LOOP;
constexpr int FIRST_VERTICES = 4;
constexpr int LAST_VERTICES = 2 * (LOOP_NUMBER - 1);

static_assert(LOOP_NUMBER >= 3);
static_assert(LOOP_NUMBER <= 8,
	"the regression oracle is deliberately capped at loop order eight");

template <int Vertices>
using graph_type = Graph<
	static_cast<Int>(Vertices),
	static_cast<Int>(LOOP_NUMBER + Vertices - 1),
	0,
	0,
	0,
	0,
	fieldType
>;

struct metadata {
	std::uint64_t automorphism_order = 0;
	GraphGeneration::final_graph_survival survival
		= GraphGeneration::final_graph_survival::none;

	bool operator==(const metadata&) const noexcept = default;
};

template <int Vertices>
using class_map = std::unordered_map<graph_type<Vertices>, metadata>;

struct counts {
	std::uint64_t unsigned_classes = 0;
	std::uint64_t odd_edges = 0;
	std::uint64_t even_edges_odd_vertices = 0;
};

template <int Vertices>
counts measure(const class_map<Vertices>& classes) {
	counts result;
	result.unsigned_classes = classes.size();
	for (const auto& [graph, info] : classes) {
		(void)graph;
		result.odd_edges += GraphGeneration::survives(
			info.survival,
			GraphGeneration::final_graph_survival::odd_edges
		);
		result.even_edges_odd_vertices += GraphGeneration::survives(
			info.survival,
			GraphGeneration::final_graph_survival::even_edges_odd_vertices
		);
	}
	return result;
}

template <int Vertices>
bool simple_adjacency(
	const graph_type<Vertices>& graph,
	std::array<std::uint64_t, Vertices>& adjacency
) noexcept {
	for (Int edge = 0; edge < graph_type<Vertices>::N_EDGES_; ++edge) {
		auto [first, second] = graph.getEdge(edge);
		if (first >= Vertices || second >= Vertices || first == second) {
			return false;
		}
		const std::uint64_t second_bit = std::uint64_t{1} << second;
		if ((adjacency[static_cast<std::size_t>(first)] & second_bit) != 0) {
			return false;
		}
		adjacency[static_cast<std::size_t>(first)] |= second_bit;
		adjacency[static_cast<std::size_t>(second)]
			|= std::uint64_t{1} << first;
	}
	return true;
}

template <int Vertices>
bool connected_after_removing(
	const std::array<std::uint64_t, Vertices>& adjacency,
	std::uint64_t removed
) noexcept {
	constexpr std::uint64_t all_vertices
		= (std::uint64_t{1} << Vertices) - 1;
	const std::uint64_t remaining = all_vertices & ~removed;
	if (remaining == 0) {
		return true;
	}
	std::uint64_t reached
		= std::uint64_t{1} << std::countr_zero(remaining);
	std::uint64_t frontier = reached;
	while (frontier != 0) {
		const unsigned vertex = std::countr_zero(frontier);
		frontier &= frontier - 1;
		const std::uint64_t discovered
			= adjacency[vertex] & remaining & ~reached;
		reached |= discovered;
		frontier |= discovered;
	}
	return reached == remaining;
}

template <int Vertices>
bool is_three_connected(const graph_type<Vertices>& graph) noexcept {
	static_assert(Vertices >= 4);
	std::array<std::uint64_t, Vertices> adjacency{};
	if (!simple_adjacency<Vertices>(graph, adjacency)) {
		return false;
	}
	for (const std::uint64_t neighbors : adjacency) {
		if (std::popcount(neighbors) < 3) {
			return false;
		}
	}
	if (!connected_after_removing<Vertices>(adjacency, 0)) {
		return false;
	}
	for (int first = 0; first < Vertices; ++first) {
		const std::uint64_t first_bit = std::uint64_t{1} << first;
		if (!connected_after_removing<Vertices>(adjacency, first_bit)) {
			return false;
		}
		for (int second = first + 1; second < Vertices; ++second) {
			if (!connected_after_removing<Vertices>(
					adjacency,
					first_bit | (std::uint64_t{1} << second)
				)) {
				return false;
			}
		}
	}
	return true;
}

template <int Vertices>
std::filesystem::path stage_path(
	const std::filesystem::path& directory,
	bool generated
) {
	std::string filename = "loop_" + std::to_string(LOOP_NUMBER)
		+ "_vertices_" + std::to_string(Vertices);
	if (generated) {
		filename += "_admissible";
	}
	filename += ".gcg";
	return directory / filename;
}

template <int Vertices>
class_map<Vertices> load_legacy(const std::filesystem::path& directory) {
	using graph = graph_type<Vertices>;
	const auto path = stage_path<Vertices>(directory, false);
	GraphGeneration::MappedGraphReader<graph> reader(path);
	GraphGeneration::final_graph_canonicalizer<graph> canonicalizer;
	class_map<Vertices> result;
	for (std::uint64_t index = 0; index < reader.size(); ++index) {
		const graph candidate = reader[index];
		if (!is_three_connected<Vertices>(candidate)) {
			continue;
		}
		const auto canonical = canonicalizer(candidate);
		const metadata info{
			.automorphism_order = canonical.automorphism_order,
			.survival = canonical.survival
		};
		const auto [position, inserted]
			= result.emplace(canonical.canonical_graph, info);
		if (!inserted && position->second != info) {
			throw std::logic_error(
				"legacy duplicates have inconsistent canonical metadata"
			);
		}
	}
	return result;
}

template <int Vertices>
class_map<Vertices> load_generated(const std::filesystem::path& directory) {
	using graph = graph_type<Vertices>;
	const auto path = stage_path<Vertices>(directory, true);
	GraphGeneration::MappedFinalGraphReader<graph> reader(path);
	class_map<Vertices> result;
	for (std::uint64_t index = 0; index < reader.size(); ++index) {
		const auto record = reader[index];
		if (!is_three_connected<Vertices>(record.graph)) {
			throw std::runtime_error(
				"generated final file contains a non-3-connected graph: "
				+ path.string()
			);
		}
		const metadata info{
			.automorphism_order = record.automorphism_order,
			.survival = record.survival
		};
		const auto [position, inserted] = result.emplace(record.graph, info);
		if (!inserted && position->second != info) {
			throw std::logic_error(
				"generated duplicates have inconsistent metadata"
			);
		}
	}
	return result;
}

template <typename GraphType>
void print_graph(const GraphType& graph) {
	for (Int edge = 0; edge < GraphType::N_EDGES_; ++edge) {
		const auto [first, second] = graph.getEdge(edge);
		if (edge != 0) {
			std::cout << ',';
		}
		std::cout << +first << '-' << +second;
	}
}

template <int Vertices>
bool report_stage(
	const std::filesystem::path* generated_directory,
	const std::filesystem::path& legacy_directory
) {
	const class_map<Vertices> expected = load_legacy<Vertices>(legacy_directory);
	const counts expected_counts = measure<Vertices>(expected);
	std::cout << "loop=" << LOOP_NUMBER
	          << " vertices=" << Vertices
	          << " edges=" << +graph_type<Vertices>::N_EDGES_
	          << " expected=" << expected_counts.unsigned_classes
	          << " odd=" << expected_counts.odd_edges
	          << " even=" << expected_counts.even_edges_odd_vertices;

	if (generated_directory == nullptr) {
		std::cout << '\n' << std::flush;
		return true;
	}

	const class_map<Vertices> actual
		= load_generated<Vertices>(*generated_directory);
	const counts actual_counts = measure<Vertices>(actual);
	std::uint64_t missing = 0;
	std::uint64_t extra = 0;
	std::uint64_t metadata_mismatches = 0;
	const graph_type<Vertices>* first_missing = nullptr;
	const graph_type<Vertices>* first_extra = nullptr;
	for (const auto& [graph, info] : expected) {
		const auto found = actual.find(graph);
		if (found == actual.end()) {
			++missing;
			if (first_missing == nullptr) {
				first_missing = &graph;
			}
		} else if (found->second != info) {
			++metadata_mismatches;
		}
	}
	for (const auto& [graph, info] : actual) {
		(void)info;
		if (!expected.contains(graph)) {
			++extra;
			if (first_extra == nullptr) {
				first_extra = &graph;
			}
		}
	}
	std::cout << " actual=" << actual_counts.unsigned_classes
	          << " actual_odd=" << actual_counts.odd_edges
	          << " actual_even=" << actual_counts.even_edges_odd_vertices
	          << " missing=" << missing
	          << " extra=" << extra
	          << " metadata_mismatches=" << metadata_mismatches
	          << '\n';
	if (first_missing != nullptr) {
		std::cout << "first_missing=";
		print_graph(*first_missing);
		std::cout << '\n';
	}
	if (first_extra != nullptr) {
		std::cout << "first_extra=";
		print_graph(*first_extra);
		std::cout << '\n';
	}
	std::cout << std::flush;
	return missing == 0 && extra == 0 && metadata_mismatches == 0;
}

template <int Vertices = FIRST_VERTICES>
bool report_all(
	const std::filesystem::path* generated_directory,
	const std::filesystem::path& legacy_directory,
	int last_vertices = LAST_VERTICES
) {
	bool success = report_stage<Vertices>(
		generated_directory, legacy_directory
	);
	if constexpr (Vertices < LAST_VERTICES) {
		if (Vertices < last_vertices) {
			success = report_all<Vertices + 1>(
				generated_directory, legacy_directory, last_vertices
			) && success;
		}
	}
	return success;
}

void usage(const char* executable) {
	std::cerr << "usage:\n"
	          << "  " << executable << " oracle <legacy_loop_directory>\n"
	          << "  " << executable
	          << " compare <generated_loop_directory> <legacy_loop_directory>"
	             " [last_vertices]\n";
}

} // namespace

int main(int argc, char** argv) {
	try {
		if (argc == 3 && std::string(argv[1]) == "oracle") {
			return report_all(nullptr, argv[2]) ? EXIT_SUCCESS : EXIT_FAILURE;
		}
		if ((argc == 4 || argc == 5) && std::string(argv[1]) == "compare") {
			const std::filesystem::path generated = argv[2];
			const int last_vertices = argc == 5
				? std::stoi(argv[4]) : LAST_VERTICES;
			if (last_vertices < FIRST_VERTICES
				|| last_vertices > LAST_VERTICES) {
				throw std::invalid_argument("last_vertices is out of range");
			}
			return report_all(&generated, argv[3], last_vertices)
				? EXIT_SUCCESS : EXIT_FAILURE;
		}
		usage(argv[0]);
		return EXIT_FAILURE;
	} catch (const std::exception& error) {
		std::cerr << "error: " << error.what() << '\n';
		return EXIT_FAILURE;
	}
}
