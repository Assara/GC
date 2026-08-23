#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

#include "GraphGeneration/FinalCanonicalization.hpp"
#include "GraphGeneration/MappedFinalGraphFile.hpp"
#include "graph.hpp"

#ifndef GC_GENERATION_LOOP
#define GC_GENERATION_LOOP 7
#endif

namespace {

constexpr int LOOP_NUMBER = GC_GENERATION_LOOP;
constexpr int MAX_VERTICES = 2 * (LOOP_NUMBER - 1);

template <int Vertices>
using GraphType = Graph<
	static_cast<Int>(Vertices),
	static_cast<Int>(LOOP_NUMBER + Vertices - 1),
	0,
	0,
	0,
	0,
	fieldType
>;

template <typename Graph>
class SequentialFinalWriter {
	public:
		explicit SequentialFinalWriter(const std::filesystem::path& path)
			: output_(path, std::ios::binary | std::ios::trunc) {
			if (!output_) {
				throw std::runtime_error("cannot open output file: " + path.string());
			}
			write_header();
		}

		void append(const GraphGeneration::final_canonicalization_result<Graph>& result) {
			output_.write(
				reinterpret_cast<const char*>(result.canonical_graph.half_edges.data()),
				static_cast<std::streamsize>(Graph::SIZE * sizeof(Int))
			);
			output_.write(
				reinterpret_cast<const char*>(&result.automorphism_order),
				sizeof(result.automorphism_order)
			);
			const auto survival = static_cast<std::uint8_t>(result.survival);
			output_.write(reinterpret_cast<const char*>(&survival), sizeof(survival));
			if (!output_) {
				throw std::runtime_error("failed writing final graph record");
			}
			++count_;
		}

		void finish() {
			output_.seekp(0);
			write_header();
			output_.flush();
			if (!output_) {
				throw std::runtime_error("failed finalizing output file");
			}
		}

		std::uint64_t size() const noexcept { return count_; }

	private:
		void write_header() {
			const GraphGeneration::MappedGraphFileHeader header{
				.vertices = Graph::N_VERTICES_,
				.edges = Graph::N_EDGES_,
				.graph_count = count_,
				.splittable_graph_count = 0,
				.format_version = GraphGeneration::MAPPED_GRAPH_FORMAT_VERSION,
				.record_size_bytes = GraphGeneration::FINAL_GRAPH_RECORD_SIZE<Graph>,
				.payload_kind = static_cast<std::uint64_t>(
					GraphGeneration::MappedGraphPayloadKind::final_graph_with_metadata
				),
				.flags = 0
			};
			output_.write(reinterpret_cast<const char*>(&header), sizeof(header));
		}

		std::ofstream output_;
		std::uint64_t count_ = 0;
};

template <typename Graph>
Graph decode_graph6(std::string_view encoded) {
	if (encoded.starts_with(">>graph6<<")) {
		encoded.remove_prefix(10);
	}
	if (encoded.empty()) {
		throw std::runtime_error("empty graph6 record");
	}
	const auto first = static_cast<unsigned char>(encoded.front());
	if (first == '~') {
		throw std::runtime_error("extended graph6 order is not supported");
	}
	const std::size_t vertices = first - 63U;
	if (vertices != Graph::N_VERTICES_) {
		throw std::runtime_error("graph6 vertex count does not match stage");
	}

	Graph graph{};
	std::size_t character = 1;
	unsigned int bits = 0;
	unsigned int remaining = 0;
	Int edge = 0;
	auto next_bit = [&]() {
		if (remaining == 0) {
			if (character >= encoded.size()) {
				throw std::runtime_error("truncated graph6 record");
			}
			const auto value = static_cast<unsigned char>(encoded[character++]);
			if (value < 63U || value > 126U) {
				throw std::runtime_error("invalid graph6 character");
			}
			bits = value - 63U;
			remaining = 6;
		}
		const bool result = (bits & (1U << (remaining - 1))) != 0;
		--remaining;
		return result;
	};

	for (Int second = 1; second < Graph::N_VERTICES_; ++second) {
		for (Int first_vertex = 0; first_vertex < second; ++first_vertex) {
			if (next_bit()) {
				if (edge >= Graph::N_EDGES_) {
					throw std::runtime_error("graph6 record has too many edges");
				}
				graph.setEdge(edge++, first_vertex, second);
			}
		}
	}
	if (edge != Graph::N_EDGES_) {
		throw std::runtime_error("graph6 record has wrong edge count");
	}
	return graph;
}

template <int Vertices>
int finalize_stage(const std::filesystem::path& output_path) {
	using StageGraph = GraphType<Vertices>;
	SequentialFinalWriter<StageGraph> writer(output_path);
	GraphGeneration::final_graph_canonicalizer<StageGraph> canonicalizer;
	std::uint64_t odd_edges = 0;
	std::uint64_t even_edges_odd_vertices = 0;
	std::string line;
	while (std::getline(std::cin, line)) {
		if (!line.empty() && line.back() == '\r') {
			line.pop_back();
		}
		if (line.empty()) {
			continue;
		}
		const auto result = canonicalizer(decode_graph6<StageGraph>(line));
		odd_edges += result.survives_odd_edges();
		even_edges_odd_vertices += result.survives_even_edges_odd_vertices();
		writer.append(result);
	}
	writer.finish();
	std::cout << writer.size() << '\t' << odd_edges << '\t'
	          << even_edges_odd_vertices << '\n';
	return 0;
}

template <int Candidate = 2>
int dispatch(int vertices, const std::filesystem::path& output_path) {
	if (vertices == Candidate) {
		return finalize_stage<Candidate>(output_path);
	}
	if constexpr (Candidate < MAX_VERTICES) {
		return dispatch<Candidate + 1>(vertices, output_path);
	}
	throw std::invalid_argument("unsupported vertex count");
}

} // namespace

int main(int argc, char** argv) {
	try {
		if (argc != 3) {
			std::cerr << "usage: " << argv[0] << " VERTICES OUTPUT.gcg\n";
			return 2;
		}
		return dispatch(std::stoi(argv[1]), argv[2]);
	} catch (const std::exception& error) {
		std::cerr << "geng finalization failed: " << error.what() << '\n';
		return 1;
	}
}
