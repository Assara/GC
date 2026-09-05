#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "GraphGeneration/MappedTransientGraph2File.hpp"
#include "GraphGeneration/TransientGraph2Standardizer.hpp"
#include "LinearProbeSet.hpp"

namespace GraphGeneration {

// Single-threaded, vertex-by-vertex generation of connected transient graphs.
// Bounds are compile-time because Graph encodes both dimensions in its type.
template <Int MaxLoopNumber, Int MaxVertices>
class GraphGenerationPipeline {
	static_assert(MaxLoopNumber >= 1, "the triangle seed has one loop");
	static_assert(MaxVertices >= 3);
	static_assert(sizeof(Int) == 1, "transient files store byte vertex labels");

	template <std::size_t V>
	static constexpr std::size_t max_loops =
		std::min<std::size_t>(MaxLoopNumber, (V - 1) * (V - 2) / 2);

	template <std::size_t V, std::size_t L>
	struct Partition {
		using graph_type = Graph<V, V - 1 + L, 0, 0, 0, 0, fieldType>;
		std::vector<graph_type> pending;
		linear_probe_set<graph_type> graphs;
		std::size_t candidates = 0;
	};

	template <std::size_t V, std::size_t... L>
	static auto stage_type(std::index_sequence<L...>) -> std::tuple<Partition<V, L>...>;
	template <std::size_t V>
	using Stage = decltype(stage_type<V>(std::make_index_sequence<max_loops<V> + 1>{}));
	template <std::size_t... I>
	static auto stages_type(std::index_sequence<I...>) -> std::tuple<Stage<I + 3>...>;
	using Stages = decltype(stages_type(std::make_index_sequence<MaxVertices - 2>{}));

	// Splitting only deposits candidates. Standardization and deduplication
	// happen when the pipeline drains these buffers after each parent.
	struct Collector {
		Stages& stages;

		template <typename G>
		void add(G child) {
			constexpr int v = G::N_VERTICES_;
			constexpr int loops = G::N_EDGES_ - v + 1;
			if constexpr (v >= 3 && v <= MaxVertices && loops >= 0
				&& loops <= MaxLoopNumber && G::N_EDGES_ <= v * (v - 1) / 2) {
				std::get<loops>(std::get<v - 3>(stages)).pending.push_back(std::move(child));
			} else {
				assert(false && "split emitted a graph outside the pipeline bounds");
			}
		}
	};

public:
	GraphGenerationPipeline() = default;
	GraphGenerationPipeline(const GraphGenerationPipeline&) = delete;
	GraphGenerationPipeline& operator=(const GraphGenerationPipeline&) = delete;
	GraphGenerationPipeline(GraphGenerationPipeline&&) = delete;
	GraphGenerationPipeline& operator=(GraphGenerationPipeline&&) = delete;

	struct StageSummary {
		Int vertices;
		std::size_t candidates;
		std::size_t unique_graphs;
	};

	static std::filesystem::path file_path(
		const std::filesystem::path& directory, Int vertices, Int edges
	) {
		return directory / ("transient_V" + std::to_string(vertices)
			+ "_E" + std::to_string(edges) + ".gcg");
	}

	// A run starts from K3 and stops at the explicit vertex bound. Files are
	// transient frontiers only; there is no final/admissible classification.
	std::vector<StageSummary> run(const std::filesystem::path& directory) {
		assert(!started_ && "pipeline instances can run only once");
		started_ = true;
		std::filesystem::create_directories(directory);
		// Check the whole output range before writing any stage.
		for (int v = 3; v <= MaxVertices; ++v) {
			const int maximum = std::min<int>(MaxLoopNumber, (v - 1) * (v - 2) / 2);
			for (int l = 0; l <= maximum; ++l) {
				const auto path = file_path(directory, v, v - 1 + l);
				if (std::filesystem::exists(path) || std::filesystem::exists(path.string() + ".tmp"))
					std::abort(); // Never overwrite existing output, including release builds.
			}
		}

		using Seed = typename Partition<3, 1>::graph_type;
		Seed seed;
		seed.setEdge(0, 0, 1);
		seed.setEdge(1, 0, 2);
		seed.setEdge(2, 1, 2);
		collector_.add(std::move(seed));
		drain<3>();

		std::vector<StageSummary> summaries;
		[&]<std::size_t... I>(std::index_sequence<I...>) {
			(process_stage<I + 3>(directory, summaries), ...);
		}(std::make_index_sequence<MaxVertices - 2>{});
		return summaries;
	}

private:
	template <typename G>
	static auto eligible_split_vertices(const G& graph,
		const std::array<Int, G::N_VERTICES_>& valences) {
		// Children reaching drain() are not yet ordered by valence.
		const Int maximum = *std::max_element(valences.begin(), valences.end());
		const auto low_valence_count = std::ranges::count_if(valences,
			[](Int valence) { return valence <= 2; });
		std::array<Int, G::N_VERTICES_> low_valence_neighbors{};
		if (low_valence_count != 0) {
			for (Int edge = 0; edge < G::N_EDGES_; ++edge) {
				const auto [a, b] = graph.getEdge(edge);
				low_valence_neighbors[a] += valences[b] <= 2;
				low_valence_neighbors[b] += valences[a] <= 2;
			}
		}
		std::array<bool, G::N_VERTICES_> eligible{};
		for (Int vertex = 0; vertex < G::N_VERTICES_; ++vertex) {
			// The splitting vertex itself need not be its own neighbour.
			eligible[vertex] = valences[vertex] == maximum
				&& low_valence_neighbors[vertex] == low_valence_count - (valences[vertex] <= 2);
		}
		return eligible;
	}

	template <std::size_t V>
	void drain() {
		std::apply([](auto&... partitions) {
			auto flush = [](auto& partition) {
				using G = typename std::remove_reference_t<decltype(partition)>::graph_type;
				transient_graph2_standardizer<G::N_VERTICES_, G::N_EDGES_> standardizer;
				partition.candidates += partition.pending.size();
				for (auto& graph : partition.pending) {
					const auto eligible = eligible_split_vertices(graph, graph.valence_array());
					if (std::ranges::none_of(eligible, [](bool value) { return value; })) continue;
					auto canonical = standardizer.standardize_no_sign(
						transient_graph2<G::N_VERTICES_, G::N_EDGES_>(std::move(graph)));
					const auto valences = canonical.graph().valence_array();
					assert(std::is_sorted(valences.begin(), valences.end()));
					partition.graphs.insert(canonical.graph());
				}
				partition.pending.clear();
			};
			(flush(partitions), ...);
		}, std::get<V - 3>(stages_));
	}

	template <typename PartitionType>
	static void write_partition(const std::filesystem::path& directory,
		const PartitionType& partition) {
		using G = typename PartitionType::graph_type;
		const auto path = file_path(directory, G::N_VERTICES_, G::N_EDGES_);
		const auto temporary = path.string() + ".tmp";
		{
			MappedTransientGraph2Writer<G> writer(temporary, partition.graphs.size());
			std::size_t record = 0;
			for (std::size_t slot = 0; slot < partition.graphs.capacity(); ++slot) {
				const auto& graph = partition.graphs.data()[slot];
				if (!graph.empty()) writer.write(record++, graph);
			}
		}
		std::filesystem::rename(temporary, path);
	}

	template <std::size_t V>
	void process_stage(const std::filesystem::path& directory,
		std::vector<StageSummary>& summaries) {
		StageSummary summary{static_cast<Int>(V), 0, 0};
		auto& stage = std::get<V - 3>(stages_);
		std::apply([&](auto&... partitions) {
			auto process = [&](auto& partition) {
				using G = typename std::remove_reference_t<decltype(partition)>::graph_type;
				summary.candidates += partition.candidates;
				summary.unique_graphs += partition.graphs.size();
				write_partition(directory, partition);
				if constexpr (V < MaxVertices) {
					for (std::size_t slot = 0; slot < partition.graphs.capacity(); ++slot) {
						const auto& graph = partition.graphs.data()[slot];
						if (graph.empty()) continue;
						const auto valences = graph.valence_array();
						const auto eligible = eligible_split_vertices(graph, valences);
						constexpr Int minimum = 1;
						const Int preserve = valences[V - 2];
						const transient_graph2<G::N_VERTICES_, G::N_EDGES_> parent(graph);
						for (Int vertex = 0; vertex < V; ++vertex) {
							if (eligible[vertex])
								parent.split(vertex, preserve, minimum, MaxLoopNumber, collector_);
						}
						drain<V + 1>();
					}
				}
				partition.graphs = {};
				decltype(partition.pending){}.swap(partition.pending);
			};
			(process(partitions), ...);
		}, stage);
		summaries.push_back(summary);
	}

	Stages stages_;
	Collector collector_{stages_};
	bool started_ = false;
};

} // namespace GraphGeneration
