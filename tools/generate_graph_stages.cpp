#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include <omp.h>

#include "GraphGeneration/FinalCanonicalization.hpp"
#include "GraphGeneration/MappedFinalGraphFile.hpp"
#include "GraphGeneration/MappedGraphFile.hpp"
#include "GraphGeneration/MappedSupportTransientFile.hpp"
#include "GraphGeneration/UnrootedSupportTransientGraph.hpp"
#include "GraphGeneration/UnrootedSupportTransientStandardizer.hpp"
#include "LinearProbeSet.hpp"
#include "graph.hpp"

#ifndef GC_GENERATION_LOOP
#define GC_GENERATION_LOOP 7
#endif

namespace {

constexpr int LOOP_NUMBER = GC_GENERATION_LOOP;
constexpr int SEED_VERTICES = 2;
constexpr int MAX_VERTICES = 2 * (LOOP_NUMBER - 1);
constexpr const char* DIMENSION_LOG_FILENAME = "generation_dimensions.tsv";
constexpr const char* TRANSIENT_LOG_KIND
	= "ordered_unrooted_support_leaf_transient";

static_assert(
	LOOP_NUMBER >= 3,
	"the unrooted support generator requires loop number >= 3"
);
static_assert(3 * LOOP_NUMBER - 3 <= std::numeric_limits<Int>::max());
static_assert(MAX_VERTICES < 64,
	"stage dispatch materializes a child support with at most 64 vertices");

template <int Vertices>
using AdmissibleGraph = Graph<
	static_cast<Int>(Vertices),
	static_cast<Int>(LOOP_NUMBER + Vertices - 1),
	0,
	0,
	0,
	0,
	fieldType
>;

template <int Vertices>
using Transient = GraphGeneration::support_transient_graph<
	static_cast<Int>(Vertices),
	static_cast<Int>(LOOP_NUMBER + Vertices - 1),
	fieldType
>;

template <typename GraphType>
using Standardizer = GraphStandardizer<
	GraphType::N_VERTICES_,
	GraphType::N_EDGES_,
	GraphType::N_OUT_HAIR_,
	GraphType::N_IN_HAIR_,
	GraphType::C_,
	GraphType::D_,
	fieldType
>;

struct ProfileTimer {
	std::uint64_t calls = 0;
	std::uint64_t nanoseconds = 0;

	void merge(const ProfileTimer& other) noexcept {
		calls += other.calls;
		nanoseconds += other.nanoseconds;
	}
};

struct GenerationStageProfile {
	ProfileTimer transient_canonicalization;
	ProfileTimer admissible_canonicalization;
	ProfileTimer transient_map_insertion;
	ProfileTimer admissible_map_insertion;
	std::uint64_t relevant_children = 0;
	std::uint64_t workers = 0;

	void merge(const GenerationStageProfile& other) noexcept {
		transient_canonicalization.merge(other.transient_canonicalization);
		admissible_canonicalization.merge(other.admissible_canonicalization);
		transient_map_insertion.merge(other.transient_map_insertion);
		admissible_map_insertion.merge(other.admissible_map_insertion);
		relevant_children += other.relevant_children;
		workers += other.workers;
	}
};

template <typename Function>
decltype(auto) profiled_call(ProfileTimer& timer, Function&& function) {
#if defined(GC_PROFILE_GRAPH_GENERATION)
	using result_type = std::invoke_result_t<Function>;
	++timer.calls;
	const auto start = std::chrono::steady_clock::now();
	if constexpr (std::is_void_v<result_type>) {
		std::invoke(std::forward<Function>(function));
		timer.nanoseconds += static_cast<std::uint64_t>(
			std::chrono::duration_cast<std::chrono::nanoseconds>(
				std::chrono::steady_clock::now() - start
			).count()
		);
	} else {
		result_type result = std::invoke(std::forward<Function>(function));
		timer.nanoseconds += static_cast<std::uint64_t>(
			std::chrono::duration_cast<std::chrono::nanoseconds>(
				std::chrono::steady_clock::now() - start
			).count()
		);
		return result;
	}
#else
	(void)timer;
	return std::invoke(std::forward<Function>(function));
#endif
}

#if defined(GC_PROFILE_GRAPH_GENERATION)
double profile_seconds(const ProfileTimer& timer) noexcept {
	return static_cast<double>(timer.nanoseconds) / 1'000'000'000.0;
}

void print_profile_timer(const char* name, const ProfileTimer& timer) {
	std::cout << "profile timer=" << name
	          << " calls=" << timer.calls
	          << " summed_worker_seconds=" << profile_seconds(timer) << '\n';
}

void print_generation_profile(
	int vertices,
	int edges,
	const GenerationStageProfile& profile,
	double generation_wall_seconds,
	double transient_map_write_wall_seconds,
	double admissible_map_write_wall_seconds,
	double admissible_map_scan_wall_seconds
) {
	std::cout << "profile stage"
	          << " loop=" << LOOP_NUMBER
	          << " vertices=" << vertices
	          << " edges=" << edges
	          << " workers=" << profile.workers
	          << " relevant_children=" << profile.relevant_children
	          << " generation_wall_seconds=" << generation_wall_seconds
	          << '\n';
	print_profile_timer(
		"transient_canonicalization", profile.transient_canonicalization
	);
	print_profile_timer(
		"admissible_canonicalization", profile.admissible_canonicalization
	);
	print_profile_timer("transient_map_insertion", profile.transient_map_insertion);
	print_profile_timer(
		"admissible_map_insertion", profile.admissible_map_insertion
	);
	std::cout << "profile timer=transient_map_traversal_and_write"
	          << " wall_seconds=" << transient_map_write_wall_seconds << '\n';
	std::cout << "profile timer=admissible_map_traversal_and_write"
	          << " wall_seconds=" << admissible_map_write_wall_seconds << '\n';
	std::cout << "profile timer=admissible_map_metadata_scan"
	          << " wall_seconds=" << admissible_map_scan_wall_seconds << '\n'
	          << std::flush;
}
#endif

template <typename GraphType>
class ConcurrentGraphSet {
	public:
		static constexpr std::size_t SHARD_COUNT = 256;
		static_assert((SHARD_COUNT & (SHARD_COUNT - 1)) == 0);

		void insert(GraphType graph) {
			const std::size_t shard_index = std::hash<GraphType>{}(graph)
				& (SHARD_COUNT - 1);
			Shard& shard = shards_[shard_index];
			std::lock_guard lock(shard.mutex);
			const bool inserted = shard.graphs.insert(std::move(graph)).second;
			shard.duplicates += !inserted;
		}

		std::uint64_t size() const noexcept {
			std::uint64_t count = 0;
			for (const Shard& shard : shards_) {
				count += shard.graphs.size();
			}
			return count;
		}

		std::uint64_t duplicates() const noexcept {
			std::uint64_t count = 0;
			for (const Shard& shard : shards_) {
				count += shard.duplicates;
			}
			return count;
		}

		const std::unordered_set<GraphType>& shard(std::size_t index) const {
			return shards_[index].graphs;
		}

	private:
		struct alignas(64) Shard {
			std::mutex mutex;
			std::unordered_set<GraphType> graphs;
			std::uint64_t duplicates = 0;
		};

		std::array<Shard, SHARD_COUNT> shards_;
};

template <typename SupportGraph>
class PartitionedTransientGraphSet {
	public:
		using set_type = linear_probe_set<SupportGraph>;
		static constexpr std::size_t PARTITION_SIDE
			= static_cast<std::size_t>(SupportGraph::N_VERTICES_) + 1;
		static constexpr std::size_t PARTITION_COUNT
			= PARTITION_SIDE * PARTITION_SIDE;
		static constexpr std::size_t SHARDS_PER_PARTITION = 16;
		static constexpr std::size_t SLOT_COUNT
			= PARTITION_COUNT * SHARDS_PER_PARTITION;

		void insert(SupportGraph graph) {
			const auto adjacency = graph.support_adjacency();
			Int maximum = 0;
			Int tied = 0;
			for (Int vertex = 0; vertex < SupportGraph::N_VERTICES_; ++vertex) {
				const Int degree = static_cast<Int>(std::popcount(
					adjacency[static_cast<std::size_t>(vertex)]
				));
				if (degree > maximum) {
					maximum = degree;
					tied = 1;
				} else if (degree == maximum) {
					++tied;
				}
			}
			insert(std::move(graph), maximum, tied);
		}

		void insert(SupportGraph graph, Int maximum, Int tied) {
			const std::size_t partition = partition_index(maximum, tied);
			const std::size_t hash = std::hash<SupportGraph>{}(graph);
			const std::size_t slot = partition * SHARDS_PER_PARTITION
				+ (hash & (SHARDS_PER_PARTITION - 1));
			Shard& shard = shards_[slot];
			std::lock_guard lock(shard.mutex);
			const bool inserted = shard.graphs.insert(std::move(graph));
			shard.duplicates += !inserted;
		}

		std::uint64_t size() const noexcept {
			std::uint64_t result = 0;
			for (const Shard& shard : shards_) {
				result += shard.graphs.size();
			}
			return result;
		}

		std::uint64_t duplicates() const noexcept {
			std::uint64_t result = 0;
			for (const Shard& shard : shards_) result += shard.duplicates;
			return result;
		}

		const set_type& slot(std::size_t index) const noexcept {
			return shards_[index].graphs;
		}

		std::uint64_t partition_size(Int maximum, Int tied) const noexcept {
			std::uint64_t result = 0;
			const std::size_t begin
				= partition_index(maximum, tied) * SHARDS_PER_PARTITION;
			for (std::size_t shard = 0; shard < SHARDS_PER_PARTITION; ++shard) {
				result += shards_[begin + shard].graphs.size();
			}
			return result;
		}

	private:
		static_assert((SHARDS_PER_PARTITION & (SHARDS_PER_PARTITION - 1)) == 0);

		static std::size_t partition_index(Int maximum, Int tied) noexcept {
			return static_cast<std::size_t>(maximum) * PARTITION_SIDE
				+ static_cast<std::size_t>(tied);
		}

		struct alignas(64) Shard {
			std::mutex mutex;
			set_type graphs;
			std::uint64_t duplicates = 0;
		};

		std::array<Shard, SLOT_COUNT> shards_;
};

template <typename GraphType>
class ConcurrentFinalGraphSet {
	public:
		using CanonicalizationResult
			= GraphGeneration::final_canonicalization_result<GraphType>;

		struct metadata_type {
			std::uint64_t automorphism_order = 0;
			GraphGeneration::final_graph_survival survival
				= GraphGeneration::final_graph_survival::none;
		};

		static constexpr std::size_t SHARD_COUNT = 256;
		static_assert((SHARD_COUNT & (SHARD_COUNT - 1)) == 0);

		void insert(CanonicalizationResult result) {
			const std::size_t shard_index
				= std::hash<GraphType>{}(result.canonical_graph) & (SHARD_COUNT - 1);
			Shard& shard = shards_[shard_index];
			std::lock_guard lock(shard.mutex);
			const metadata_type incoming{
				.automorphism_order = result.automorphism_order,
				.survival = result.survival
			};
			auto [iterator, inserted] = shard.graphs.try_emplace(
				std::move(result.canonical_graph), incoming
			);
			if (!inserted) {
				++shard.duplicates;
				if (iterator->second.automorphism_order != incoming.automorphism_order
					|| iterator->second.survival != incoming.survival) {
					metadata_conflict_.store(true, std::memory_order_relaxed);
				}
			}
		}

		std::uint64_t size() const noexcept {
			std::uint64_t count = 0;
			for (const Shard& shard : shards_) {
				count += shard.graphs.size();
			}
			return count;
		}

		std::uint64_t duplicates() const noexcept {
			std::uint64_t count = 0;
			for (const Shard& shard : shards_) {
				count += shard.duplicates;
			}
			return count;
		}

		std::uint64_t surviving(
			GraphGeneration::final_graph_survival sign_case
		) const noexcept {
			std::uint64_t count = 0;
			for (const Shard& shard : shards_) {
				for (const auto& [graph, metadata] : shard.graphs) {
					(void)graph;
					count += GraphGeneration::survives(metadata.survival, sign_case);
				}
			}
			return count;
		}

		bool has_metadata_conflict() const noexcept {
			return metadata_conflict_.load(std::memory_order_relaxed);
		}

		const std::unordered_map<GraphType, metadata_type>& shard(
			std::size_t index
		) const {
			return shards_[index].graphs;
		}

	private:
		struct alignas(64) Shard {
			std::mutex mutex;
			std::unordered_map<GraphType, metadata_type> graphs;
			std::uint64_t duplicates = 0;
		};

		std::array<Shard, SHARD_COUNT> shards_;
		std::atomic<bool> metadata_conflict_{false};
};

template <typename GraphType>
void write_transient_graph_set(
	const std::filesystem::path& path,
	const PartitionedTransientGraphSet<GraphType>& graphs
) {
	std::array<
		std::uint64_t,
		PartitionedTransientGraphSet<GraphType>::SLOT_COUNT + 1
	> offsets{};
	for (std::size_t shard_index = 0;
		 shard_index < PartitionedTransientGraphSet<GraphType>::SLOT_COUNT;
		 ++shard_index) {
		offsets[shard_index + 1]
			= offsets[shard_index] + graphs.slot(shard_index).size();
	}

	const std::uint64_t graph_count = offsets.back();
	GraphGeneration::MappedSupportTransientWriter<GraphType> writer(path, graph_count);

	#pragma omp parallel for schedule(static)
	for (std::int32_t shard_index = 0;
		 shard_index < static_cast<std::int32_t>(
			PartitionedTransientGraphSet<GraphType>::SLOT_COUNT
		 );
		 ++shard_index) {
		std::uint64_t output_index = offsets[static_cast<std::size_t>(shard_index)];
		const auto& slot = graphs.slot(static_cast<std::size_t>(shard_index));
		for (std::size_t index = 0; index < slot.capacity(); ++index) {
			const GraphType& graph = slot.data()[index];
			if (!graph.empty()) writer.write(output_index++, graph);
		}
	}
}

template <typename GraphType>
void write_final_graph_set(
	const std::filesystem::path& path,
	const ConcurrentFinalGraphSet<GraphType>& graphs
) {
	std::array<
		std::uint64_t,
		ConcurrentFinalGraphSet<GraphType>::SHARD_COUNT + 1
	> offsets{};
	for (std::size_t shard_index = 0;
		 shard_index < ConcurrentFinalGraphSet<GraphType>::SHARD_COUNT;
		 ++shard_index) {
		offsets[shard_index + 1]
			= offsets[shard_index] + graphs.shard(shard_index).size();
	}

	GraphGeneration::MappedFinalGraphWriter<GraphType> writer(path, offsets.back());
	#pragma omp parallel for schedule(static)
	for (std::int32_t shard_index = 0;
		 shard_index < static_cast<std::int32_t>(
			ConcurrentFinalGraphSet<GraphType>::SHARD_COUNT
		 );
		 ++shard_index) {
		std::uint64_t output_index = offsets[static_cast<std::size_t>(shard_index)];
		for (const auto& [graph, metadata]
			: graphs.shard(static_cast<std::size_t>(shard_index))) {
			writer.write(
				output_index++,
				GraphGeneration::final_graph_record<GraphType>{
					.graph = graph,
					.automorphism_order = metadata.automorphism_order,
					.survival = metadata.survival
				}
			);
		}
	}
}

template <int Vertices>
bool is_simple(const AdmissibleGraph<Vertices>& graph) {
	std::array<std::array<bool, Vertices>, Vertices> seen{};
	for (Int edge = 0; edge < AdmissibleGraph<Vertices>::N_EDGES_; ++edge) {
		auto [first, second] = graph.getEdge(edge);
		if (first == second) {
			return false;
		}
		if (second < first) {
			std::swap(first, second);
		}
		if (seen[static_cast<std::size_t>(first)][static_cast<std::size_t>(second)]) {
			return false;
		}
		seen[static_cast<std::size_t>(first)][static_cast<std::size_t>(second)] = true;
	}
	return true;
}

template <std::size_t Vertices>
bool adjacency_is_vertex_irreducible(
	const std::array<std::uint64_t, Vertices>& adjacency
) {
	static_assert(Vertices <= 64);
	if constexpr (Vertices <= 1) {
		return true;
	} else {
	const std::uint64_t all_vertices = Vertices == 64
		? ~std::uint64_t{0}
		: (std::uint64_t{1} << Vertices) - 1;
	for (Int removed = 0; removed < Vertices; ++removed) {
		const std::uint64_t remaining
			= all_vertices & ~(std::uint64_t{1} << removed);
		std::uint64_t reached
			= std::uint64_t{1} << std::countr_zero(remaining);
		std::uint64_t frontier = reached;
		while (frontier != 0) {
			const std::size_t vertex = static_cast<std::size_t>(
				std::countr_zero(frontier)
			);
			frontier &= frontier - 1;
			const std::uint64_t discovered
				= adjacency[vertex] & remaining & ~reached;
			reached |= discovered;
			frontier |= discovered;
		}
		if (reached != remaining) {
			return false;
		}
	}
	return true;
	}
}

template <std::size_t Vertices>
bool adjacency_is_triconnected(
	const std::array<std::uint64_t, Vertices>& adjacency
) {
	static_assert(Vertices <= 64);
	if constexpr (Vertices < 4) {
		return false;
	} else {
	if (!adjacency_is_vertex_irreducible(adjacency)) {
		return false;
	}
	const std::uint64_t all_vertices = Vertices == 64
		? ~std::uint64_t{0}
		: (std::uint64_t{1} << Vertices) - 1;
	for (Int first = 0; first < Vertices; ++first) {
		for (Int second = static_cast<Int>(first + 1);
			 second < Vertices;
			 ++second) {
			const std::uint64_t remaining = all_vertices
				& ~(std::uint64_t{1} << first)
				& ~(std::uint64_t{1} << second);
			std::uint64_t reached
				= std::uint64_t{1} << std::countr_zero(remaining);
			std::uint64_t frontier = reached;
			while (frontier != 0) {
				const std::size_t vertex = static_cast<std::size_t>(
					std::countr_zero(frontier)
				);
				frontier &= frontier - 1;
				const std::uint64_t discovered
					= adjacency[vertex] & remaining & ~reached;
				reached |= discovered;
				frontier |= discovered;
			}
			if (reached != remaining) {
				return false;
			}
		}
	}
	return true;
	}
}

template <int Vertices>
bool is_vertex_irreducible(const AdmissibleGraph<Vertices>& graph) {
	std::array<std::uint64_t, Vertices> adjacency{};
	for (Int edge = 0; edge < AdmissibleGraph<Vertices>::N_EDGES_; ++edge) {
		const auto [first, second] = graph.getEdge(edge);
		if (first == second) {
			continue;
		}
		adjacency[static_cast<std::size_t>(first)]
			|= std::uint64_t{1} << second;
		adjacency[static_cast<std::size_t>(second)]
			|= std::uint64_t{1} << first;
	}
	return adjacency_is_vertex_irreducible(adjacency);
}

template <int Vertices>
bool is_triconnected(const AdmissibleGraph<Vertices>& graph) {
	std::array<std::uint64_t, Vertices> adjacency{};
	for (Int edge = 0; edge < AdmissibleGraph<Vertices>::N_EDGES_; ++edge) {
		const auto [first, second] = graph.getEdge(edge);
		if (first == second) {
			continue;
		}
		adjacency[static_cast<std::size_t>(first)]
			|= std::uint64_t{1} << second;
		adjacency[static_cast<std::size_t>(second)]
			|= std::uint64_t{1} << first;
	}
	return adjacency_is_triconnected(adjacency);
}

template <typename Support>
bool support_is_triconnected(const Support& graph) {
	return adjacency_is_triconnected(graph.support_adjacency());
}

template <typename Support>
bool support_is_vertex_irreducible(const Support& graph) {
	return adjacency_is_vertex_irreducible(graph.support_adjacency());
}

void initialize_dimension_log(const std::filesystem::path& path) {
	std::ofstream output(path, std::ios::trunc);
	if (!output) {
		throw std::runtime_error("could not create dimension log: " + path.string());
	}
	output << "loop\tvertices\tedges\tunsigned_classes"
	          "\todd_edge_classes\teven_edge_odd_vertex_classes\n"
	       << std::flush;
}

void append_dimension_log(
	const std::filesystem::path& path,
	int vertices,
	std::uint64_t unsigned_classes,
	std::uint64_t odd_edge_classes,
	std::uint64_t even_edge_odd_vertex_classes
) {
	std::ofstream output(path, std::ios::app);
	if (!output) {
		throw std::runtime_error("could not open dimension log: " + path.string());
	}
	output << LOOP_NUMBER << '\t'
	       << vertices << '\t'
	       << LOOP_NUMBER + vertices - 1 << '\t'
	       << unsigned_classes << '\t'
	       << odd_edge_classes << '\t'
	       << even_edge_odd_vertex_classes << '\n'
	       << std::flush;
	if (!output) {
		throw std::runtime_error("could not append dimension log: " + path.string());
	}
}

std::filesystem::path transient_path(
	const std::filesystem::path& directory,
	int vertices
) {
	return directory
		/ ("loop_" + std::to_string(LOOP_NUMBER)
			+ "_vertices_" + std::to_string(vertices)
			+ "_transient.gcg");
}

std::filesystem::path admissible_path(
	const std::filesystem::path& directory,
	int vertices
) {
	return directory
		/ ("loop_" + std::to_string(LOOP_NUMBER)
			+ "_vertices_" + std::to_string(vertices)
			+ "_admissible.gcg");
}

void check_header_for_build(const GraphGeneration::MappedGraphFileHeader& header) {
	if (header.vertices == 0 || header.edges + 1 != header.vertices + LOOP_NUMBER) {
		throw std::runtime_error("input graph dimensions do not have this build's loop number");
	}
}

void require_payload_kind(
	const GraphGeneration::MappedGraphFileHeader& header,
	GraphGeneration::MappedGraphPayloadKind expected,
	const std::string& description
) {
	if (header.format_version != GraphGeneration::MAPPED_GRAPH_FORMAT_VERSION
		|| header.payload_kind != static_cast<std::uint64_t>(expected)) {
		throw std::runtime_error("input is not a " + description + " file");
	}
}

struct GeneratedStageResult {
	std::uint64_t transient_count = 0;
	std::uint64_t admissible_count = 0;
	std::uint64_t odd_edge_count = 0;
	std::uint64_t even_edge_odd_vertex_count = 0;
};

template <int Vertices>
GeneratedStageResult generate_next_stage(
	const std::filesystem::path& input_path,
	const std::filesystem::path& output_transient_path,
	const std::filesystem::path& output_admissible_path,
	const std::filesystem::path* dimension_log
) {
	using InputTransient = Transient<Vertices>;
	using OutputTransient = Transient<Vertices + 1>;
	using OutputAdmissible = AdmissibleGraph<Vertices + 1>;

	const auto start = std::chrono::steady_clock::now();
	GraphGeneration::MappedSupportTransientReader<InputTransient> input(input_path);
	if (input.splittable_size() != input.size()) {
		throw std::runtime_error(
			"support-transient input must contain only active records"
		);
	}
	if (input.size() > static_cast<std::uint64_t>(
			std::numeric_limits<std::int64_t>::max())) {
		throw std::overflow_error("too many transient graphs for the OpenMP loop index");
	}

	PartitionedTransientGraphSet<OutputTransient> output_transients;
	ConcurrentFinalGraphSet<OutputAdmissible> output_admissible;
	[[maybe_unused]] GenerationStageProfile stage_profile;
	const std::int64_t input_size = static_cast<std::int64_t>(input.size());
	std::atomic<bool> generation_failed{false};
	std::exception_ptr generation_error;
	const auto generation_start = std::chrono::steady_clock::now();

	#pragma omp parallel
	{
		GenerationStageProfile local_profile;
#if defined(GC_PROFILE_GRAPH_GENERATION)
		local_profile.workers = 1;
#endif
		GraphGeneration::unrooted_support_transient_standardizer<
			OutputTransient::N_VERTICES_, OutputTransient::N_EDGES_
		> transient_standardizer;
		GraphGeneration::final_graph_canonicalizer<OutputAdmissible>
			admissible_standardizer;

		#pragma omp for schedule(dynamic, 8)
		for (std::int64_t i = 0; i < input_size; ++i) {
			if (generation_failed.load(std::memory_order_relaxed)) {
				continue;
			}
			try {
				const InputTransient parent = input[static_cast<std::uint64_t>(i)];
				GraphGeneration::unrooted_support_transient_graph<
					InputTransient::N_VERTICES_,
					InputTransient::N_EDGES_
				> transient(parent);
				transient.for_each_allowed_root_split(
					[&](auto split, Int /* split_root */) {
#if defined(GC_PROFILE_GRAPH_GENERATION)
						++local_profile.relevant_children;
#endif
						auto& child_state = split.child;
						if (child_state.support().has_simple_expansion()) {
							const OutputAdmissible simple_graph
								= child_state.support().expand_simple().to_hairless_graph();
							if (is_vertex_irreducible<OutputTransient::N_VERTICES_>(
									simple_graph
								)) {
								auto canonical = profiled_call(
									local_profile.admissible_canonicalization,
									[&]() {
										return admissible_standardizer(simple_graph);
									}
								);
								profiled_call(
									local_profile.admissible_map_insertion,
									[&]() {
										output_admissible.insert(std::move(canonical));
									}
								);
							}
						}

						if (child_state.has_active_root_choice()) {
							auto canonical = profiled_call(
								local_profile.transient_canonicalization,
								[&]() {
									return transient_standardizer.standardize_with_info(
										child_state
									);
								}
							);
							profiled_call(
								local_profile.transient_map_insertion,
								[&]() {
									output_transients.insert(
										std::move(canonical.canonical_graph.support()),
										canonical.maximum_valence,
										canonical.maximum_valence_count
									);
								}
							);
						}
					}
				);
			} catch (...) {
				generation_failed.store(true, std::memory_order_relaxed);
				#pragma omp critical(gc_graph_generation_error)
				{
					if (generation_error == nullptr) {
						generation_error = std::current_exception();
					}
				}
			}
		}
#if defined(GC_PROFILE_GRAPH_GENERATION)
		#pragma omp critical(gc_graph_generation_profile)
		stage_profile.merge(local_profile);
#endif
	}
	const double generation_wall_seconds = std::chrono::duration<double>(
		std::chrono::steady_clock::now() - generation_start
	).count();
	if (generation_error != nullptr) {
		std::rethrow_exception(generation_error);
	}
	if (output_admissible.has_metadata_conflict()) {
		throw std::logic_error(
			"duplicate final graphs produced inconsistent canonical metadata"
		);
	}

	const auto transient_map_write_start = std::chrono::steady_clock::now();
	write_transient_graph_set(output_transient_path, output_transients);
	const double transient_map_write_wall_seconds = std::chrono::duration<double>(
		std::chrono::steady_clock::now() - transient_map_write_start
	).count();
	const double transient_seconds = std::chrono::duration<double>(
		std::chrono::steady_clock::now() - start
	).count();
	std::cout << "completed file=" << output_transient_path.string()
	          << " kind=" << TRANSIENT_LOG_KIND
	          << " loop=" << LOOP_NUMBER
	          << " vertices=" << Vertices + 1
	          << " edges=" << +OutputTransient::N_EDGES_
	          << " graphs=" << output_transients.size()
	          << " duplicates=" << output_transients.duplicates()
	          << " stage_elapsed_seconds=" << transient_seconds << '\n'
	          << std::flush;
	for (Int maximum = 0; maximum <= OutputTransient::N_VERTICES_; ++maximum) {
		for (Int tied = 1; tied <= OutputTransient::N_VERTICES_; ++tied) {
			const std::uint64_t count
				= output_transients.partition_size(maximum, tied);
			if (count != 0) {
				std::cout << "transient_partition"
				          << " vertices=" << Vertices + 1
				          << " edges=" << +OutputTransient::N_EDGES_
				          << " maximum_support_valence=" << +maximum
				          << " tied_maximum_vertices=" << +tied
				          << " graphs=" << count << '\n';
			}
		}
	}
	std::cout << std::flush;

	const auto admissible_map_write_start = std::chrono::steady_clock::now();
	write_final_graph_set(
		output_admissible_path,
		output_admissible
	);
	const double admissible_map_write_wall_seconds = std::chrono::duration<double>(
		std::chrono::steady_clock::now() - admissible_map_write_start
	).count();
	const auto admissible_map_scan_start = std::chrono::steady_clock::now();
	const std::uint64_t odd_edge_count = output_admissible.surviving(
		GraphGeneration::final_graph_survival::odd_edges
	);
	const std::uint64_t even_edge_odd_vertex_count = output_admissible.surviving(
		GraphGeneration::final_graph_survival::even_edges_odd_vertices
	);
	const double admissible_map_scan_wall_seconds = std::chrono::duration<double>(
		std::chrono::steady_clock::now() - admissible_map_scan_start
	).count();
	const double total_seconds = std::chrono::duration<double>(
		std::chrono::steady_clock::now() - start
	).count();
	std::cout << "completed file=" << output_admissible_path.string()
	          << " kind=admissible"
	          << " loop=" << LOOP_NUMBER
	          << " vertices=" << Vertices + 1
	          << " edges=" << +OutputAdmissible::N_EDGES_
	          << " graphs=" << output_admissible.size()
	          << " odd_edge_graphs=" << odd_edge_count
	          << " even_edge_odd_vertex_graphs=" << even_edge_odd_vertex_count
	          << " duplicates=" << output_admissible.duplicates()
	          << " stage_elapsed_seconds=" << total_seconds << '\n'
	          << std::flush;

#if defined(GC_PROFILE_GRAPH_GENERATION)
	print_generation_profile(
		Vertices + 1,
		OutputTransient::N_EDGES_,
		stage_profile,
		generation_wall_seconds,
		transient_map_write_wall_seconds,
		admissible_map_write_wall_seconds,
		admissible_map_scan_wall_seconds
	);
#else
	(void)generation_wall_seconds;
	(void)transient_map_write_wall_seconds;
	(void)admissible_map_write_wall_seconds;
	(void)admissible_map_scan_wall_seconds;
#endif

	if (dimension_log != nullptr) {
		append_dimension_log(
			*dimension_log,
			Vertices + 1,
			output_admissible.size(),
			odd_edge_count,
			even_edge_odd_vertex_count
		);
	}
	return {
		output_transients.size(),
		output_admissible.size(),
		odd_edge_count,
		even_edge_odd_vertex_count
	};
}

template <int CandidateVertices = SEED_VERTICES>
GeneratedStageResult dispatch_next_stage(
	std::uint64_t vertices,
	const std::filesystem::path& input_path,
	const std::filesystem::path& output_transient_path,
	const std::filesystem::path& output_admissible_path,
	const std::filesystem::path* dimension_log = nullptr
) {
	if (vertices == CandidateVertices) {
		return generate_next_stage<CandidateVertices>(
			input_path,
			output_transient_path,
			output_admissible_path,
			dimension_log
		);
	}
	if constexpr (CandidateVertices < MAX_VERTICES) {
		return dispatch_next_stage<CandidateVertices + 1>(
			vertices,
			input_path,
			output_transient_path,
			output_admissible_path,
			dimension_log
		);
	}
	throw std::runtime_error("unsupported vertex count for this loop-number build");
}

struct StageStats {
	std::uint64_t transient_classes = 0;
	std::uint64_t admissible_classes = 0;
	std::uint64_t odd_edge_classes = 0;
	std::uint64_t even_edge_odd_vertex_classes = 0;
};

template <int Vertices>
StageStats measure_stage(
	const std::filesystem::path& transient_file,
	const std::filesystem::path& admissible_file
) {
	using Support = Transient<Vertices>;
	using Hairless = AdmissibleGraph<Vertices>;
	GraphGeneration::MappedSupportTransientReader<Support> transients(
		transient_file
	);
	GraphGeneration::MappedFinalGraphReader<Hairless> admissible(admissible_file);

	StageStats stats;
	stats.transient_classes = transients.size();
	stats.admissible_classes = admissible.size();
	for (std::uint64_t i = 0; i < admissible.size(); ++i) {
		const auto record = admissible[i];
		stats.odd_edge_classes += GraphGeneration::survives(
			record.survival,
			GraphGeneration::final_graph_survival::odd_edges
		);
		stats.even_edge_odd_vertex_classes += GraphGeneration::survives(
			record.survival,
			GraphGeneration::final_graph_survival::even_edges_odd_vertices
		);
	}
	return stats;
}

template <int CandidateVertices = SEED_VERTICES>
StageStats dispatch_measure_stage(
	std::uint64_t vertices,
	const std::filesystem::path& transient_file,
	const std::filesystem::path& admissible_file
) {
	if (vertices == CandidateVertices) {
		return measure_stage<CandidateVertices>(transient_file, admissible_file);
	}
	if constexpr (CandidateVertices < MAX_VERTICES) {
		return dispatch_measure_stage<CandidateVertices + 1>(
			vertices, transient_file, admissible_file
		);
	}
	throw std::runtime_error("unsupported vertex count for this loop-number build");
}

void print_stage_stats(std::uint64_t vertices, const StageStats& stats) {
	std::cout << "loop number = " << LOOP_NUMBER << '\n';
	std::cout << "vertices = " << vertices << '\n';
	std::cout << "transient classes = " << stats.transient_classes << '\n';
	std::cout << "admissible classes = " << stats.admissible_classes << '\n';
	std::cout << "odd-edge classes = " << stats.odd_edge_classes << '\n';
	std::cout << "even-edge/odd-vertex classes = "
	          << stats.even_edge_odd_vertex_classes << '\n';
}

template <int Vertices>
bool same_simple_classes(
	const std::filesystem::path& final_path,
	const std::filesystem::path& reference_path
) {
	using GraphType = AdmissibleGraph<Vertices>;
	GraphGeneration::MappedFinalGraphReader<GraphType> final_graphs(final_path);
	std::unordered_set<GraphType> final_classes;
	for (std::uint64_t i = 0; i < final_graphs.size(); ++i) {
		final_classes.insert(final_graphs[i].graph);
	}

	GraphGeneration::MappedGraphReader<GraphType> reference(reference_path);
	Standardizer<GraphType> standardizer;
	std::unordered_set<GraphType> reference_classes;
	for (std::uint64_t i = 0; i < reference.size(); ++i) {
		const GraphType graph = reference[i];
		if (is_simple<Vertices>(graph)) {
			reference_classes.insert(standardizer.standardize_no_sign(graph));
		}
	}
	if (final_classes == reference_classes) {
		return true;
	}
	std::cout << "generated=" << final_classes.size()
	          << " expected=" << reference_classes.size() << '\n';
	for (const GraphType& graph : reference_classes) {
		if (!final_classes.contains(graph)) {
			std::cout << "missing=";
			for (Int edge = 0; edge < GraphType::N_EDGES_; ++edge) {
				const auto [first, second] = graph.getEdge(edge);
				if (edge != 0) std::cout << ',';
				std::cout << +first << '-' << +second;
			}
			std::cout << '\n';
		}
	}
	for (const GraphType& graph : final_classes) {
		if (!reference_classes.contains(graph)) {
			std::cout << "extra=";
			for (Int edge = 0; edge < GraphType::N_EDGES_; ++edge) {
				const auto [first, second] = graph.getEdge(edge);
				if (edge != 0) std::cout << ',';
				std::cout << +first << '-' << +second;
			}
			std::cout << '\n';
		}
	}
	return false;
}

template <int Vertices>
bool same_vertex_irreducible_classes(
	const std::filesystem::path& final_path,
	const std::filesystem::path& reference_path
) {
	using GraphType = AdmissibleGraph<Vertices>;
	GraphGeneration::MappedFinalGraphReader<GraphType> final_graphs(final_path);
	std::unordered_set<GraphType> final_classes;
	for (std::uint64_t i = 0; i < final_graphs.size(); ++i) {
		const GraphType graph = final_graphs[i].graph;
		if (!is_vertex_irreducible<Vertices>(graph)) {
			return false;
		}
		final_classes.insert(graph);
	}

	Standardizer<GraphType> standardizer;
	std::unordered_set<GraphType> reference_classes;
	const auto insert_reference = [&](const GraphType& graph) {
		if (is_simple<Vertices>(graph)
			&& is_vertex_irreducible<Vertices>(graph)) {
			reference_classes.insert(standardizer.standardize_no_sign(graph));
		}
	};
	const auto header
		= GraphGeneration::read_mapped_graph_file_header(reference_path);
	if (header.format_version == GraphGeneration::MAPPED_GRAPH_FORMAT_VERSION
		&& header.payload_kind == static_cast<std::uint64_t>(
			GraphGeneration::MappedGraphPayloadKind::final_graph_with_metadata
		)) {
		GraphGeneration::MappedFinalGraphReader<GraphType> reference(reference_path);
		for (std::uint64_t i = 0; i < reference.size(); ++i) {
			insert_reference(reference[i].graph);
		}
	} else {
		GraphGeneration::MappedGraphReader<GraphType> reference(reference_path);
		for (std::uint64_t i = 0; i < reference.size(); ++i) {
			insert_reference(reference[i]);
		}
	}
	if (final_classes == reference_classes) {
		return true;
	}
	std::cout << "generated=" << final_classes.size()
	          << " expected=" << reference_classes.size() << '\n';
	for (const GraphType& graph : reference_classes) {
		if (!final_classes.contains(graph)) {
			std::cout << "missing=";
			for (Int edge = 0; edge < GraphType::N_EDGES_; ++edge) {
				const auto [first, second] = graph.getEdge(edge);
				if (edge != 0) std::cout << ',';
				std::cout << +first << '-' << +second;
			}
			std::cout << '\n';
		}
	}
	for (const GraphType& graph : final_classes) {
		if (!reference_classes.contains(graph)) {
			std::cout << "extra=";
			for (Int edge = 0; edge < GraphType::N_EDGES_; ++edge) {
				const auto [first, second] = graph.getEdge(edge);
				if (edge != 0) std::cout << ',';
				std::cout << +first << '-' << +second;
			}
			std::cout << '\n';
		}
	}
	return false;
}

template <int Vertices>
bool same_triconnected_classes(
	const std::filesystem::path& final_path,
	const std::filesystem::path& reference_path
) {
	using GraphType = AdmissibleGraph<Vertices>;
	GraphGeneration::MappedFinalGraphReader<GraphType> final_graphs(final_path);
	std::unordered_set<GraphType> final_classes;
	for (std::uint64_t i = 0; i < final_graphs.size(); ++i) {
		const GraphType graph = final_graphs[i].graph;
		if (!is_simple<Vertices>(graph)
			|| !is_triconnected<Vertices>(graph)) {
			return false;
		}
		final_classes.insert(graph);
	}

	Standardizer<GraphType> standardizer;
	std::unordered_set<GraphType> reference_classes;
	const auto insert_reference = [&](const GraphType& graph) {
		if (is_simple<Vertices>(graph) && is_triconnected<Vertices>(graph)) {
			reference_classes.insert(standardizer.standardize_no_sign(graph));
		}
	};
	const auto header
		= GraphGeneration::read_mapped_graph_file_header(reference_path);
	if (header.format_version == GraphGeneration::MAPPED_GRAPH_FORMAT_VERSION
		&& header.payload_kind == static_cast<std::uint64_t>(
			GraphGeneration::MappedGraphPayloadKind::final_graph_with_metadata
		)) {
		GraphGeneration::MappedFinalGraphReader<GraphType> reference(reference_path);
		for (std::uint64_t i = 0; i < reference.size(); ++i) {
			insert_reference(reference[i].graph);
		}
	} else {
		GraphGeneration::MappedGraphReader<GraphType> reference(reference_path);
		for (std::uint64_t i = 0; i < reference.size(); ++i) {
			insert_reference(reference[i]);
		}
	}
	return final_classes == reference_classes;
}

template <int CandidateVertices = SEED_VERTICES>
bool dispatch_same_simple_classes(
	std::uint64_t vertices,
	const std::filesystem::path& left_path,
	const std::filesystem::path& right_path
) {
	if (vertices == CandidateVertices) {
		return same_simple_classes<CandidateVertices>(left_path, right_path);
	}
	if constexpr (CandidateVertices < MAX_VERTICES) {
		return dispatch_same_simple_classes<CandidateVertices + 1>(
			vertices, left_path, right_path
		);
	}
	throw std::runtime_error("unsupported vertex count for this loop-number build");
}

template <int CandidateVertices = SEED_VERTICES>
bool dispatch_same_vertex_irreducible_classes(
	std::uint64_t vertices,
	const std::filesystem::path& left_path,
	const std::filesystem::path& right_path
) {
	if (vertices == CandidateVertices) {
		return same_vertex_irreducible_classes<CandidateVertices>(
			left_path, right_path
		);
	}
	if constexpr (CandidateVertices < MAX_VERTICES) {
		return dispatch_same_vertex_irreducible_classes<CandidateVertices + 1>(
			vertices, left_path, right_path
		);
	}
	throw std::runtime_error("unsupported vertex count for this loop-number build");
}

template <int CandidateVertices = 4>
bool dispatch_same_triconnected_classes(
	std::uint64_t vertices,
	const std::filesystem::path& left_path,
	const std::filesystem::path& right_path
) {
	if (vertices == CandidateVertices) {
		return same_triconnected_classes<CandidateVertices>(
			left_path, right_path
		);
	}
	if constexpr (CandidateVertices < MAX_VERTICES) {
		return dispatch_same_triconnected_classes<CandidateVertices + 1>(
			vertices, left_path, right_path
		);
	}
	throw std::runtime_error("unsupported vertex count for this loop-number build");
}

template <int Vertices>
bool same_transient_classes(
	const std::filesystem::path& left_path,
	const std::filesystem::path& right_path
) {
	using GraphType = Transient<Vertices>;
	auto load = [](const std::filesystem::path& path) {
		GraphGeneration::MappedSupportTransientReader<GraphType> input(path);
		std::unordered_set<GraphType> classes;
		for (std::uint64_t i = 0; i < input.size(); ++i) {
			classes.insert(input[i]);
		}
		return classes;
	};
	return load(left_path) == load(right_path);
}

template <int CandidateVertices = SEED_VERTICES>
bool dispatch_same_transient_classes(
	std::uint64_t vertices,
	const std::filesystem::path& left_path,
	const std::filesystem::path& right_path
) {
	if (vertices == CandidateVertices) {
		return same_transient_classes<CandidateVertices>(left_path, right_path);
	}
	if constexpr (CandidateVertices < MAX_VERTICES) {
		return dispatch_same_transient_classes<CandidateVertices + 1>(
			vertices, left_path, right_path
		);
	}
	throw std::runtime_error("unsupported vertex count for this loop-number build");
}

template <int Vertices>
void validate_stage_files(
	const std::filesystem::path& transient_file,
	const std::filesystem::path& final_file
) {
	using Support = Transient<Vertices>;
	using Hairless = AdmissibleGraph<Vertices>;
	GraphGeneration::MappedSupportTransientReader<Support> transients(
		transient_file
	);
	GraphGeneration::MappedFinalGraphReader<Hairless> final_graphs(final_file);
	if (transients.splittable_size() != transients.size()) {
		throw std::runtime_error("transient file header violates active-frontier invariants");
	}

	GraphGeneration::unrooted_support_transient_standardizer<
		Support::N_VERTICES_, Support::N_EDGES_
	> transient_standardizer;
	for (std::uint64_t i = 0; i < transients.size(); ++i) {
		const Support graph = transients[i];
		GraphGeneration::unrooted_support_transient_graph<
			Support::N_VERTICES_, Support::N_EDGES_
		> state(graph);
		if (transient_standardizer.standardize_no_sign(state).support() != graph) {
			throw std::runtime_error(
				"transient record is not canonical at index " + std::to_string(i)
			);
		}
		if (!state.has_active_root_choice()) {
			throw std::runtime_error(
				"transient record has no splittable maximum-valence vertex at index "
					+ std::to_string(i)
			);
		}
	}

	GraphGeneration::final_graph_canonicalizer<Hairless> finalizer;
	for (std::uint64_t i = 0; i < final_graphs.size(); ++i) {
		const auto stored = final_graphs[i];
		if (!is_simple<Vertices>(stored.graph)
			|| !is_vertex_irreducible<Vertices>(stored.graph)) {
			throw std::runtime_error(
				"final graph is not simple and vertex-irreducible at index "
					+ std::to_string(i)
			);
		}
		const auto computed = finalizer(stored.graph);
		if (computed.canonical_graph != stored.graph
			|| computed.automorphism_order != stored.automorphism_order
			|| computed.survival != stored.survival) {
			throw std::runtime_error(
				"final graph metadata is inconsistent at index " + std::to_string(i)
			);
		}
	}
}

template <int CandidateVertices = SEED_VERTICES>
void dispatch_validate_stage_files(
	std::uint64_t vertices,
	const std::filesystem::path& transient_file,
	const std::filesystem::path& final_file
) {
	if (vertices == CandidateVertices) {
		validate_stage_files<CandidateVertices>(transient_file, final_file);
		return;
	}
	if constexpr (CandidateVertices < MAX_VERTICES) {
		dispatch_validate_stage_files<CandidateVertices + 1>(
			vertices, transient_file, final_file
		);
		return;
	}
	throw std::runtime_error("unsupported vertex count for this loop-number build");
}

void write_seed(
	const std::filesystem::path& transient_file,
	const std::filesystem::path& admissible_file,
	const std::filesystem::path* dimension_log = nullptr
) {
	using SeedTransient = Transient<SEED_VERTICES>;
	using SeedAdmissible = AdmissibleGraph<SEED_VERTICES>;
	const auto start = std::chrono::steady_clock::now();
	SeedTransient edge;
	edge.set_support_edge(0, 1);
	edge.validate();
	constexpr std::uint64_t transient_count = 1;
	{
		GraphGeneration::MappedSupportTransientWriter<SeedTransient> writer(
			transient_file, transient_count
		);
		writer.write(0, edge);
	}
	const double transient_seconds = std::chrono::duration<double>(
		std::chrono::steady_clock::now() - start
	).count();
	std::cout << "completed file=" << transient_file.string()
	          << " kind=" << TRANSIENT_LOG_KIND
	          << " loop=" << LOOP_NUMBER
	          << " vertices=" << SEED_VERTICES
	          << " edges=" << +SeedTransient::N_EDGES_
	          << " graphs=" << transient_count
	          << " duplicates=0"
	          << " stage_elapsed_seconds=" << transient_seconds << '\n'
	          << std::flush;

	{
		GraphGeneration::MappedFinalGraphWriter<SeedAdmissible> writer(
			admissible_file, 0
		);
	}
	const double total_seconds = std::chrono::duration<double>(
		std::chrono::steady_clock::now() - start
	).count();
	std::cout << "completed file=" << admissible_file.string()
	          << " kind=admissible"
	          << " loop=" << LOOP_NUMBER
	          << " vertices=" << SEED_VERTICES
	          << " edges=" << +SeedAdmissible::N_EDGES_
	          << " graphs=0"
	          << " odd_edge_graphs=0"
	          << " even_edge_odd_vertex_graphs="
	          << 0
	          << " duplicates=0"
	          << " stage_elapsed_seconds=" << total_seconds << '\n'
	          << std::flush;

	if (dimension_log != nullptr) {
		append_dimension_log(
			*dimension_log,
			SEED_VERTICES,
			0,
			0,
			0
		);
	}
}

struct FinalDimensionCounts {
	std::uint64_t unsigned_classes = 0;
	std::uint64_t odd_edge_classes = 0;
	std::uint64_t even_edge_odd_vertex_classes = 0;
};

template <int Vertices>
FinalDimensionCounts measure_final_dimensions(
	const std::filesystem::path& path
) {
	using GraphType = AdmissibleGraph<Vertices>;
	GraphGeneration::MappedFinalGraphReader<GraphType> input(path);
	FinalDimensionCounts counts;
	counts.unsigned_classes = input.size();
	for (std::uint64_t i = 0; i < input.size(); ++i) {
		const auto record = input[i];
		counts.odd_edge_classes += GraphGeneration::survives(
			record.survival,
			GraphGeneration::final_graph_survival::odd_edges
		);
		counts.even_edge_odd_vertex_classes += GraphGeneration::survives(
			record.survival,
			GraphGeneration::final_graph_survival::even_edges_odd_vertices
		);
	}
	return counts;
}

template <int CandidateVertices = SEED_VERTICES>
FinalDimensionCounts dispatch_final_dimensions(
	std::uint64_t vertices,
	const std::filesystem::path& path
) {
	if (vertices == CandidateVertices) {
		return measure_final_dimensions<CandidateVertices>(path);
	}
	if constexpr (CandidateVertices < MAX_VERTICES) {
		return dispatch_final_dimensions<CandidateVertices + 1>(vertices, path);
	}
	throw std::runtime_error("unsupported vertex count for this loop-number build");
}

int print_dimension_rows(const std::filesystem::path& directory) {
	std::cout << "loop\tvertices\tedges\tunsigned_classes"
	          << "\todd_edge_classes\teven_edge_odd_vertex_classes\n";
	for (int vertices = SEED_VERTICES; vertices <= MAX_VERTICES; ++vertices) {
		const std::filesystem::path path = admissible_path(directory, vertices);
		if (!std::filesystem::exists(path)) {
			if (vertices == SEED_VERTICES) {
				throw std::runtime_error("dimension directory contains no final stages");
			}
			break;
		}
		const FinalDimensionCounts counts = dispatch_final_dimensions(
			vertices, path
		);
		std::cout << LOOP_NUMBER << '\t'
		          << vertices << '\t'
		          << LOOP_NUMBER + vertices - 1 << '\t'
		          << counts.unsigned_classes << '\t'
		          << counts.odd_edge_classes << '\t'
		          << counts.even_edge_odd_vertex_classes << '\n';
	}
	return EXIT_SUCCESS;
}

int run_generate(int argc, char** argv) {
	const std::filesystem::path output_directory = argv[2];
	const int max_vertices = argc == 4 ? std::stoi(argv[3]) : MAX_VERTICES;
	if (max_vertices < SEED_VERTICES || max_vertices > MAX_VERTICES) {
		throw std::runtime_error(
			"max_vertices must be between " + std::to_string(SEED_VERTICES)
				+ " and " + std::to_string(MAX_VERTICES)
		);
	}

	std::filesystem::create_directories(output_directory);
	const std::filesystem::path dimension_log
		= output_directory / DIMENSION_LOG_FILENAME;
	initialize_dimension_log(dimension_log);
	std::filesystem::path current_transient = transient_path(
		output_directory, SEED_VERTICES
	);
	write_seed(
		current_transient,
		admissible_path(output_directory, SEED_VERTICES),
		&dimension_log
	);
	for (int vertices = SEED_VERTICES; vertices < max_vertices; ++vertices) {
		const std::filesystem::path next_transient
			= transient_path(output_directory, vertices + 1);
		const GeneratedStageResult result = dispatch_next_stage(
			vertices,
			current_transient,
			next_transient,
			admissible_path(output_directory, vertices + 1),
			&dimension_log
		);
		current_transient = next_transient;
		if (result.transient_count == 0 && vertices + 1 < max_vertices) {
			throw std::runtime_error(
				"transient frontier became empty before the requested terminal stage"
			);
		}
	}
	return EXIT_SUCCESS;
}

void print_usage(const char* argv0) {
	std::cerr << "graph-stage generator compiled for loop number " << LOOP_NUMBER << "\n";
	std::cerr << "usage:\n";
	std::cerr << "  " << argv0
	          << " seed <transient_output> <admissible_output>\n";
	std::cerr << "  " << argv0
	          << " next <transient_input> <transient_output> <admissible_output>\n";
	std::cerr << "  " << argv0
	          << " stats <transient_file> <admissible_file>\n";
	std::cerr << "  " << argv0
	          << " generate <output_directory> [max_vertices]\n";
	std::cerr << "  " << argv0 << " dimensions <output_directory>\n";
	std::cerr << "  " << argv0
	          << " compare-simple <admissible_file> <reference_file>\n";
	std::cerr << "  " << argv0
	          << " compare-irreducible <admissible_file> <reference_file>\n";
	std::cerr << "  " << argv0
	          << " compare-triconnected <admissible_file> <reference_file>\n";
	std::cerr << "  " << argv0
	          << " compare-transient <left_file> <right_file>\n";
	std::cerr << "  " << argv0
	          << " validate <transient_file> <admissible_file>\n";
}

} // namespace

int main(int argc, char** argv) {
	try {
		if (argc == 4 && std::string(argv[1]) == "seed") {
			write_seed(argv[2], argv[3]);
			return EXIT_SUCCESS;
		}
		if (argc == 5 && std::string(argv[1]) == "next") {
			const std::filesystem::path input_path = argv[2];
			const auto header = GraphGeneration::read_mapped_graph_file_header(input_path);
			check_header_for_build(header);
			require_payload_kind(
				header,
				GraphGeneration::MAPPED_SUPPORT_TRANSIENT_PAYLOAD_KIND,
				"unrooted maximum-valence transient graph"
			);
			if (header.vertices >= MAX_VERTICES) {
				throw std::runtime_error("input is already the terminal trivalent stage");
			}
			dispatch_next_stage(header.vertices, input_path, argv[3], argv[4]);
			return EXIT_SUCCESS;
		}
		if (argc == 4 && std::string(argv[1]) == "stats") {
			const auto transient_header
				= GraphGeneration::read_mapped_graph_file_header(argv[2]);
			const auto admissible_header
				= GraphGeneration::read_mapped_graph_file_header(argv[3]);
			check_header_for_build(transient_header);
			check_header_for_build(admissible_header);
			require_payload_kind(
				transient_header,
				GraphGeneration::MAPPED_SUPPORT_TRANSIENT_PAYLOAD_KIND,
				"unrooted maximum-valence transient graph"
			);
			require_payload_kind(
				admissible_header,
				GraphGeneration::MappedGraphPayloadKind::final_graph_with_metadata,
				"final graph metadata"
			);
			if (transient_header.vertices != admissible_header.vertices
				|| transient_header.edges != admissible_header.edges) {
				throw std::runtime_error("stage dimensions differ");
			}
			print_stage_stats(
				transient_header.vertices,
				dispatch_measure_stage(
					transient_header.vertices, argv[2], argv[3]
				)
			);
			return EXIT_SUCCESS;
		}
		if ((argc == 3 || argc == 4) && std::string(argv[1]) == "generate") {
			return run_generate(argc, argv);
		}
		if (argc == 3 && std::string(argv[1]) == "dimensions") {
			return print_dimension_rows(argv[2]);
		}
		if (argc == 4 && std::string(argv[1]) == "compare-simple") {
			const auto left_header
				= GraphGeneration::read_mapped_graph_file_header(argv[2]);
			const auto right_header
				= GraphGeneration::read_mapped_graph_file_header(argv[3]);
			check_header_for_build(left_header);
			check_header_for_build(right_header);
			require_payload_kind(
				left_header,
				GraphGeneration::MappedGraphPayloadKind::final_graph_with_metadata,
				"final graph metadata"
			);
			if (left_header.vertices != right_header.vertices
				|| left_header.edges != right_header.edges) {
				throw std::runtime_error("stage dimensions differ");
			}
			return dispatch_same_simple_classes(
				left_header.vertices, argv[2], argv[3]
			) ? EXIT_SUCCESS : EXIT_FAILURE;
		}
		if (argc == 4 && std::string(argv[1]) == "compare-irreducible") {
			const auto left_header
				= GraphGeneration::read_mapped_graph_file_header(argv[2]);
			const auto right_header
				= GraphGeneration::read_mapped_graph_file_header(argv[3]);
			check_header_for_build(left_header);
			check_header_for_build(right_header);
			require_payload_kind(
				left_header,
				GraphGeneration::MappedGraphPayloadKind::final_graph_with_metadata,
				"final graph metadata"
			);
			if (left_header.vertices != right_header.vertices
				|| left_header.edges != right_header.edges) {
				throw std::runtime_error("stage dimensions differ");
			}
			return dispatch_same_vertex_irreducible_classes(
				left_header.vertices, argv[2], argv[3]
			) ? EXIT_SUCCESS : EXIT_FAILURE;
		}
		if (argc == 4 && std::string(argv[1]) == "compare-triconnected") {
			const auto left_header
				= GraphGeneration::read_mapped_graph_file_header(argv[2]);
			const auto right_header
				= GraphGeneration::read_mapped_graph_file_header(argv[3]);
			check_header_for_build(left_header);
			check_header_for_build(right_header);
			require_payload_kind(
				left_header,
				GraphGeneration::MappedGraphPayloadKind::final_graph_with_metadata,
				"final graph metadata"
			);
			if (left_header.vertices != right_header.vertices
				|| left_header.edges != right_header.edges) {
				throw std::runtime_error("stage dimensions differ");
			}
			return dispatch_same_triconnected_classes(
				left_header.vertices, argv[2], argv[3]
			) ? EXIT_SUCCESS : EXIT_FAILURE;
		}
		if (argc == 4 && std::string(argv[1]) == "compare-transient") {
			const auto left_header
				= GraphGeneration::read_mapped_graph_file_header(argv[2]);
			const auto right_header
				= GraphGeneration::read_mapped_graph_file_header(argv[3]);
			check_header_for_build(left_header);
			check_header_for_build(right_header);
			require_payload_kind(
				left_header,
				GraphGeneration::MAPPED_SUPPORT_TRANSIENT_PAYLOAD_KIND,
				"unrooted maximum-valence transient graph"
			);
			require_payload_kind(
				right_header,
				GraphGeneration::MAPPED_SUPPORT_TRANSIENT_PAYLOAD_KIND,
				"unrooted maximum-valence transient graph"
			);
			if (left_header.vertices != right_header.vertices
				|| left_header.edges != right_header.edges) {
				throw std::runtime_error("stage dimensions differ");
			}
			return dispatch_same_transient_classes(
				left_header.vertices, argv[2], argv[3]
			) ? EXIT_SUCCESS : EXIT_FAILURE;
		}
		if (argc == 4 && std::string(argv[1]) == "validate") {
			const auto transient_header
				= GraphGeneration::read_mapped_graph_file_header(argv[2]);
			const auto final_header
				= GraphGeneration::read_mapped_graph_file_header(argv[3]);
			check_header_for_build(transient_header);
			check_header_for_build(final_header);
			require_payload_kind(
				transient_header,
				GraphGeneration::MAPPED_SUPPORT_TRANSIENT_PAYLOAD_KIND,
				"unrooted maximum-valence transient graph"
			);
			require_payload_kind(
				final_header,
				GraphGeneration::MappedGraphPayloadKind::final_graph_with_metadata,
				"final graph metadata"
			);
			if (transient_header.vertices != final_header.vertices
				|| transient_header.edges != final_header.edges) {
				throw std::runtime_error("stage dimensions differ");
			}
			dispatch_validate_stage_files(
				transient_header.vertices, argv[2], argv[3]
			);
			std::cout << "validated loop=" << LOOP_NUMBER
			          << " vertices=" << transient_header.vertices
			          << " edges=" << transient_header.edges
			          << " transient_graphs=" << transient_header.graph_count
			          << " final_graphs=" << final_header.graph_count << '\n';
			return EXIT_SUCCESS;
		}

		print_usage(argv[0]);
		return EXIT_FAILURE;
	} catch (const std::exception& error) {
		std::cerr << "error: " << error.what() << '\n';
		return EXIT_FAILURE;
	}
}
