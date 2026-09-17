#include "GraphGeneration/TriangleSplitCommon.hpp"
#include "GraphGeneration/ValenceSplitPlan.hpp"

namespace {
using namespace GraphGeneration;
using namespace GraphGeneration::triangle_split_detail;
using Groups = std::map<ValenceArray, std::size_t>;

std::string group_name(const ValenceArray& degrees) {
    std::string name;
    for (auto degree : degrees) {
        if (!name.empty()) name += '-';
        name += std::to_string(degree);
    }
    return name + ".bin";
}

template <int V>
void write_record(std::ostream& out, const SplitStageGraph<V>& graph) {
    out.write(reinterpret_cast<const char*>(graph.half_edges.data()), graph.half_edges.size());
}

// Read bounded batches, even if one parent group is larger than available RAM.
// These private scratch files contain only canonical half-edge bytes.
template <int V, class Function>
void read_group(const std::filesystem::path& path, std::size_t count, Function function) {
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open parent group: " + path.string());
    constexpr std::size_t batch_size = 4096;
    std::vector<SplitStageGraph<V>> batch(batch_size);
    while (count) {
        const auto size = std::min(count, batch_size);
        for (std::size_t i = 0; i < size; ++i) {
            auto& edges = batch[i].half_edges;
            if (!input.read(reinterpret_cast<char*>(edges.data()), edges.size()))
                throw std::runtime_error("truncated parent group: " + path.string());
        }
        std::exception_ptr error;
        std::atomic<bool> failed{false};
        #pragma omp parallel for schedule(dynamic, 8)
        for (std::size_t i = 0; i < size; ++i) {
            if (failed.load(std::memory_order_relaxed)) continue;
            try { function(batch[i]); }
            catch (...) {
                failed.store(true, std::memory_order_relaxed);
                #pragma omp critical(valence_split_error)
                { if (!error) error = std::current_exception(); }
            }
        }
        if (error) std::rethrow_exception(error);
        count -= size;
    }
    if (input.peek() != std::char_traits<char>::eof() || input.bad())
        throw std::runtime_error("parent group count mismatch: " + path.string());
}

class Pipeline {
    std::filesystem::path input_, output_, scratch_;
    std::map<int, std::size_t> seed_counts_;
    std::ofstream counts_, group_counts_;

    std::filesystem::path seed_path(int v) const {
        return triangle_seed_path(input_, v);
    }
    std::filesystem::path stage_path(int v) const {
        return scratch_ / ("V" + std::to_string(v));
    }

    template <int V>
    Groups partition_seeds(const std::filesystem::path& directory) {
        Groups groups;
        if (!seed_counts_.contains(V)) return groups;
        std::ifstream input(seed_path(V));
        if (!input) throw std::runtime_error("cannot open seed file");
        std::size_t loaded = 0;
        // Seed sets are small. Keep only group streams, never all seed graphs.
        std::map<ValenceArray, std::ofstream> streams;
        for (std::string line; std::getline(input, line);) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            auto graph = read_seed<V>(line);
            auto degrees = graph.valence_array();
            std::sort(degrees.begin(), degrees.end());
            ValenceArray key(degrees.begin(), degrees.end());
            auto [it, inserted] = streams.try_emplace(key);
            if (inserted) {
                it->second.exceptions(std::ios::badbit | std::ios::failbit);
                it->second.open(directory / group_name(key), std::ios::binary);
            }
            write_record<V>(it->second, graph);
            ++groups[key];
            ++loaded;
        }
        if (!input.eof() || loaded != seed_counts_.at(V))
            throw std::runtime_error("seed count does not match counts.tsv");
        for (auto& [key, stream] : streams) stream.close();
        return groups;
    }

    template <int V>
    void process(const Groups& parents) {
        const auto started = std::chrono::steady_clock::now();
        const auto directory = stage_path(V);
        std::filesystem::create_directory(directory);
        const auto seed_directory = directory / "seeds";
        std::filesystem::create_directory(seed_directory);
        const auto seeds = partition_seeds<V>(seed_directory);
        struct Source {
            ValenceArray degrees;
            std::size_t count;
            std::vector<ValenceSplit> splits;
        };
        std::map<ValenceArray, std::vector<Source>> targets;
        for (const auto& [degrees, count] : parents)
            for (auto& [target, plan] : valence_split_plans(degrees))
                targets[target].push_back({degrees, count, std::move(plan)});
        for (const auto& [degrees, count] : seeds) targets.try_emplace(degrees);

        std::ofstream graphs;
        graphs.exceptions(std::ios::badbit | std::ios::failbit);
        graphs.open(output_ / ("graphs_L" + std::to_string(loop_number) + "_V" + std::to_string(V) + ".g6"));
        Groups completed;
        std::size_t candidates = 0, output_count = 0, stored_count = 0;
        for (const auto& [target, sources] : targets) {
            const auto group_started = std::chrono::steady_clock::now();
            Stage<V, true> destination;
            if (auto seed = seeds.find(target); seed != seeds.end())
                read_group<V>(seed_directory / group_name(target), seed->second,
                    [&](const auto& graph) { destination.add(graph); });
            if constexpr (V > first_vertices) {
                for (const auto& source : sources) {
                    read_group<V - 1>(stage_path(V - 1) / group_name(source.degrees), source.count,
                        [&](const SplitStageGraph<V - 1>& graph) {
                            const CutVertexSplitRule<V - 1> cut(graph);
                            std::optional<SplitReductionRule<SplitStageGraph<V - 1>>> reduction;
                            if (cut.vertex < 0) reduction.emplace(graph);
                            // Vertex indices and moved sizes were computed once,
                            // from the sorted parent/target arrays, before I/O.
                            int previous_vertex = -1;
                            bool skip = false;
                            for (const auto& split : source.splits) {
                                const int vertex = split.vertex;
                                if (vertex != previous_vertex) {
                                    previous_vertex = vertex;
                                    skip = (cut.vertex >= 0 && vertex != cut.vertex)
                                        || (reduction && reduction->skip_vertex(vertex));
                                }
                                if (skip) continue;
                                const auto adjacent = graph.adjacent(vertex);
                                const Int last = adjacent.size() - 1;
                                auto subset = combutils::firstSubset(1, split.moved);
                                do {
                                    std::uint64_t neighbours = 0;
                                    for (auto index : subset)
                                        neighbours |= std::uint64_t{1}
                                            << graph.half_edges[adjacent[index] ^ 1];
                                    if (cut.vertex >= 0 && !cut.resolves(neighbours)) continue;
                                    if (reduction && reduction->redundant(vertex, neighbours)) continue;
                                    destination.add(graph.splitGraph(vertex, adjacent, subset));
                                } while (combutils::nextSubset(subset, last));
                            }
                        });
                }
            }
            // Drain directly to disk: no second, packed copy of the target set.
            std::ofstream stored;
            stored.exceptions(std::ios::badbit | std::ios::failbit);
            if constexpr (V < max_vertices)
                stored.open(directory / group_name(target), std::ios::binary);
            std::size_t group_size = 0, group_output = 0;
            const auto group_candidates = destination.candidates();
            for (std::size_t b = 0; b < destination.bucket_count; ++b) {
                auto& set = destination.buckets[b].graphs;
                for (std::size_t i = 0; i < set.capacity(); ++i) {
                    const auto& entry = set.data()[i];
                    if (entry.empty()) continue;
                    const auto& graph = entry.graph;
                    #ifndef NDEBUG
                    const auto degrees = graph.valence_array();
                    assert(std::equal(degrees.begin(), degrees.end(), target.begin(), target.end()));
                    #endif
                    if constexpr (V < max_vertices) write_record<V>(stored, graph);
                    ++group_size;
                    if (CutVertexSplitRule<V>(graph).vertex < 0) {
                        write_graph6<V>(graphs, graph);
                        ++group_output;
                    }
                }
                set = {};
            }
            if constexpr (V < max_vertices) stored.close();
            if (group_size) completed.emplace(target, group_size);
            candidates += group_candidates;
            output_count += group_output;
            stored_count += group_size;
            const auto seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - group_started).count();
            group_counts_ << V << '\t' << group_name(target) << '\t' << sources.size()
                << '\t' << group_candidates << '\t' << group_size << '\t' << group_output
                << '\t' << seconds << std::endl;
            std::cout << "V=" << V << " valences=" << group_name(target)
                << " graphs=" << group_size << " seconds=" << seconds << std::endl;
        }
        graphs.close();
        std::filesystem::remove_all(seed_directory);
        if constexpr (V > first_vertices) std::filesystem::remove_all(stage_path(V - 1));
        counts_ << loop_number << '\t' << V << '\t' << V - 1 + loop_number << '\t'
            << (seed_counts_.contains(V) ? seed_counts_.at(V) : 0) << '\t'
            << candidates << '\t' << output_count << '\t'
            << std::chrono::duration<double>(std::chrono::steady_clock::now() - started).count()
            << std::endl;
        std::cout << "Finished V=" << V << " graphs=" << output_count
                  << " pending_cut_graphs=" << stored_count - output_count << std::endl;
        if constexpr (V < max_vertices) process<V + 1>(completed);
    }

public:
    void run(const std::filesystem::path& input, const std::filesystem::path& output) {
        input_ = input;
        output_ = output;
        seed_counts_ = read_triangle_seed_counts(input_);
        if (!output_.parent_path().empty()) std::filesystem::create_directories(output_.parent_path());
        if (!std::filesystem::create_directory(output_))
            throw std::runtime_error("use a fresh output directory: " + output_.string());
        scratch_ = output_ / "working";
        std::filesystem::create_directory(scratch_);
        counts_.exceptions(std::ios::badbit | std::ios::failbit);
        counts_.open(output_ / "counts.tsv");
        counts_ << "loop\tvertices\tedges\tseed_records\tcandidates\tgraphs\tseconds\n";
        group_counts_.exceptions(std::ios::badbit | std::ios::failbit);
        group_counts_.open(output_ / "valence_counts.tsv");
        group_counts_ << "vertices\tvalences\tparent_groups\tcandidates\tstored_graphs\tgraphs\tseconds\n";
        process<first_vertices>({});
        counts_.close();
        group_counts_.close();
        std::filesystem::remove_all(scratch_);
    }
};
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " TRIANGLE_SEED_DIRECTORY NEW_OUTPUT_DIRECTORY\n";
        return 2;
    }
    try {
        std::cout << "Valence-group splitting: loop=" << loop_number
            << " workers=" << omp_get_max_threads() << std::endl;
        Pipeline{}.run(argv[1], argv[2]);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
