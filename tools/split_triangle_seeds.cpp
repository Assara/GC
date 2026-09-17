#include "GraphGeneration/TriangleSplitCommon.hpp"

namespace {
using namespace GraphGeneration::triangle_split_detail;

class Pipeline {
    std::filesystem::path input_, output_;
    std::map<int, std::size_t> seed_counts_;
    std::ofstream counts_;
    GraphGeneration::SplitStageEstimator estimator_;

    std::filesystem::path seed_path(int vertices) const {
        return triangle_seed_path(input_, vertices);
    }

    template <int V>
    void process(Stage<V> stage) {
        const auto started = std::chrono::steady_clock::now();
        auto last_progress = started;
        std::size_t loaded = 0;
        if (seed_counts_.contains(V)) {
            std::ifstream input(seed_path(V));
            if (!input) throw std::runtime_error("cannot open seed file");
            for (std::string line; std::getline(input, line);) {
                if (!line.empty() && line.back() == '\r') line.pop_back();
                stage.add(read_seed<V>(line));
                ++loaded;
            }
            if (!input.eof() || loaded != seed_counts_.at(V))
                throw std::runtime_error("seed file count does not match counts.tsv");
        }
        auto parents = stage.take_sorted();
        std::cout << "Starting L=" << loop_number << " V=" << V
                  << " seeds=" << loaded << " graphs=" << parents.size() << std::endl;
        std::ofstream graphs(output_ / ("graphs_L" + std::to_string(loop_number)
            + "_V" + std::to_string(V) + ".g6"));
        graphs.exceptions(std::ios::badbit | std::ios::failbit);
        std::size_t output_count = 0;
        Stage<V + 1> next;
        if constexpr (V < max_vertices) {
            const auto estimate = estimator_.next(loop_number, V, parents.size());
            const auto budget = GraphGeneration::split_reservation_budget();
            const auto reserved = next.reserve(estimate, budget);
            std::cout << "Reserving V=" << V + 1 << " estimated_graphs=" << estimate
                      << " reserved_graphs=" << reserved << " budget_bytes=" << budget << std::endl;
            if (reserved < estimate)
                std::cout << "Reservation stepped down to fit memory headroom" << std::endl;
        }
        const auto expand = [&](const SplitStageGraph<V>& graph, Stage<V + 1>& destination) {
            const GraphGeneration::CutVertexSplitRule<V> cut(graph);
            if constexpr (V < max_vertices) {
                const auto valences = graph.valence_array();
                std::optional<GraphGeneration::SplitReductionRule<SplitStageGraph<V>>> reduction;
                if (cut.vertex < 0) reduction.emplace(graph);
                for (int vertex = V - 1; vertex >= 0; --vertex) {
                    // Resolve one cut vertex completely before any other split.
                    // Recompute on the child to resolve any remaining cuts.
                    if (cut.vertex >= 0 && vertex != cut.vertex) continue;
                    if (valences[vertex] <= 3) continue;
                    if (reduction && reduction->skip_vertex(vertex)) continue;
                    const auto adjacent = graph.adjacent(vertex);
                    const Int last = adjacent.size() - 1;
                    // Keep incidence zero on the old vertex to identify
                    // complementary splits. Each side needs two incidences
                    // plus the connecting edge to reach valence three.
                    for (Int moved = 2; moved < last; ++moved) {
                        auto subset = combutils::firstSubset(1, moved);
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
                }
            }
            return cut.vertex < 0;
        };
#ifdef GC_TRIANGLE_PARALLEL_SPLITS
        if constexpr (V < max_vertices) {
            // Stable parents and separate output flags preserve serial file order.
            std::vector<unsigned char> output_flags(parents.size());
            std::atomic<bool> failed{false};
            std::atomic<std::size_t> processed{0};
            std::exception_ptr error;
            #pragma omp parallel
            {
                const int worker = omp_get_thread_num();
                #pragma omp for schedule(dynamic, 8)
                for (std::size_t i = 0; i < parents.size(); ++i) {
                    if (failed.load(std::memory_order_relaxed)) continue;
                    try {
                        output_flags[i] = expand(parents[i], next);
                        const auto done = processed.fetch_add(1, std::memory_order_relaxed) + 1;
                        if (worker == 0) {
                            const auto now = std::chrono::steady_clock::now();
                            if (now - last_progress >= std::chrono::seconds(5)) {
                                std::cout << "Expanding V=" << V << " parents=" << done
                                          << '/' << parents.size() << std::endl;
                                last_progress = now;
                            }
                        }
                    } catch (...) {
                        failed.store(true, std::memory_order_relaxed);
                        #pragma omp critical(triangle_split_error)
                        {
                            if (!error) error = std::current_exception();
                        }
                    }
                }
            }
            if (error) std::rethrow_exception(error);
            for (std::size_t i = 0; i < parents.size(); ++i) {
                if (!output_flags[i]) continue;
                write_graph6<V>(graphs, parents[i]);
                ++output_count;
            }
        } else
#endif
        {
            std::size_t processed = 0;
            for (const auto& graph : parents) {
                if (expand(graph, next)) {
                    write_graph6<V>(graphs, graph);
                    ++output_count;
                }
                ++processed;
                const auto now = std::chrono::steady_clock::now();
                if (now - last_progress >= std::chrono::seconds(5)) {
                    std::cout << "Expanding V=" << V << " parents=" << processed
                              << '/' << parents.size() << std::endl;
                    last_progress = now;
                }
            }
        }
        graphs.close();
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        counts_ << loop_number << '\t' << V << '\t' << V - 1 + loop_number
                << '\t' << loaded << '\t' << stage.candidates() << '\t' << output_count
                << '\t' << seconds << std::endl;
        std::cout << "Finished V=" << V << " candidates=" << stage.candidates()
                  << " graphs=" << output_count << " pending_cut_graphs=" << parents.size() - output_count
                  << " seconds=" << seconds << std::endl;
        std::vector<SplitStageGraph<V>>{}.swap(parents);
        if constexpr (V < max_vertices) process<V + 1>(std::move(next));
    }

public:
    void run(const std::filesystem::path& input, const std::filesystem::path& output) {
        estimator_ = {};
        input_ = input;
        output_ = output;
        seed_counts_ = read_triangle_seed_counts(input_);
        if (!output_.parent_path().empty()) std::filesystem::create_directories(output_.parent_path());
        if (!std::filesystem::create_directory(output_))
            throw std::runtime_error("use a fresh output directory: " + output_.string());
        counts_.open(output_ / "counts.tsv");
        counts_.exceptions(std::ios::badbit | std::ios::failbit);
        counts_ << "loop\tvertices\tedges\tseed_records\tcandidates\tgraphs\tseconds\n";
        const auto started = std::chrono::steady_clock::now();
        process<first_vertices>({});
        counts_.close();
        std::cout << "Total seconds=" << std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count() << std::endl;
    }
};
} // namespace

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " TRIANGLE_SEED_DIRECTORY NEW_OUTPUT_DIRECTORY\n";
        return 2;
    }
    try {
        std::cout << "Triangle seed splitting: loop=" << loop_number
                  << " max_vertices=" << max_vertices << "; all vertices of valence > 3, no duplication\n";
#ifdef GC_TRIANGLE_PARALLEL_SPLITS
        std::cout << "OpenMP workers=" << omp_get_max_threads() << std::endl;
#endif
        Pipeline{}.run(argv[1], argv[2]);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
