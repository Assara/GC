#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include "GraphGeneration/TransientGraph2Standardizer.hpp"

#ifndef GC_TRIANGLE_SPLIT_LOOP
#define GC_TRIANGLE_SPLIT_LOOP 6
#endif

namespace {
constexpr int loop_number = GC_TRIANGLE_SPLIT_LOOP;
constexpr int max_vertices = 2 * (loop_number - 1); // Minimum degree three.
constexpr int first_vertices = [] {
    int v = 4;
    while (v * (v - 1) / 2 < v - 1 + loop_number) ++v;
    return v;
}();
static_assert(loop_number >= 3 && max_vertices <= 62);

template <int V>
using Transient = GraphGeneration::transient_graph2<V, V - 1 + loop_number>;

template <int V>
struct Stage {
    struct Less {
        bool operator()(const Transient<V>& a, const Transient<V>& b) const {
            return a.graph().half_edges < b.graph().half_edges;
        }
    };
    std::set<Transient<V>, Less> graphs;
    std::size_t candidates = 0;

    void add(Transient<V> graph) {
        ++candidates;
        GraphGeneration::transient_graph2_standardizer<V, V - 1 + loop_number> standardizer;
        graphs.insert(standardizer.standardize_no_sign(graph));
    }
};

template <int V>
struct Collector {
    Stage<V>& stage;
    template <typename G>
    void add(G graph) {
        if constexpr (G::N_VERTICES_ == V && G::N_EDGES_ == V - 1 + loop_number)
            stage.add(Transient<V>(std::move(graph)));
        else
            throw std::logic_error("split changed loop order");
    }
};

template <int V>
Transient<V> read_seed(const std::string& line) {
    constexpr int edges = V - 1 + loop_number;
    constexpr int bits = V * (V - 1) / 2;
    if (line.size() != 1 + (bits + 5) / 6 || line[0] != V + 63)
        throw std::runtime_error("invalid seed graph6 dimensions");
    for (unsigned char c : line)
        if (c < 63 || c > 126) throw std::runtime_error("invalid graph6 character");
    typename Transient<V>::graph_type graph;
    std::array<std::array<bool, V>, V> adjacent{};
    int bit = 0, e = 0;
    for (int b = 1; b < V; ++b)
        for (int a = 0; a < b; ++a, ++bit)
            if (((line[1 + bit / 6] - 63) >> (5 - bit % 6)) & 1) {
                if (e == edges) throw std::runtime_error("seed has too many edges");
                graph.setEdge(e++, a, b);
                adjacent[a][b] = adjacent[b][a] = true;
            }
    if (e != edges) throw std::runtime_error("seed has wrong loop order");
    const auto valences = graph.valence_array();
    if (std::ranges::any_of(valences, [](Int d) { return d < 3; }))
        throw std::runtime_error("use the minimum-degree-three seed files");
    for (int a = 0; a < V; ++a)
        for (int b = a + 1; b < V; ++b) {
            if (!adjacent[a][b]) continue;
            bool covered = false;
            for (int c = 0; c < V; ++c) covered |= adjacent[a][c] && adjacent[b][c];
            if (!covered) throw std::runtime_error("seed edge is not covered by a triangle");
        }
    return Transient<V>(std::move(graph));
}

template <int V>
void write_graph6(std::ostream& out, const Transient<V>& transient) {
    std::array<std::array<bool, V>, V> adjacent{};
    const auto& graph = transient.graph();
    for (int e = 0; e < graph.N_EDGES_; ++e) {
        const auto [a, b] = graph.getEdge(e);
        adjacent[a][b] = adjacent[b][a] = true;
    }
    out.put(static_cast<char>(V + 63));
    int value = 0, bits = 0;
    for (int b = 1; b < V; ++b)
        for (int a = 0; a < b; ++a) {
            value = (value << 1) | adjacent[a][b];
            if (++bits == 6) {
                out.put(static_cast<char>(value + 63));
                value = bits = 0;
            }
        }
    if (bits) out.put(static_cast<char>((value << (6 - bits)) + 63));
    out.put('\n');
}

class Pipeline {
    std::filesystem::path input_, output_;
    std::map<int, std::size_t> seed_counts_;
    std::ofstream counts_;

    std::filesystem::path seed_path(int vertices) const {
        return input_ / ("seeds_L" + std::to_string(loop_number)
            + "_V" + std::to_string(vertices) + ".g6");
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
        std::cout << "Starting L=" << loop_number << " V=" << V
                  << " seeds=" << loaded << " graphs=" << stage.graphs.size() << std::endl;
        std::ofstream graphs(output_ / ("graphs_L" + std::to_string(loop_number)
            + "_V" + std::to_string(V) + ".g6"));
        graphs.exceptions(std::ios::badbit | std::ios::failbit);
        for (const auto& graph : stage.graphs) write_graph6<V>(graphs, graph);
        graphs.close();

        Stage<V + 1> next;
        if constexpr (V < max_vertices) {
            Collector<V + 1> collector{next};
            std::size_t processed = 0;
            for (const auto& graph : stage.graphs) {
                const auto valences = graph.graph().valence_array();
                assert(std::is_sorted(valences.begin(), valences.end()));
                const Int maximum = valences.back();
                for (int vertex = V - 1; vertex >= 0 && valences[vertex] > 3; --vertex)
                    // Equal loop bounds prohibit every shared incidence.
                    graph.split(vertex, 0, 3, maximum, loop_number, loop_number, collector);
                ++processed;
                const auto now = std::chrono::steady_clock::now();
                if (now - last_progress >= std::chrono::seconds(5)) {
                    std::cout << "Expanding V=" << V << " parents=" << processed
                              << '/' << stage.graphs.size() << std::endl;
                    last_progress = now;
                }
            }
        }
        const double seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count();
        counts_ << loop_number << '\t' << V << '\t' << V - 1 + loop_number
                << '\t' << loaded << '\t' << stage.candidates << '\t' << stage.graphs.size()
                << '\t' << seconds << std::endl;
        std::cout << "Finished V=" << V << " candidates=" << stage.candidates
                  << " graphs=" << stage.graphs.size() << " seconds=" << seconds << std::endl;
        stage.graphs.clear();
        if constexpr (V < max_vertices) process<V + 1>(std::move(next));
    }

public:
    void run(const std::filesystem::path& input, const std::filesystem::path& output) {
        input_ = input;
        output_ = output;
        std::ifstream manifest(input_ / "counts.tsv");
        std::string line;
        if (!std::getline(manifest, line)
            || (line != "loop\tvertices\tedges\tcandidates\ttriangle_graphs\tseed_graphs\tseconds"
                && line != "loop\tvertices\tedges\tcandidates\ttriangle_graphs\tmin_degree_3\tseconds"))
            throw std::runtime_error("expected triangle seed generator counts.tsv");
        while (std::getline(manifest, line)) {
            int loop, vertices, edges;
            std::size_t candidates, total, seeds;
            double seconds;
            std::istringstream row(line);
            if (!(row >> loop >> vertices >> edges >> candidates >> total >> seeds >> seconds))
                throw std::runtime_error("invalid triangle seed count row");
            if (loop == loop_number && vertices >= first_vertices && vertices <= max_vertices) {
                if (edges != vertices - 1 + loop || seeds > total
                    || !seed_counts_.emplace(vertices, seeds).second)
                    throw std::runtime_error("inconsistent triangle seed counts");
            }
        }
        if (!manifest.eof()) throw std::runtime_error("failed reading seed counts");
        for (int v = first_vertices; v <= std::min(max_vertices, loop_number + 2); ++v)
            if (!seed_counts_.contains(v) || !std::filesystem::is_regular_file(seed_path(v)))
                throw std::runtime_error("missing triangle seed stage V=" + std::to_string(v));
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
        Pipeline{}.run(argv[1], argv[2]);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
