#include <algorithm>
#include <array>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef GC_TRIANGLE_KEY_STANDARDIZER
#include "GraphGeneration/TriangleSeedKeyStandardizer.hpp"
#else
#include "GraphGeneration/TriangleSeedStandardizer.hpp"
#endif
#include "GraphGeneration/TriangleComponentStack.hpp"
#include "GraphGeneration/TriangleCompletionBudget.hpp"
#include "GraphGeneration/TriangleSeedEntry.hpp"
#include "LinearProbeSet.hpp"

#ifndef GC_TRIANGLE_MAX_LOOP
#define GC_TRIANGLE_MAX_LOOP 6
#endif

namespace {
constexpr int max_loop = GC_TRIANGLE_MAX_LOOP;
// A triangle attached at one vertex adds two vertices and one loop.
constexpr int max_vertices = 2 * max_loop + 1;
static_assert(max_loop >= 1);
static_assert(max_vertices <= 62); // Single-byte graph6 vertex count.
static_assert(6 * max_loop <= std::numeric_limits<Int>::max());

using Edges = std::vector<Int>; // Consecutive endpoint pairs.
using Adjacency = std::vector<std::vector<bool>>;

Adjacency adjacency(int vertices, const Edges& edges) {
    Adjacency result(vertices, std::vector<bool>(vertices));
    for (std::size_t i = 0; i < edges.size(); i += 2)
        result[edges[i]][edges[i + 1]] = result[edges[i + 1]][edges[i]] = true;
    return result;
}

// Every vertex in a triangle-grown graph has degree at least two.
std::uint64_t bivalent_vertices(int vertices, const Edges& edges) {
    std::array<int, max_vertices> degrees{};
    for (Int v : edges) ++degrees[v];
    std::uint64_t mask = 0;
    for (int v = 0; v < vertices; ++v)
        if (degrees[v] == 2) mask |= std::uint64_t{1} << v;
    return mask;
}

void write_graph6(std::ostream& out, int vertices, const Edges& edges) {
    const auto adjacent = adjacency(vertices, edges);
    out.put(static_cast<char>(vertices + 63));
    int value = 0, bits = 0;
    for (int b = 1; b < vertices; ++b)
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

#ifdef GC_TRIANGLE_KEY_STANDARDIZER
struct StoredGraph : GraphGeneration::TriangleSeedEntry<max_loop> {
    std::array<Int, 6 * max_loop> original{};
};
#else
using StoredGraph = GraphGeneration::TriangleSeedEntry<max_loop>;
#endif
class Generator {
    struct LiveGraph {
        Edges edges;
        // Components are relabelled with the graph, preserving stack order.
        GraphGeneration::TriangleComponentStack<max_loop> components;
    };
    struct Bucket {
        linear_probe_set<StoredGraph> graphs;
        std::size_t candidates = 0;
    };
    // Every augmentation strictly increases loop order, including those with no new vertex.
    std::map<std::pair<int, int>, Bucket> pending_; // (loop order, vertex count)

    void add(int vertices, LiveGraph graph) {
        const auto& edges = graph.edges;
        const int loop = static_cast<int>(edges.size() / 2) - vertices + 1;
        if (loop > max_loop || vertices > max_vertices) return;
        const auto bivalent = bivalent_vertices(vertices, edges);
        if (!GraphGeneration::triangle_completion_fits(
                vertices, std::popcount(bivalent), loop, max_loop)) return;
        if (bivalent == 0) graph.components.forget();
#ifdef GC_TRIANGLE_KEY_STANDARDIZER
        if (graph.components.reversed_sizes_are_greater()) return;
#else
        if (graph.components.last_exceeds_first()) return;
#endif
        // Forgetting the construction stack does not remove actual cut vertices.
#ifdef GC_TRIANGLE_KEY_STANDARDIZER
        const auto key = GraphGeneration::TriangleSeedKeyStandardizer<max_loop>{}
            .canonical_key(vertices, edges);
        StoredGraph entry;
        entry.endpoint_count = edges.size();
        std::copy(key.begin(), key.end(), entry.canonical.begin());
        std::copy(edges.begin(), edges.end(), entry.original.begin());
        entry.components = graph.components;
#else
        const auto entry = GraphGeneration::TriangleSeedStandardizer<max_loop>{}
            .standardize(vertices, edges, graph.components);
        if (entry.empty()) return; // Noncanonical component orientation.
#endif
        auto& bucket = pending_[{loop, vertices}];
        ++bucket.candidates;
        bucket.graphs.insert(entry);
    }

    void expand(int vertices, const LiveGraph& graph) {
        const auto& edges = graph.edges;
        // New labels are consecutive: vertices, then vertices + 1.
        // At least one existing vertex keeps every intermediate graph connected.
        const auto adjacent = adjacency(vertices + 2, edges);
        const int limit = std::min(vertices + 2, max_vertices);
        for (int a = 0; a < vertices; ++a)
            for (int b = a + 1; b < limit; ++b)
                for (int c = b + 1; c < limit; ++c) {
                    if (c == vertices + 1 && b != vertices) continue;
                    std::uint64_t existing = 0, introduced = 0;
                    for (int v : {a, b, c}) {
                        if (v < vertices) existing |= std::uint64_t{1} << v;
                        else introduced |= std::uint64_t{1} << v;
                    }
                    if (!graph.components.touches(existing)) continue;
                    const std::array<std::pair<int, int>, 3> triangle{{{a, b}, {a, c}, {b, c}}};
                    const int added_vertices = (b >= vertices) + (c >= vertices);
                    int added_edges = 0;
                    for (auto [u, v] : triangle) added_edges += !adjacent[u][v];
                    if (added_edges == 0) continue;
                    const int child_loop = static_cast<int>(edges.size() / 2) - vertices + 1
                        + added_edges - added_vertices;
                    if (child_loop > max_loop) continue;
                    LiveGraph child = graph;
                    child.components.add_triangle(existing, introduced, added_vertices);
                    for (auto [u, v] : triangle)
                        if (!adjacent[u][v]) {
                            child.edges.push_back(static_cast<Int>(u));
                            child.edges.push_back(static_cast<Int>(v));
                        }
                    add(vertices + added_vertices, std::move(child));
                }
    }

public:
    void run(const std::filesystem::path& directory) {
        // Refuse existing directories so previous runs cannot be mixed or overwritten.
        if (!directory.parent_path().empty())
            std::filesystem::create_directories(directory.parent_path());
        if (!std::filesystem::create_directory(directory))
            throw std::runtime_error("use a fresh output directory: " + directory.string());
        std::ofstream counts(directory / "counts.tsv");
        counts.exceptions(std::ios::badbit | std::ios::failbit);
        counts << "loop\tvertices\tedges\tcandidates\ttriangle_graphs\tseed_graphs\tseconds\n";
        const auto started = std::chrono::steady_clock::now();
        auto last_progress = started;
        add(3, {{0, 1, 0, 2, 1, 2}, {}});
        while (!pending_.empty()) {
            auto current = pending_.extract(pending_.begin());
            const auto [loop, vertices] = current.key();
            const auto& bucket = current.mapped();
            const int edges = vertices - 1 + loop;
            std::cout << "Starting L=" << loop << " V=" << vertices
                      << " graphs=" << bucket.graphs.size() << std::endl;
            const auto stage_start = std::chrono::steady_clock::now();
            const std::string suffix = "_L" + std::to_string(loop) + "_V" + std::to_string(vertices) + ".g6";
            std::ofstream all(directory / ("triangles" + suffix));
            std::ofstream seeds(directory / ("seeds" + suffix));
            all.exceptions(std::ios::badbit | std::ios::failbit);
            seeds.exceptions(std::ios::badbit | std::ios::failbit);
            std::size_t seed_count = 0, processed = 0;
            for (std::size_t slot = 0; slot < bucket.graphs.capacity(); ++slot) {
                const auto& entry = bucket.graphs.data()[slot];
                if (entry.empty()) continue;
                const Edges canonical(entry.canonical.begin(),
                                      entry.canonical.begin() + entry.endpoint_count);
#ifdef GC_TRIANGLE_KEY_STANDARDIZER
                const LiveGraph graph{
                    Edges(entry.original.begin(), entry.original.begin() + entry.endpoint_count),
                    entry.components};
#else
                const LiveGraph graph{canonical, entry.components};
#endif
                write_graph6(all, vertices, canonical);
                std::vector<int> degrees(vertices);
                for (Int v : canonical) ++degrees[v];
                if (std::ranges::all_of(degrees, [](int d) { return d >= 3; })) {
                    write_graph6(seeds, vertices, canonical);
                    ++seed_count;
                }
                if (loop < max_loop) expand(vertices, graph);
                ++processed;
                const auto now = std::chrono::steady_clock::now();
                if (now - last_progress >= std::chrono::seconds(5)) {
                    std::cout << "Expanding L=" << loop << " V=" << vertices
                              << " parents=" << processed << '/' << bucket.graphs.size() << std::endl;
                    last_progress = now;
                }
            }
            all.close();
            seeds.close();
            const double seconds = std::chrono::duration<double>(
                std::chrono::steady_clock::now() - stage_start).count();
            counts << loop << '\t' << vertices << '\t' << edges << '\t' << bucket.candidates
                   << '\t' << bucket.graphs.size() << '\t' << seed_count << '\t' << seconds << std::endl;
            std::cout << "Finished L=" << loop << " V=" << vertices
                      << " candidates=" << bucket.candidates << " triangle_graphs=" << bucket.graphs.size()
                      << " seed_graphs=" << seed_count << " seconds=" << seconds << std::endl;
        }
        counts.close();
        std::cout << "Total seconds=" << std::chrono::duration<double>(
            std::chrono::steady_clock::now() - started).count() << std::endl;
    }
};
} // namespace

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " NEW_OUTPUT_DIRECTORY\n";
        return 2;
    }
    try {
#ifdef GC_TRIANGLE_KEY_STANDARDIZER
        std::cout << "Standardizer: full-graph key, original working graph; first/last rule\n";
#else
        std::cout << "Standardizer: custom component-aware; ordered stack, first/last rule\n";
#endif
        std::cout << "Triangle seed proof of concept: K3, max_loop=" << max_loop
                  << " max_vertices=" << max_vertices << std::endl;
        Generator{}.run(argv[1]);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
