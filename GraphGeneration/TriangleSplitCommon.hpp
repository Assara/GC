#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifdef GC_TRIANGLE_PARALLEL_SPLITS
#include <atomic>
#include <exception>
#include <omp.h>
#include <mutex>
#endif

#include "graph.hpp"
#include "GraphStandardizer.hpp"
#include "LinearProbeSet.hpp"
#include "GraphGeneration/CutVertexSplitRule.hpp"
#include "GraphGeneration/SplitReductionRule.hpp"
#include "GraphGeneration/SplitStageEstimate.hpp"
#include "GraphGeneration/SplitReservation.hpp"
#include "GraphGeneration/ValenceGraphStandardizer.hpp"

#ifndef GC_TRIANGLE_SPLIT_LOOP
#define GC_TRIANGLE_SPLIT_LOOP 6
#endif

namespace GraphGeneration::triangle_split_detail {
constexpr int loop_number = GC_TRIANGLE_SPLIT_LOOP;
constexpr int max_vertices = 2 * (loop_number - 1); // Minimum degree three.
constexpr int first_vertices = [] {
    int v = 4;
    while (v * (v - 1) / 2 < v - 1 + loop_number) ++v;
    return v;
}();
static_assert(loop_number >= 3 && max_vertices <= 62);

template <int V>
using SplitStageGraph = Graph<V, V - 1 + loop_number, 0, 0, 0, 0, fieldType>;

template <int V, bool IncreasingValency = false>
struct Stage {
    struct Entry {
        SplitStageGraph<V> graph;
        Entry() noexcept : graph{} {}
        explicit Entry(const SplitStageGraph<V>& value) noexcept : graph(value) {}
        // Simple graphs have distinct endpoints on every edge. The default
        // all-zero graph is therefore an unused slot, with no extra flag.
        bool empty() const noexcept { return graph.half_edges[0] == graph.half_edges[1]; }
        std::size_t hash() const noexcept { return graph.hash(); }
        bool operator==(const Entry& other) const noexcept {
            return graph.half_edges == other.graph.half_edges;
        }
    };
    // Sharding bounds resize peaks and lets workers deduplicate immediately.
    static constexpr std::size_t bucket_count = 64;
    struct Bucket {
        linear_probe_set<Entry> graphs;
        std::size_t candidates = 0;
#ifdef GC_TRIANGLE_PARALLEL_SPLITS
        std::mutex mutex;
#endif
    };
    std::unique_ptr<Bucket[]> buckets = std::make_unique<Bucket[]>(bucket_count);

    // Called only on a fresh destination stage, before any workers start.
    std::size_t reserve(std::size_t expected, std::size_t budget) {
        expected = GraphGeneration::step_down_split_reservation<Entry>(expected, bucket_count, budget);
        while (expected) {
            const auto per_bucket = expected / bucket_count + (expected % bucket_count != 0);
            std::size_t b = 0;
            for (; b < bucket_count; ++b)
                if (!buckets[b].graphs.try_reserve(per_bucket)) break;
            if (b == bucket_count) return expected;
            // Release a partial reservation before retrying a smaller one.
            for (std::size_t i = 0; i < bucket_count; ++i) buckets[i].graphs = {};
            expected /= 2;
        }
        return 0;
    }

    std::size_t candidates() const {
        std::size_t count = 0;
        for (std::size_t i = 0; i < bucket_count; ++i) count += buckets[i].candidates;
        return count;
    }

    void add(const SplitStageGraph<V>& graph) {
        GraphStandardizer<V, V - 1 + loop_number, 0, 0, 0, 0, fieldType> standardizer;
        const Entry entry([&] {
            if constexpr (IncreasingValency)
                return GraphGeneration::standardize_sorted_valences(graph, standardizer);
            else return standardizer.standardize_no_sign(graph);
        }());
        // Use high hash bits: the probe table uses low bits for its slot.
        auto& bucket = buckets[entry.hash() >> (sizeof(std::size_t) * 8 - 6)];
#ifdef GC_TRIANGLE_PARALLEL_SPLITS
        const std::lock_guard lock(bucket.mutex);
#endif
        ++bucket.candidates;
        bucket.graphs.insert(entry);
    }

    std::vector<SplitStageGraph<V>> take_sorted() {
        std::vector<SplitStageGraph<V>> result;
        std::size_t size = 0;
        for (std::size_t i = 0; i < bucket_count; ++i) size += buckets[i].graphs.size();
        result.reserve(size);
        for (std::size_t b = 0; b < bucket_count; ++b) {
            auto& graphs = buckets[b].graphs;
            for (std::size_t i = 0; i < graphs.capacity(); ++i)
                if (!graphs.data()[i].empty()) result.push_back(graphs.data()[i].graph);
            graphs = {};
        }
        std::sort(result.begin(), result.end(), [](const auto& a, const auto& b) {
            return a.half_edges < b.half_edges;
        });
        return result;
    }
};

template <int V>
SplitStageGraph<V> read_seed(const std::string& line) {
    constexpr int edges = V - 1 + loop_number;
    constexpr int bits = V * (V - 1) / 2;
    if (line.size() != 1 + (bits + 5) / 6 || line[0] != V + 63)
        throw std::runtime_error("invalid seed graph6 dimensions");
    for (unsigned char c : line)
        if (c < 63 || c > 126) throw std::runtime_error("invalid graph6 character");
    SplitStageGraph<V> graph;
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
    return graph;
}

template <int V>
void write_graph6(std::ostream& out, const SplitStageGraph<V>& graph) {
    std::array<std::array<bool, V>, V> adjacent{};
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

} // namespace GraphGeneration::triangle_split_detail
