#include "GraphGeneration/GraphGenerationPipeline.hpp"
#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>

using namespace GraphGeneration;
using Pipeline = GraphGenerationPipeline<3, 6>;

template <typename F>
void must_throw(F&& action) {
    bool threw = false;
    try { action(); } catch (const std::exception&) { threw = true; }
    assert(threw);
}

std::vector<char> bytes(const std::filesystem::path& path) {
    std::ifstream input(path, std::ios::binary);
    return {std::istreambuf_iterator<char>(input), {}};
}

template <typename G>
void check_canonicalization(const G& graph) {
    using T = transient_graph2<G::N_VERTICES_, G::N_EDGES_>;
    transient_graph2_standardizer<G::N_VERTICES_, G::N_EDGES_> standardizer;
    const auto canonical = standardizer.standardize_no_sign(T(graph));
    const auto valences = canonical.graph().valence_array();
    assert(std::is_sorted(valences.begin(), valences.end()));
    assert(standardizer.standardize_no_sign(canonical).graph() == canonical.graph());
    Permutation<G::N_VERTICES_> permutation;
    do {
        G permuted;
        permuted.assignPermutedDirectedSortedEdgesNoSign(graph, permutation);
        assert(standardizer.standardize_no_sign(T(permuted)).graph() == canonical.graph());
    } while (std::next_permutation(permutation.p.begin(), permutation.p.end()));
}

template <std::size_t V, std::size_t L>
std::size_t check_file(const std::filesystem::path& directory) {
    using G = Graph<V, V - 1 + L, 0, 0, 0, 0, fieldType>;
    const auto path = Pipeline::file_path(directory, V, G::N_EDGES_);
    MappedTransientGraph2Reader<G> reader(path);
    assert(reader.header().record_size_bytes == 2 * G::N_EDGES_);
    assert(reader.header().graph_count == reader.header().splittable_graph_count);
    const auto contents = bytes(path);
    assert(contents.size() == 64 + reader.size() * 2 * G::N_EDGES_);
    std::set<std::array<Int, G::SIZE>> unique;
    transient_graph2_standardizer<V, G::N_EDGES_> standardizer;
    for (std::size_t i = 0; i < reader.size(); ++i) {
        const G graph = reader[i];
        assert(unique.insert(graph.half_edges).second);
        assert(standardizer.standardize_no_sign(transient_graph2<V, G::N_EDGES_>(graph)).graph() == graph);
        const auto valences = graph.valence_array();
        assert(std::is_sorted(valences.begin(), valences.end()));
        assert(valences.front() >= 1);
        std::set<std::pair<Int, Int>> edges;
        for (Int e = 0; e < G::N_EDGES_; ++e) {
            const auto [a, b] = graph.getEdge(e);
            assert(a < b && b < V);
            assert(edges.emplace(a, b).second);
            assert(static_cast<unsigned char>(contents[64 + i * G::SIZE + 2 * e]) == a);
            assert(static_cast<unsigned char>(contents[64 + i * G::SIZE + 2 * e + 1]) == b);
        }
        // One and the same maximum-valence vertex must support all leaves
        // and bivalent vertices together, excluding itself as in the seed.
        bool has_common_maximum = false;
        for (Int candidate = 0; candidate < V; ++candidate) {
            if (valences[candidate] != valences.back()) continue;
            bool supports_all = true;
            for (Int vertex = 0; vertex < V; ++vertex)
                if (vertex != candidate && valences[vertex] <= 2)
                    supports_all &= edges.contains(std::minmax(candidate, vertex));
            has_common_maximum |= supports_all;
        }
        assert(has_common_maximum);
    }
    must_throw([&] { (void)reader[reader.size()]; });
    return reader.size();
}

template <std::size_t V>
std::size_t check_stage(const std::filesystem::path& directory) {
    constexpr auto last_loop = std::min<std::size_t>(3, (V - 1) * (V - 2) / 2);
    return [&]<std::size_t... L>(std::index_sequence<L...>) {
        return (check_file<V, L>(directory) + ...);
    }(std::make_index_sequence<last_loop + 1>{});
}

int main() {
    using Star = Graph<5, 4, 0, 0, 0, 0, fieldType>;
    Star star;
    for (Int e = 0; e < 4; ++e) star.setEdge(e, 0, e + 1);
    check_canonicalization(star);
    using Diamond = Graph<4, 5, 0, 0, 0, 0, fieldType>;
    Diamond diamond;
    diamond.setEdge(0, 0, 1); diamond.setEdge(1, 0, 2);
    diamond.setEdge(2, 0, 3); diamond.setEdge(3, 1, 2); diamond.setEdge(4, 1, 3);
    check_canonicalization(diamond);

    char name[] = "/tmp/gc-new-pipeline-test-XXXXXX";
    const char* created = ::mkdtemp(name);
    assert(created);
    const std::filesystem::path directory(created);
    Pipeline pipeline;
    const auto summary = pipeline.run(directory / "first");
    assert(summary.size() == 4);
    assert(summary[0].unique_graphs == 1); // K3
    assert(summary[1].unique_graphs == 3); // Triangle with a leaf, diamond, K4.
    assert(summary[1].candidates > summary[1].unique_graphs);
    assert(check_stage<3>(directory / "first") == summary[0].unique_graphs);
    assert(check_stage<4>(directory / "first") == summary[1].unique_graphs);
    assert(check_stage<5>(directory / "first") == summary[2].unique_graphs);
    assert(check_stage<6>(directory / "first") == summary[3].unique_graphs);
    // The triangle survives, but the four-cycle has no vertex adjacent to
    // all three other bivalent vertices and is discarded before storage.
    assert((check_file<4, 1>(directory / "first") == 1));
    using Paw = Graph<4, 4, 0, 0, 0, 0, fieldType>;
    Paw paw;
    paw.setEdge(0, 0, 1); paw.setEdge(1, 0, 2);
    paw.setEdge(2, 1, 2); paw.setEdge(3, 0, 3);
    const auto canonical_paw = transient_graph2_standardizer<4, 4>{}.standardize_no_sign(
        transient_graph2<4, 4>(paw));
    const auto stored_paw = MappedTransientGraph2Reader<Paw>(Pipeline::file_path(directory / "first", 4, 4))[0];
    assert(stored_paw == canonical_paw.graph());

    Pipeline again;
    again.run(directory / "second");
    for (const auto& file : std::filesystem::directory_iterator(directory / "first")) {
        assert(file.path().extension() == ".gcg");
        assert(bytes(file.path()) == bytes(directory / "second" / file.path().filename()));
    }

    using Triangle = Graph<3, 3, 0, 0, 0, 0, fieldType>;
    const auto original = Pipeline::file_path(directory / "first", 3, 3);
    const auto corrupt = directory / "corrupt.gcg";
    std::filesystem::copy_file(original, corrupt);
    {
        auto header = read_mapped_graph_file_header(corrupt);
        header.payload_kind = 1; // Reject legacy transient files even with the same dimensions.
        std::fstream file(corrupt, std::ios::in | std::ios::out | std::ios::binary);
        file.write(reinterpret_cast<const char*>(&header), sizeof(header));
    }
    must_throw([&] { MappedTransientGraph2Reader<Triangle> reader(corrupt); });
    std::filesystem::resize_file(corrupt, 10);
    must_throw([&] { MappedTransientGraph2Reader<Triangle> reader(corrupt); });
    must_throw([&] { MappedTransientGraph2Reader<Diamond> reader(original); });

    for (const auto& stage : summary)
        std::cout << "V=" << +stage.vertices << " candidates=" << stage.candidates
                  << " unique=" << stage.unique_graphs << '\n';
    std::filesystem::remove_all(directory);
    std::cout << "graph generation pipeline tests passed\n";
}
