#include "GraphGeneration/GraphGenerationPipeline.hpp"
#include <cassert>
#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <sstream>

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
        assert(valences.front() >= 2);
        std::set<std::pair<Int, Int>> edges;
        for (Int e = 0; e < G::N_EDGES_; ++e) {
            const auto [a, b] = graph.getEdge(e);
            assert(a < b && b < V);
            assert(edges.emplace(a, b).second);
            assert(static_cast<unsigned char>(contents[64 + i * G::SIZE + 2 * e]) == a);
            assert(static_cast<unsigned char>(contents[64 + i * G::SIZE + 2 * e + 1]) == b);
        }

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
    using Cycle = Graph<4, 4, 0, 0, 0, 0, fieldType>;
    Cycle cycle;
    cycle.setEdge(0, 0, 1); cycle.setEdge(1, 1, 2);
    cycle.setEdge(2, 2, 3); cycle.setEdge(3, 0, 3);
    assert((bivalent_split_vertices(cycle, cycle.valence_array())
            == std::array<bool, 4>{false, false, false, false}));
    using Path = Graph<3, 2, 0, 0, 0, 0, fieldType>;
    Path path;
    path.setEdge(0, 0, 2); path.setEdge(1, 1, 2);
    assert((bivalent_split_vertices(path, path.valence_array())
            == std::array<bool, 3>{true, true, false}));
    using OrderedDiamond = Graph<4, 5, 0, 0, 0, 0, fieldType>;
    OrderedDiamond ordered;
    ordered.setEdge(0, 0, 2); ordered.setEdge(1, 0, 3);
    ordered.setEdge(2, 1, 2); ordered.setEdge(3, 1, 3); ordered.setEdge(4, 2, 3);
    assert((bivalent_split_vertices(ordered, ordered.valence_array())
            == std::array<bool, 4>{false, false, true, true}));
    using ClawRule = Graph<4, 3, 0, 0, 0, 0, fieldType>;
    ClawRule claw;
    for (Int i = 0; i < 3; ++i) claw.setEdge(i, i, 3);
    assert((bivalent_split_vertices(claw, claw.valence_array())
            == std::array<bool, 4>{true, true, true, true}));
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
    std::ostringstream progress;
    const auto summary = pipeline.run(directory / "first", &progress);
    // Splitting K3 emits K4 and two equivalent diamonds; extended-edge creation is disabled.
    assert(progress.str().find("COUNT\t3\t4\t6\t1\n") != std::string::npos);
    assert(progress.str().find("COUNT\t2\t4\t5\t0\n") != std::string::npos);
    assert(progress.str().find("COUNT\t1\t4\t4\t0\n") != std::string::npos);
    assert(progress.str().find("Starting V=3") < progress.str().find("Finished V=3"));
    assert(summary.size() == 4);
    assert(summary[0].unique_graphs == 1); // K3
    assert(summary[1].unique_graphs == 2); // K4 and diamond.
    assert(summary[1].candidates == 3); // One K4 and two diamond split assignments.
    assert(check_stage<3>(directory / "first") == summary[0].unique_graphs);
    assert(check_stage<4>(directory / "first") == summary[1].unique_graphs);
    assert(check_stage<5>(directory / "first") == summary[2].unique_graphs);
    assert(check_stage<6>(directory / "first") == summary[3].unique_graphs);
    // The four-cycle assignment is rejected before creating a child.
    assert((check_file<4, 1>(directory / "first") == 0));
    assert((check_file<4, 2>(directory / "first") == 1)); // A bivalent split produces the diamond.
    assert((check_file<3, 0>(directory / "first") == 0)); // No star seeds.

    // Even at a larger loop budget, only the triangle seed is inserted.
    GraphGenerationPipeline<6, 4> higher_budget;
    higher_budget.run(directory / "higher_budget");
    assert((check_file<3, 0>(directory / "higher_budget") == 0));
    assert((check_file<3, 1>(directory / "higher_budget") == 1));
    assert((check_file<4, 0>(directory / "higher_budget") == 0));
    assert((check_file<4, 3>(directory / "higher_budget") == 1));

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
