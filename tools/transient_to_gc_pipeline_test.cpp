#include "GraphGeneration/TransientToGCPipeline.hpp"
#include <cassert>
#include <fstream>
#include <iostream>
#include <vector>

using namespace GraphGeneration;

template <typename F>
void must_throw(F&& action) {
    bool threw = false;
    try { action(); } catch (const std::exception&) { threw = true; }
    assert(threw);
}

template <typename G>
G make_graph(std::initializer_list<std::pair<Int, Int>> edges) {
    assert(edges.size() == G::N_EDGES_);
    G graph;
    Int e = 0;
    for (auto [a, b] : edges) graph.setEdge(e++, a, b);
    return graph;
}

template <signedInt C, signedInt D, typename G>
bool signed_survival(const G& graph) {
    using Signed = typename G::template RebindDegree<C, D>;
    Signed signed_graph;
    signed_graph.half_edges = graph.half_edges;
    return GraphStandardizer<G::N_VERTICES_, G::N_EDGES_, 0, 0, C, D,
        typename G::Field>{}.standardize4(typename Signed::Basis(signed_graph)).getCoefficient()
        != typename G::Field{};
}

template <typename G>
void write_input(const std::filesystem::path& path, const std::vector<G>& graphs) {
    MappedTransientGraph2Writer<G> writer(path, graphs.size());
    for (std::size_t i = 0; i < graphs.size(); ++i) writer.write(i, graphs[i]);
}

int main() {
    char name[] = "/tmp/gc-convert-test-XXXXXX";
    const char* created = ::mkdtemp(name);
    assert(created);
    const std::filesystem::path directory(created);

    // K4: an odd vertex permutation also reverses an odd number of edges.
    // Their signs cancel, so the tetrahedron survives both conventions.
    using Tetrahedron = Graph<4, 6, 0, 0, 0, 0, fieldType>;
    const auto tetrahedron = make_graph<Tetrahedron>({{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}});
    write_input<Tetrahedron>(directory / "tetrahedron.gcg", {tetrahedron});
    const auto tetra_counts = TransientToGCPipeline<4, 6>{}.run(
        directory / "tetrahedron.gcg", directory / "tetrahedron_out.gcg");
    assert(tetra_counts.total == 1 && tetra_counts.odd_gc == 1 && tetra_counts.even_gc == 1);
    assert((signed_survival<1, 1>(tetrahedron) && signed_survival<0, 1>(tetrahedron)));
    std::array<Int, 4> permutation{0, 1, 2, 3};
    do {
        auto relabeled = tetrahedron;
        for (Int e = 0; e < 6; ++e) {
            const auto [a, b] = tetrahedron.getEdge(e);
            relabeled.setEdge(e, permutation[b], permutation[a]);
        }
        const auto result = final_graph_canonicalizer<Tetrahedron, true>{}(relabeled);
        assert(result.canonical_graph == tetrahedron);
        assert(result.survives_odd_edges() && result.survives_even_edges_odd_vertices());
    } while (std::next_permutation(permutation.begin(), permutation.end()));

    using G = Graph<6, 10, 0, 0, 0, 0, fieldType>;
    const auto wheel = make_graph<G>({{0,1},{0,2},{0,3},{0,4},{0,5},
                                     {1,2},{1,3},{2,4},{3,5},{4,5}});
    const auto even = make_graph<G>({{0,1},{0,2},{0,3},{0,4},{1,2},
                                     {1,3},{1,5},{2,4},{3,5},{4,5}});
    auto bivalent = wheel;
    bivalent.setEdge(4, 1, 4); // Vertex 5 now has valence two.
    // A leaf with no bivalent vertices must also be excluded from GC files.
    const auto monovalent = make_graph<G>({{0,1},{0,2},{0,3},{0,4},{0,5},
                                          {1,2},{1,3},{1,4},{2,3},{2,4}});
    const auto input = directory / "input.gcg";
    const auto output = directory / "output.gcg";
    write_input<G>(input, {wheel, bivalent, monovalent, even});
    const auto counts = TransientToGCPipeline<6, 10>{}.run(input, output);
    assert(counts.total == 2 && counts.odd_gc == 0 && counts.even_gc == 2);
    MappedGCGraphReader<G> reader(output);
    assert(reader.counts().total == 2 && reader.counts().odd_gc == 0 && reader.counts().even_gc == 2);
    assert(std::filesystem::file_size(output) == 64 + 2 * (20 + 1));
    for (std::size_t i = 0; i < reader.size(); ++i) {
        const auto record = reader[i];
        assert((record.odd_gc == signed_survival<1, 1>(record.graph)));
        assert((record.even_gc == signed_survival<0, 1>(record.graph)));
        assert(record.graph == canonicalize_final_graph(record.graph).canonical_graph);
    }
    assert(!reader[0].odd_gc && reader[0].even_gc);
    assert(!reader[1].odd_gc && reader[1].even_gc);
    // Verify little-endian count fields and byte flags in the mmap output.
    {
        std::ifstream file(output, std::ios::binary);
        std::vector<unsigned char> bytes{std::istreambuf_iterator<char>(file), {}};
        assert(bytes[8] == 2 && bytes[32] == 2 && bytes[40] == 0 && bytes[48] == 2);
        assert(bytes[64 + 20] == 2 && bytes[64 + 21 + 20] == 2);
    }
    must_throw([&] { (void)reader[2]; });

    // A graph that vanishes in both complexes is retained in total.
    using K33 = Graph<6, 9, 0, 0, 0, 0, fieldType>;
    const auto k33 = make_graph<K33>({{0,3},{0,4},{0,5},{1,3},{1,4},
                                       {1,5},{2,3},{2,4},{2,5}});
    write_input<K33>(directory / "k33.gcg", {k33});
    const auto zero = TransientToGCPipeline<6, 9>{}.run(directory / "k33.gcg", directory / "zero.gcg");
    assert(zero.total == 1 && zero.odd_gc == 0 && zero.even_gc == 0);
    MappedGCGraphReader<K33> zero_reader(directory / "zero.gcg");
    assert(!zero_reader[0].odd_gc && !zero_reader[0].even_gc);
    assert((!signed_survival<1, 1>(k33) && !signed_survival<0, 1>(k33)));

    // The two flags are independent: this graph survives only oddGC.
    using Odd = Graph<7, 12, 0, 0, 0, 0, fieldType>;
    const auto odd = make_graph<Odd>({{0,1},{0,2},{0,3},{0,4},{0,5},{1,2},
                                     {1,3},{1,6},{2,4},{3,5},{4,6},{5,6}});
    write_input<Odd>(directory / "odd.gcg", {odd});
    const auto odd_counts = TransientToGCPipeline<7, 12>{}.run(directory / "odd.gcg", directory / "odd_out.gcg");
    assert(odd_counts.total == 1 && odd_counts.odd_gc == 1 && odd_counts.even_gc == 0);
    const auto odd_record = MappedGCGraphReader<Odd>(directory / "odd_out.gcg")[0];
    assert(odd_record.odd_gc && !odd_record.even_gc);
    assert((signed_survival<1, 1>(odd) && !signed_survival<0, 1>(odd)));

    // All-filtered and already-empty input files still produce valid headers.
    write_input<G>(directory / "filtered.gcg", {bivalent, monovalent});
    const auto filtered = TransientToGCPipeline<6, 10>{}.run(directory / "filtered.gcg", directory / "empty.gcg");
    assert(filtered.total == 0 && filtered.odd_gc == 0 && filtered.even_gc == 0);
    assert(std::filesystem::file_size(directory / "empty.gcg") == 64);
    assert(MappedGCGraphReader<G>(directory / "empty.gcg").size() == 0);
    write_input<G>(directory / "none.gcg", {});
    assert((TransientToGCPipeline<6, 10>{}.run(directory / "none.gcg", directory / "none_out.gcg").total == 0));
    must_throw([&] { MappedGCGraphReader<G> wrong(input); });
    must_throw([&] { MappedGCGraphReader<K33> wrong(output); });
    const auto corrupt = directory / "corrupt.gcg";
    std::filesystem::copy_file(output, corrupt);
    {
        std::fstream file(corrupt, std::ios::in | std::ios::out | std::ios::binary);
        file.seekp(8); file.put(char(1)); // Reject old vertex-only sign files.
    }
    must_throw([&] { MappedGCGraphReader<G> wrong(corrupt); });
    {
        std::fstream file(corrupt, std::ios::in | std::ios::out | std::ios::binary);
        file.seekp(8); file.put(char(2));
        file.seekp(64 + 20); file.put(char(4));
    }
    must_throw([&] { (void)MappedGCGraphReader<G>(corrupt)[0]; });
    std::filesystem::resize_file(corrupt, 12);
    must_throw([&] { MappedGCGraphReader<G> wrong(corrupt); });
    std::filesystem::remove_all(directory);
    std::cout << "transient to GC pipeline tests passed\n";
}
