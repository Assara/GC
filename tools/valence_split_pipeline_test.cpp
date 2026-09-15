#include "GraphGeneration/TriangleSplitCommon.hpp"
#include "GraphGeneration/ValenceSplitPlan.hpp"

using namespace GraphGeneration;
using namespace GraphGeneration::triangle_split_detail;
std::size_t checked_graphs = 0, checked_splits = 0;

template <int V>
SplitStageGraph<V> parse(const std::string& line) {
    if (line.size() != 1 + (V * (V - 1) / 2 + 5) / 6 || line[0] != V + 63)
        throw std::runtime_error("bad graph6 dimensions");
    SplitStageGraph<V> graph;
    int bit = 0, edges = 0;
    for (int b = 1; b < V; ++b)
        for (int a = 0; a < b; ++a, ++bit)
            if (((line[1 + bit / 6] - 63) >> (5 - bit % 6)) & 1) {
                assert(edges < graph.N_EDGES_);
                graph.setEdge(edges++, a, b);
            }
    assert(edges == graph.N_EDGES_);
    return graph;
}

template <int V>
void compare(const std::filesystem::path& reference, const std::filesystem::path& result) {
    GraphStandardizer<V, V - 1 + loop_number, 0, 0, 0, 0, fieldType> standardizer;
    const auto read = [&](const auto& path, bool grouped) {
        std::ifstream input(path / ("graphs_L" + std::to_string(loop_number) + "_V" + std::to_string(V) + ".g6"));
        if (!input) throw std::runtime_error("missing test input");
        std::vector<SplitStageGraph<V>> graphs;
        for (std::string line; std::getline(input, line);) {
            auto graph = parse<V>(line);
            if (grouped) {
                const auto degrees = graph.valence_array();
                assert(std::is_sorted(degrees.begin(), degrees.end()));
            }
            graphs.push_back(standardizer.standardize_no_sign(graph));
        }
        assert(input.eof());
        std::sort(graphs.begin(), graphs.end(), [](auto& a, auto& b) { return a.half_edges < b.half_edges; });
        for (std::size_t i = 1; i < graphs.size(); ++i)
            assert(graphs[i - 1].half_edges != graphs[i].half_edges);
        return graphs;
    };
    const auto expected = read(reference, false), actual = read(result, true);
    assert(expected.size() == actual.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        assert(expected[i].half_edges == actual[i].half_edges);
        ++checked_graphs;
        if (i >= 100) continue;
        const auto graph = standardize_sorted_valences(expected[i], standardizer);
        assert(standardize_sorted_valences(graph, standardizer).half_edges == graph.half_edges);
        Permutation<V> reverse;
        for (int v = 0; v < V; ++v) reverse[v] = V - 1 - v;
        SplitStageGraph<V> relabelled;
        relabelled.assignPermutedDirectedSortedEdgesNoSign(graph, reverse);
        assert(standardize_sorted_valences(relabelled, standardizer).half_edges == graph.half_edges);
        if constexpr (V < max_vertices) {
            const auto degrees = graph.valence_array();
            const auto plans = valence_split_plans(ValenceArray(degrees.begin(), degrees.end()));
            std::size_t plan_pairs = 0, expected_pairs = 0;
            for (const auto& [key, splits] : plans) plan_pairs += splits.size();
            for (int vertex = 0; vertex < V; ++vertex) {
                const auto adjacent = graph.adjacent(vertex);
                for (Int moved = 2; moved + 1 < adjacent.size(); ++moved) {
                    ++expected_pairs;
                    auto subset = combutils::firstSubset(1, moved);
                    do {
                        auto child_degrees = graph.splitGraph(vertex, adjacent, subset).valence_array();
                        std::sort(child_degrees.begin(), child_degrees.end());
                        const auto& splits = plans.at(ValenceArray(child_degrees.begin(), child_degrees.end()));
                        assert(std::count(splits.begin(), splits.end(), ValenceSplit{vertex, moved}) == 1);
                        ++checked_splits;
                    } while (combutils::nextSubset(subset, Int(adjacent.size() - 1)));
                }
            }
            assert(plan_pairs == expected_pairs);
        }
    }
    std::cout << "V=" << V << " exact graphs=" << actual.size() << '\n';
    if constexpr (V < max_vertices) compare<V + 1>(reference, result);
}
int main(int argc, char** argv) {
    if (argc != 3) return 2;
    compare<first_vertices>(argv[1], argv[2]);
    std::cout << "Verified " << checked_graphs << " graphs and " << checked_splits << " split routings\n";
}
