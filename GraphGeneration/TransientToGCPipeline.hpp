#pragma once

#include <cassert>
#include <cstdlib>

#include "GraphGeneration/FinalCanonicalization.hpp"
#include "GraphGeneration/MappedGCGraphFile.hpp"
#include "GraphGeneration/MappedTransientGraph2File.hpp"

namespace GraphGeneration {

template <Int Vertices, Int Edges>
class TransientToGCPipeline {
public:
    using graph_type = Graph<Vertices, Edges, 0, 0, 0, 0, fieldType>;

    GCGraphCounts run(const std::filesystem::path& input,
                      const std::filesystem::path& output) const {
        if (std::filesystem::exists(output))
            std::abort(); // Never overwrite existing output, including release builds.
        MappedTransientGraph2Reader<graph_type> reader(input);
        std::uint64_t retained = 0;
        for (std::uint64_t i = 0; i < reader.size(); ++i)
            retained += is_at_least_trivalent(reader[i]);
        if (output.has_parent_path()) std::filesystem::create_directories(output.parent_path());
        const auto temporary = output.string() + ".tmp";
        GCGraphCounts counts;
        {
            // Count first, then standardize once per retained record directly into mmap.
            MappedGCGraphWriter<graph_type> writer(temporary, retained);
            final_graph_canonicalizer<graph_type, true> canonicalize;
            for (std::uint64_t i = 0; i < reader.size(); ++i) {
                const auto graph = reader[i];
                if (!is_at_least_trivalent(graph)) continue;
                auto result = canonicalize(graph);
                writer.append({std::move(result.canonical_graph),
                    result.survives_even_edges_odd_vertices(), result.survives_odd_edges()});
            }
            counts = writer.finish();
        }
        std::filesystem::rename(temporary, output);
        return counts;
    }
private:
    static bool is_at_least_trivalent(const graph_type& graph) {
        for (Int vertex : graph.half_edges)
            assert(vertex < Vertices);
        const auto valences = graph.valence_array();
        return std::ranges::all_of(valences, [](Int valence) { return valence >= 3; });
    }
};

} // namespace GraphGeneration
