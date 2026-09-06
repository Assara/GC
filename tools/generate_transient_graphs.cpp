#include "GraphGeneration/GraphGenerationPipeline.hpp"
#include <iostream>

#ifndef GC_PIPELINE_MAX_LOOP
#define GC_PIPELINE_MAX_LOOP 3
#endif
#ifndef GC_PIPELINE_MAX_VERTICES
#define GC_PIPELINE_MAX_VERTICES 6
#endif

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " OUTPUT_DIRECTORY\n";
        return 2;
    }
    GraphGeneration::GraphGenerationPipeline<GC_PIPELINE_MAX_LOOP,
        GC_PIPELINE_MAX_VERTICES> pipeline;
    std::cout << "TransientGraph2 (bivalent common-neighbour selection, root-preserving splits, triangle seed only, minimum split valence 2"
              << ", preserve second-largest valence): max_loop=" << GC_PIPELINE_MAX_LOOP
              << " max_vertices=" << GC_PIPELINE_MAX_VERTICES << std::endl;
    std::cout << "COUNT_SCOPE biconnected_min_degree_3" << std::endl;
    pipeline.run(argv[1], &std::cout);
}
