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
    for (const auto& stage : pipeline.run(argv[1])) {
        std::cout << "V=" << +stage.vertices << " candidates=" << stage.candidates
                  << " unique=" << stage.unique_graphs << '\n';
    }
}
