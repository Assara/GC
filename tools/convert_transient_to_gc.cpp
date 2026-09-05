#include "GraphGeneration/TransientToGCPipeline.hpp"
#include <iostream>

#ifndef GC_CONVERT_VERTICES
#define GC_CONVERT_VERTICES 4
#endif
#ifndef GC_CONVERT_EDGES
#define GC_CONVERT_EDGES 6
#endif

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " INPUT_TRANSIENT_FILE OUTPUT_GC_FILE\n";
        return 2;
    }
    const auto counts = GraphGeneration::TransientToGCPipeline<GC_CONVERT_VERTICES,
        GC_CONVERT_EDGES>{}.run(argv[1], argv[2]);
    std::cout << "total=" << counts.total << " oddGC=" << counts.odd_gc
              << " evenGC=" << counts.even_gc << '\n';
}
