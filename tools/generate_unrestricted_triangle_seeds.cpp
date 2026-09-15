#include "GraphGeneration/UnrestrictedTriangleSeedGenerator.hpp"

#ifndef GC_TRIANGLE_MAX_LOOP
#define GC_TRIANGLE_MAX_LOOP 6
#endif

int main(int argc, char** argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " NEW_OUTPUT_DIRECTORY\n";
        return 2;
    }
    try {
        std::cout << "Standardizer: full-graph key; no component stack; cubic completion budget\n";
        std::cout << "Triangle seed proof of concept: K3, max_loop=" << GC_TRIANGLE_MAX_LOOP
                  << " max_vertices=" << (2 * GC_TRIANGLE_MAX_LOOP + 1) << std::endl;
        GraphGeneration::UnrestrictedTriangleSeedGenerator<GC_TRIANGLE_MAX_LOOP>{}.run(argv[1]);
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
