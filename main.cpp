#include <filesystem>
#include <fstream>
#include <iostream>
#include "graph.hpp"
#include "GC.hpp"
#include "examplegraphs.hpp"

#define GC_CONTRACTION_PLAYGROUND_NO_MAIN
#include "GC_contraction_playground.cpp"

int main() {
    using WheelGC = OddGCdegZero<20>;
    constexpr const char* kOutputPath = "/home/assar/Projects/GC/output/classes/W19_class.txt";

    WheelGC start(wheel_graph<19>());
    auto rep = try_find_quadratic_contraction_representative_via_min_triangle_split_solver(start);
    if (!rep.has_value()) {
        std::cout << "no solution" << std::endl;
        return 1;
    }

    auto gamma = *rep;
    gamma.standardize_all();
    gamma.sort_elements();

    std::filesystem::create_directories("/home/assar/Projects/GC/output/classes");
    std::ofstream out(kOutputPath, std::ios::trunc);
    if (!out) {
        std::cerr << "failed to open output file: " << kOutputPath << std::endl;
        return 1;
    }

    out << "graph_size: (20,38)\n";
    out << "field_type: " << fieldType::name() << "\n";
    out << "number_of_graphs: " << gamma.size() << "\n";
    for (const auto& be : gamma.data()) {
        out << be.getCoefficient() << "; ";
        for (Int e = 0; e < WheelGC::GraphType::N_EDGES_; ++e) {
            auto edge = be.getValue().getEdge(e);
            out << "(" << +edge.first << "," << +edge.second << ")";
            if (e + 1 < WheelGC::GraphType::N_EDGES_) {
                out << ", ";
            }
        }
        out << "\n";
    }
    out.close();

    std::cout << "saved class to " << kOutputPath << std::endl;
    std::cout << "representative size = " << gamma.size() << std::endl;

    auto d_final = gamma.d_contraction();
    std::cout << "d_contraction size = " << d_final.size() << std::endl;
    return 0;
}
