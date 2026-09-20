#include <chrono>
#include <iostream>
#include "GraphHomology/NaturalAdjoint.hpp"
#ifndef GC_HOMOLOGY_LOOP
#define GC_HOMOLOGY_LOOP 9
#endif
#ifndef GC_HOMOLOGY_VERTICES
#define GC_HOMOLOGY_VERTICES 13
#endif

// Independent cross-multiplication, with no division or finite-field arithmetic.
std::size_t verify(const GraphHomology::ContractionMatrix& c,
                   const GraphHomology::NaturalAdjointMatrix& s,
                   std::span<const std::uint64_t> source_aut,
                   std::span<const std::uint64_t> target_aut) {
    if(s.nonzeros()!=c.nonzeros()) throw std::runtime_error("adjoint sparsity mismatch");
    for(std::size_t col=0;col<c.columns();++col) {
        for(auto i=c.offsets[col];i<c.offsets[col+1];++i) {
            const auto row=c.row_indices[i];
            const auto first=s.row_indices.begin()+s.offsets[row];
            const auto last=s.row_indices.begin()+s.offsets[row+1];
            const auto it=std::lower_bound(first,last,col);
            if(it==last || *it!=col) throw std::runtime_error("missing adjoint entry");
            const auto j=it-s.row_indices.begin();
            const __int128 left=__int128(s.coefficients[j])*source_aut[col];
            const __int128 right=__int128(c.coefficients[i])*target_aut[row];
            if(left!=right) throw std::runtime_error("integer adjoint identity failed");
        }
    }
    return c.nonzeros();
}

template<GraphHomology::Parity P>
void run(const char* directory) {
    const auto start=std::chrono::steady_clock::now();
    GraphHomology::ContractionWindow<GC_HOMOLOGY_LOOP,GC_HOMOLOGY_VERTICES,P> w(directory);
    GraphHomology::NaturalAdjoints s(w,2);
    const auto count=verify(w.down,s.down,s.middle_aut,s.lower_aut)
                    +verify(w.up,s.up,s.upper_aut,s.middle_aut);
    std::cout << "L=" << GC_HOMOLOGY_LOOP << " V=" << GC_HOMOLOGY_VERTICES
        << " parity=" << (P==GraphHomology::Parity::even?"even":"odd")
        << " checked_nonzeros=" << count << " all_integral=true exact_identity=true seconds="
        << std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count() << std::endl;
}
int main(int argc,char** argv) {
    if(argc!=2) return 2;
    try {
        run<GraphHomology::Parity::even>(argv[1]);
        run<GraphHomology::Parity::odd>(argv[1]);
    } catch(const std::exception& e) { std::cerr << e.what() << std::endl; return 1; }
}
