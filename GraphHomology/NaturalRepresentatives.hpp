#pragma once
#include "NaturalAdjoint.hpp"

namespace GraphHomology {
template<VectorSpace::Field K>
void check_automorphism_units(const NaturalAdjoints& adjoints) {
    for(const auto* sizes:{&adjoints.lower_aut,&adjoints.middle_aut,&adjoints.upper_aut})
        for(auto a:*sizes)
            if(a % std::uint64_t(K::characteristic()) == 0)
                throw std::runtime_error("automorphism order is not invertible in this field; use another prime");
}

// Convert z to x=D_middle^-1 z, then normalize in original graph coordinates.
// A Laplacian kernel over a finite field can contain noncycles: verify BOTH
// original contraction and natural splitting residuals before saving anything.
template<VectorSpace::Field K>
void convert_natural_representatives(std::vector<std::vector<K>>& basis,
    const ContractionMatrix& down, const NaturalAdjoints& adjoints) {
    check_automorphism_units<K>(adjoints);
    VectorSpace::OwnedArray<K> inverse_aut(adjoints.middle_aut.size());
    for(std::size_t i=0;i<inverse_aut.size();++i)
        inverse_aut[i]=K(adjoints.middle_aut[i] % std::uint64_t(K::characteristic())).inv();
    for(auto& x:basis) {
        if(x.size()!=inverse_aut.size())throw std::invalid_argument("representative dimension mismatch");
        for(std::size_t i=0;i<x.size();++i)x[i]*=inverse_aut[i];
        for(auto value:down.template apply<K>(x))
            if(value!=K{})throw std::runtime_error("natural composition kernel contains a noncycle; use another prime");
        for(auto value:adjoints.up.template apply<K>(x))
            if(value!=K{})throw std::runtime_error("natural composition kernel fails the adjoint residual; use another prime");
        const auto pivot=std::find_if(x.begin(),x.end(),[](const auto& value){return value!=K{};});
        if(pivot==x.end())throw std::runtime_error("zero natural representative");
        const auto scale=pivot->inv();
        for(auto& value:x)value*=scale;
    }
}
}
