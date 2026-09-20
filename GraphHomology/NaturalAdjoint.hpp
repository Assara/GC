#pragma once
#include <numeric>
#include "ContractionMatrices.hpp"

namespace GraphHomology {
// Natural splitting coefficients need not obey the contraction edge-count bound.
using NaturalAdjointMatrix = ContractionMatrixStorage<VectorSpace::OwnedArray, std::int32_t>;
using AutomorphismSizes = VectorSpace::OwnedArray<std::uint64_t>;

template<int L, int V, Parity P>
AutomorphismSizes automorphism_sizes(const EnumeratedBasis<L,V,P>& basis, int threads = 8) {
    if (threads < 1) throw std::invalid_argument("thread count must be positive");
    AutomorphismSizes sizes(basis.size());
#pragma omp parallel for schedule(dynamic, 16) num_threads(threads) if(basis.size() >= 128)
    for (std::size_t i = 0; i < basis.size(); ++i) {
        GraphStandardizer<V,L+V-1,0,0,P == Parity::odd ? 1 : 0,1,fieldType> standardizer;
        using G = typename EnumeratedBasis<L,V,P>::G;
        sizes[i] = standardizer.standardize4_with_aut_count(typename G::Basis(basis.graphs[i])).second;
    }
    return sizes;
}

// C: source -> target. S[source,target] = C[target,source] * aut(target)/aut(source).
// Cancel the ratio first, so no product of two large automorphism sizes is needed.
inline NaturalAdjointMatrix natural_adjoint(const ContractionMatrix& contraction,
    std::span<const std::uint64_t> source_aut, std::span<const std::uint64_t> target_aut) {
    if (source_aut.size() != contraction.columns() || target_aut.size() != contraction.rows)
        throw std::invalid_argument("automorphism dimensions do not match contraction");
    for (auto a : source_aut) if (!a) throw std::invalid_argument("zero automorphism size");
    for (auto a : target_aut) if (!a) throw std::invalid_argument("zero automorphism size");
    NaturalAdjointMatrix result;
    result.rows = contraction.columns();
    result.offsets = VectorSpace::OwnedArray<std::size_t>(std::size_t(contraction.rows)+1);
    result.row_indices = VectorSpace::OwnedArray<std::uint32_t>(contraction.nonzeros());
    result.coefficients = VectorSpace::OwnedArray<std::int32_t>(contraction.nonzeros());
    for (auto row : contraction.row_indices) ++result.offsets[std::size_t(row)+1];
    std::partial_sum(result.offsets.begin(), result.offsets.end(), result.offsets.begin());
    VectorSpace::OwnedArray<std::size_t> next{std::span<const std::size_t>(result.offsets)};
    for (std::size_t c = 0; c < contraction.columns(); ++c) {
        for (auto i = contraction.offsets[c]; i < contraction.offsets[c+1]; ++i) {
            const auto r = contraction.row_indices[i];
            const auto common = std::gcd(source_aut[c], target_aut[r]);
            const auto denominator = source_aut[c]/common;
            const auto numerator = target_aut[r]/common;
            const int coefficient = contraction.coefficients[i];
            const auto magnitude = std::uint64_t(coefficient < 0 ? -coefficient : coefficient);
            if (magnitude % denominator)
                throw std::runtime_error("natural adjoint coefficient is not integral: c=" + std::to_string(c)
                    + " r=" + std::to_string(r) + " coefficient=" + std::to_string(coefficient)
                    + " source_aut=" + std::to_string(source_aut[c]) + " target_aut=" + std::to_string(target_aut[r]));
            const auto quotient = magnitude / denominator;
            if (quotient && numerator > std::uint64_t(INT32_MAX)/quotient)
                throw std::overflow_error("natural adjoint coefficient exceeds int32");
            const auto value = std::int32_t(quotient * numerator);
            const auto j = next[r]++;
            result.row_indices[j] = c;
            result.coefficients[j] = coefficient < 0 ? -value : value;
        }
    }
    return result;
}

// Calculate each basis's automorphism sizes once, sharing the middle sizes.
struct NaturalAdjoints {
    AutomorphismSizes lower_aut, middle_aut, upper_aut;
    NaturalAdjointMatrix down, up;
    NaturalAdjoints() = default;
    template<int L, int V, Parity P>
    explicit NaturalAdjoints(const ContractionWindow<L,V,P>& window, int threads = 8)
        : lower_aut(automorphism_sizes(window.lower, threads)),
          middle_aut(automorphism_sizes(window.middle, threads)),
          upper_aut(automorphism_sizes(window.upper, threads)),
          down(natural_adjoint(window.down, middle_aut, lower_aut)),
          up(natural_adjoint(window.up, upper_aut, middle_aut)) {}
    std::size_t allocated_bytes() const {
        return (lower_aut.size()+middle_aut.size()+upper_aut.size())*sizeof(std::uint64_t)
            + down.allocated_bytes()+up.allocated_bytes();
    }
};
}
