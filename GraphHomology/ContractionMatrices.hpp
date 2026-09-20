#pragma once

#include <algorithm>
#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <span>
#include <stdexcept>
#include <string>
#include <vector>
#include "graph.hpp"
#include "VectorSpace/Field.hpp"
#include "VectorSpace/OwnedArray.hpp"
#include "VectorSpace/mmap_archive.hpp"
#include "GraphGeneration/CutVertexSplitRule.hpp"

namespace GraphHomology {

// Names refer to GC parity, not the historical OddGC alias (odd edges).
enum class Parity { even, odd };
template <int L, int V, Parity P>
using BasisGraph = Graph<V, L + V - 1, 0, 0, P == Parity::odd ? 1 : 0, 1, fieldType>;

// Separate arrays avoid padding an (uint32_t, int8_t) pair to eight bytes.
template<template<class...> class Storage, class Coefficient = SmallSignedInt>
struct ContractionMatrixStorage {
    std::uint32_t rows = 0;
    Storage<std::size_t> offsets = Storage<std::size_t>(1);
    Storage<std::uint32_t> row_indices;
    Storage<Coefficient> coefficients;

    std::size_t columns() const { return offsets.size() - 1; }
    std::size_t nonzeros() const { return coefficients.size(); }
    std::size_t allocated_bytes() const {
        return offsets.size() * sizeof(std::size_t)
            + row_indices.size() * sizeof(std::uint32_t)
            + coefficients.size() * sizeof(Coefficient);
    }
    void save(VectorSpace::serialization::archive_output& out) const {
        out.section("sparse-matrix", [&] {
            out.word(1); // CSC layout
            out.word(rows); out.word(columns()); out.word(nonzeros());
            out.integers<std::size_t>("offsets", offsets);
            out.integers<std::uint32_t>("row-indices", row_indices);
            out.integers<Coefficient>("coefficients", coefficients);
        });
    }
    static ContractionMatrixStorage load(VectorSpace::serialization::archive_input& in) {
        ContractionMatrixStorage matrix;
        in.section("sparse-matrix", [&] {
            in.expect(1);
            const auto rows=in.word(), columns=in.word(), entries=in.word();
            if(rows>UINT32_MAX || columns>=std::numeric_limits<std::size_t>::max()
                || entries>std::numeric_limits<std::size_t>::max()/sizeof(std::size_t))
                throw std::runtime_error("invalid sparse matrix dimensions");
            matrix.rows=rows;
            matrix.offsets=Storage<std::size_t>(columns+1);
            matrix.row_indices=Storage<std::uint32_t>(entries);
            matrix.coefficients=Storage<Coefficient>(entries);
            in.integers<std::size_t>("offsets", matrix.offsets);
            in.integers<std::uint32_t>("row-indices", matrix.row_indices);
            in.integers<Coefficient>("coefficients", matrix.coefficients);
            if(matrix.offsets[0]!=0 || matrix.offsets.back()!=entries
                || !std::is_sorted(matrix.offsets.begin(),matrix.offsets.end())
                || std::any_of(matrix.row_indices.begin(),matrix.row_indices.end(),
                    [&](auto row){return row>=rows;}))
                throw std::runtime_error("invalid sparse matrix structure");
        });
        return matrix;
    }
    template <VectorSpace::Field Field>
    std::vector<Field> apply(std::span<const Field> input) const {
        if (input.size() != columns()) throw std::invalid_argument("matrix input dimension mismatch");
        std::vector<Field> result(rows);
        for (std::size_t c = 0; c < columns(); ++c)
            for (auto i = offsets[c]; i < offsets[c + 1]; ++i)
                result[row_indices[i]] += input[c] * coefficients[i];
        return result;
    }
};

using ContractionMatrix = ContractionMatrixStorage<VectorSpace::OwnedArray>;
// Growth is confined to construction; the solver receives fixed owning arrays.
struct ContractionMatrixBuilder : ContractionMatrixStorage<std::vector> {
    // Input may contain repeated rows. Accumulate in int, cancel, then narrow.
    void append_column(std::vector<std::pair<std::uint32_t, int>>& entries) {
        std::sort(entries.begin(), entries.end());
        for (std::size_t i = 0; i < entries.size();) {
            auto row = entries[i].first;
            if (row >= rows) throw std::out_of_range("matrix row out of range");
            std::int64_t sum = 0;
            do { sum += entries[i++].second; } while (i < entries.size() && entries[i].first == row);
            if (!sum) continue;
            row_indices.push_back(row);
            coefficients.push_back(static_cast<SmallSignedInt>(sum));
        }
        offsets.push_back(nonzeros());
    }
    ContractionMatrix finish() && {
        ContractionMatrix result;
        result.rows = rows;
        result.offsets = VectorSpace::OwnedArray<std::size_t>(std::span<const std::size_t>(offsets));
        std::vector<std::size_t>().swap(offsets);
        result.row_indices = VectorSpace::OwnedArray<std::uint32_t>(std::span<const std::uint32_t>(row_indices));
        std::vector<std::uint32_t>().swap(row_indices);
        result.coefficients = VectorSpace::OwnedArray<SmallSignedInt>(std::span<const SmallSignedInt>(coefficients));
        std::vector<SmallSignedInt>().swap(coefficients);
        return result;
    }
};

template <int L, int V, Parity P>
struct EnumeratedBasis {
    static_assert(L >= 3 && V >= 3 && V <= 62);
    using G = BasisGraph<L, V, P>;
    std::vector<G> graphs;
    std::size_t input_graphs = 0, zero_graphs = 0;

    std::uint32_t size() const { return static_cast<std::uint32_t>(graphs.size()); }
    std::size_t allocated_bytes() const { return graphs.capacity() * sizeof(G); }
    std::uint32_t index(const G& graph) const {
        auto it = std::lower_bound(graphs.begin(), graphs.end(), graph,
            [](const G& a, const G& b) { return a.half_edges < b.half_edges; });
        if (it == graphs.end() || it->half_edges != graph.half_edges)
            throw std::runtime_error("nonzero contraction missing from target basis at L="
                + std::to_string(L) + " V=" + std::to_string(V)
                + (G(graph).has_double_edge() ? " (parallel-edge graph; simple graph6 basis is not closed)" : ""));
        return static_cast<std::uint32_t>(it - graphs.begin());
    }

    static EnumeratedBasis load(const std::filesystem::path& directory) {
        EnumeratedBasis result;
        // Adjacent chain groups beyond the simple minimum-degree-three range are zero.
        if constexpr (2 * G::N_EDGES_ < 3 * V || G::N_EDGES_ > V * (V - 1) / 2) return result;
        auto path = directory / ("graphs_L" + std::to_string(L) + "_V" + std::to_string(V) + ".g6");
        std::ifstream input(path);
        if (!input) throw std::runtime_error("cannot open basis input: " + path.string());
        for (std::string line; std::getline(input, line);) {
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.size() != 1 + (V * (V - 1) / 2 + 5) / 6 || line[0] != V + 63
                || std::ranges::any_of(line, [](unsigned char c) { return c < 63 || c > 126; }))
                throw std::runtime_error("invalid graph6 record in " + path.string());
            G graph;
            int bit = 0, edge = 0;
            for (int b = 1; b < V; ++b)
                for (int a = 0; a < b; ++a, ++bit)
                    if (((line[1 + bit / 6] - 63) >> (5 - bit % 6)) & 1) {
                        if (edge == G::N_EDGES_) throw std::runtime_error("wrong graph loop number");
                        graph.setEdge(edge++, a, b);
                    }
            if (edge != G::N_EDGES_) throw std::runtime_error("wrong graph loop number");
            if (std::ranges::any_of(graph.valence_array(), [](Int d) { return d < 3; })
                || GraphGeneration::CutVertexSplitRule<V>(graph).vertex >= 0)
                throw std::runtime_error("expected minimum-degree-three cut-free input graphs");
            typename G::Basis canonical(graph);
            G::std(canonical);
            ++result.input_graphs;
            if (canonical.getCoefficient() == fieldType{0}) { ++result.zero_graphs; continue; }
            result.graphs.push_back(canonical.getValue());
        }
        if (!input.eof()) throw std::runtime_error("failed reading basis input");
        std::sort(result.graphs.begin(), result.graphs.end(),
            [](const G& a, const G& b) { return a.half_edges < b.half_edges; });
        for (std::size_t i = 1; i < result.graphs.size(); ++i)
            if (result.graphs[i - 1].half_edges == result.graphs[i].half_edges)
                throw std::runtime_error("duplicate graph in basis input");
        if (result.graphs.size() > std::numeric_limits<std::uint32_t>::max())
            throw std::overflow_error("basis requires more than 32-bit indices");
        result.graphs.shrink_to_fit();
        return result;
    }
};

template <int L, int V, Parity P>
ContractionMatrix contraction_matrix(const EnumeratedBasis<L, V, P>& source,
                                    const EnumeratedBasis<L, V - 1, P>& target) {
    using G = BasisGraph<L, V, P>;
    ContractionMatrixBuilder matrix;
    matrix.rows = target.size();
    matrix.offsets.reserve(std::size_t(source.size()) + 1);
    std::vector<std::pair<std::uint32_t, int>> entries;
    entries.reserve(G::N_EDGES_);
    std::vector<std::pair<typename G::ContGraph, int>> contractions;
    contractions.reserve(G::N_EDGES_);
    for (const auto& graph : source.graphs) {
        entries.clear();
        contractions.clear();
        // In a simple parent, contracting an edge creates parallel edges
        // exactly when its endpoints share a neighbour. These targets are
        // zero in the simple-graph quotient, for either parity.
        std::array<std::uint64_t, V> adjacent{};
        for (Int e = 0; e < G::N_EDGES_; ++e) {
            const auto [a, b] = graph.getEdge(e);
            adjacent[a] |= std::uint64_t{1} << b;
            adjacent[b] |= std::uint64_t{1} << a;
        }
        for (Int e = 0; e < G::N_EDGES_; ++e) {
            const auto [a, b] = graph.getEdge(e);
            if (adjacent[a] & adjacent[b]) continue;
            auto term = graph.contract_edge(e, fieldType{1});
            if (term.getCoefficient() == fieldType{0}) continue;
            G::ContGraph::std(term);
            auto coefficient = term.getCoefficient();
            if (coefficient == fieldType{0}) continue;
            int sign;
            if (coefficient == fieldType{1}) sign = 1;
            else if (coefficient == fieldType{-1}) sign = -1;
            else throw std::runtime_error("individual contraction is not a signed unit");
            contractions.emplace_back(term.getValue(), sign);
        }
        std::sort(contractions.begin(), contractions.end(), [](const auto& a, const auto& b) {
            return a.first.half_edges < b.first.half_edges;
        });
        for (std::size_t i = 0; i < contractions.size();) {
            const auto& target_graph = contractions[i].first;
            int coefficient = 0;
            do { coefficient += contractions[i++].second; }
            while (i < contractions.size() && contractions[i].first.half_edges == target_graph.half_edges);
            if (coefficient) entries.emplace_back(target.index(target_graph), coefficient);
        }
        matrix.append_column(entries);
    }
    return std::move(matrix).finish();
}

// Verify every column of d_V d_(V+1), without allocating a product matrix.
inline void check_chain(const ContractionMatrix& down, const ContractionMatrix& up) {
    if (down.columns() != up.rows) throw std::invalid_argument("incompatible consecutive differentials");
    std::vector<std::pair<std::uint32_t, int>> terms;
    for (std::size_t c = 0; c < up.columns(); ++c) {
        terms.clear();
        for (auto i = up.offsets[c]; i < up.offsets[c + 1]; ++i) {
            auto mid = up.row_indices[i];
            for (auto j = down.offsets[mid]; j < down.offsets[mid + 1]; ++j)
                terms.emplace_back(down.row_indices[j], int(up.coefficients[i]) * int(down.coefficients[j]));
        }
        std::sort(terms.begin(), terms.end());
        for (std::size_t i = 0; i < terms.size();) {
            const auto row = terms[i].first;
            std::int64_t sum = 0;
            do { sum += terms[i++].second; } while (i < terms.size() && terms[i].first == row);
            if (sum) throw std::runtime_error("contraction differential does not square to zero");
        }
    }
}

template <int L, int V, Parity P>
struct ContractionWindow {
    EnumeratedBasis<L, V - 1, P> lower;
    EnumeratedBasis<L, V, P> middle;
    EnumeratedBasis<L, V + 1, P> upper;
    ContractionMatrix down, up;
    explicit ContractionWindow(const std::filesystem::path& directory)
        : lower(decltype(lower)::load(directory)), middle(decltype(middle)::load(directory)),
          upper(decltype(upper)::load(directory)),
          down(contraction_matrix(middle, lower)), up(contraction_matrix(upper, middle)) {}
    std::size_t allocated_bytes() const {
        return lower.allocated_bytes() + middle.allocated_bytes() + upper.allocated_bytes()
            + down.allocated_bytes() + up.allocated_bytes();
    }
};
} // namespace GraphHomology
