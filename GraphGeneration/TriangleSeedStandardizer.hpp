#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <numeric>
#include <span>
#include <stdexcept>
#include <vector>
#include "GraphGeneration/TriangleSeedEntry.hpp"

namespace GraphGeneration {

// Unsigned individualization/refinement, following standardize4's colour
// updates, partition splitting, best-attempt selection and final edge comparison.
// Component incidence supplies the initial vertex separations and refinement.
// No Graph, TransientGraph2 or general-purpose standardizer dependency.
template <int MaxLoop>
class TriangleSeedStandardizer {
public:
    enum class Orientation { Keep, Discard, Both };

private:
    static constexpr int MaxVertices = 2 * MaxLoop + 1;
    static_assert(MaxLoop >= 1 && MaxVertices <= 62);
    using Mask = std::uint64_t;
    using Entry = TriangleSeedEntry<MaxLoop>;
    using Components = TriangleComponentStack<MaxLoop>;
    struct Attempt {
        std::array<Mask, MaxVertices> colours{};
        std::array<Mask, MaxLoop> component_colours{};
        std::array<int, MaxVertices> order{}, ends{};
        int groups = 0;
        bool reversed = false;
    };
    static Mask hash(Mask n) noexcept {
        n += 0x9e3779b97f4a7c15ULL;
        n = (n ^ (n >> 30)) * 0xbf58476d1ce4e5b9ULL;
        n = (n ^ (n >> 27)) * 0x94d049bb133111ebULL;
        return n ^ (n >> 31);
    }
    static void split_groups(Attempt& a) {
        auto old_ends = a.ends;
        const int old_groups = a.groups;
        a.groups = 0;
        int begin = 0;
        for (int g = 0; g < old_groups; ++g) {
            const int end = old_ends[g];
            std::sort(a.order.begin() + begin, a.order.begin() + end,
                      [&](int u, int v) { return a.colours[u] > a.colours[v]; });
            for (int i = begin + 1; i < end; ++i)
                if (a.colours[a.order[i]] != a.colours[a.order[i - 1]]) a.ends[a.groups++] = i;
            a.ends[a.groups++] = end;
            begin = end;
        }
    }
    static int compare(const Attempt& a, const Attempt& b) {
        if (a.groups != b.groups) return a.groups > b.groups ? 1 : -1;
        for (int i = 0; i < a.groups; ++i)
            if (a.ends[i] != b.ends[i]) return a.ends[i] > b.ends[i] ? 1 : -1;
        int begin = 0;
        for (int i = 0; i < a.groups; ++i) {
            const auto x = a.colours[a.order[begin]], y = b.colours[b.order[begin]];
            if (x != y) return x > y ? 1 : -1;
            begin = a.ends[i];
        }
        return 0;
    }
    static void refine(Attempt& a, int vertices, std::span<const Int> edges,
                       std::span<const Mask> components) {
        std::array<Mask, MaxVertices> next{};
        std::array<Mask, MaxLoop> next_components{};
        for (int v = 0; v < vertices; ++v) next[v] = hash(a.colours[v]);
        for (std::size_t e = 0; e < edges.size(); e += 2) {
            next[edges[e]] += a.colours[edges[e + 1]];
            next[edges[e + 1]] += a.colours[edges[e]];
        }
        for (std::size_t c = 0; c < components.size(); ++c) {
            next_components[c] = hash(a.component_colours[c]);
            auto mask = components[c];
            while (mask) {
                const int v = std::countr_zero(mask);
                next[v] += a.component_colours[c];
                next_components[c] += a.colours[v];
                mask &= mask - 1;
            }
        }
        a.colours = next;
        a.component_colours = next_components;
        split_groups(a);
    }
    static Entry materialize(const Attempt& a, int vertices, std::span<const Int> edges,
                             const Components& components) {
        std::array<Int, MaxVertices> labels{};
        for (int i = 0; i < vertices; ++i) labels[a.order[i]] = i;
        std::array<std::pair<Int, Int>, 3 * MaxLoop> pairs{};
        for (std::size_t e = 0; e < edges.size() / 2; ++e) {
            const auto u = labels[edges[2 * e]], v = labels[edges[2 * e + 1]];
            pairs[e] = {std::min(u, v), std::max(u, v)};
        }
        std::sort(pairs.begin(), pairs.begin() + edges.size() / 2);
        Entry result;
        result.endpoint_count = edges.size();
        for (std::size_t e = 0; e < edges.size() / 2; ++e) {
            result.canonical[2 * e] = pairs[e].first;
            result.canonical[2 * e + 1] = pairs[e].second;
        }
        std::array<Mask, MaxLoop> masks{};
        std::size_t count = 0;
        for (auto mask : components.masks()) {
            Mask renamed = 0;
            while (mask) {
                renamed |= Mask{1} << labels[std::countr_zero(mask)];
                mask &= mask - 1;
            }
            masks[count++] = renamed;
        }
        result.components.assign({masks.data(), count});
        return result;
    }

    static Orientation surviving_orientation(const std::vector<Attempt>& attempts) {
        bool forward = false, reverse = false;
        for (const auto& attempt : attempts) {
            forward |= !attempt.reversed;
            reverse |= attempt.reversed;
        }
        return forward ? (reverse ? Orientation::Both : Orientation::Keep) : Orientation::Discard;
    }

    static void keep_best(std::vector<Attempt>& attempts) {
        std::vector<Attempt> kept;
        for (auto& attempt : attempts) {
            const int cmp = kept.empty() ? 1 : compare(attempt, kept.front());
            if (cmp > 0) kept.clear();
            if (cmp >= 0) kept.push_back(std::move(attempt));
        }
        attempts = std::move(kept);
    }

    static Entry search(std::vector<Attempt> attempts, int vertices, std::span<const Int> edges,
                        const Components& components, Orientation& orientation) {
        // Both directions compete under the same weights. The watermark is
        // carried through branching, but never enters refinement or comparison.
        keep_best(attempts);
        if (orientation == Orientation::Both) {
            orientation = surviving_orientation(attempts);
            if (orientation == Orientation::Discard) return {};
        }
        while (attempts.front().groups < vertices) {
            for (int reload = 0; reload < std::max(1, vertices / 3); ++reload) {
                for (auto& attempt : attempts)
                    refine(attempt, vertices, edges, components.masks());
                keep_best(attempts);
                if (orientation == Orientation::Both) {
                    orientation = surviving_orientation(attempts);
                    if (orientation == Orientation::Discard) return {};
                }
            }
            if (attempts.front().groups == vertices) break;
            int begin = 0, end = 0;
            for (int g = 0; g < attempts.front().groups; ++g) {
                end = attempts.front().ends[g];
                if (end - begin > 1) break;
                begin = end;
            }
            std::vector<Attempt> children;
            for (const auto& attempt : attempts)
                for (int i = begin; i < end; ++i) {
                    children.push_back(attempt);
                    ++children.back().colours[attempt.order[i]];
                }
            attempts = std::move(children);
        }
        Entry best;
        for (const auto& attempt : attempts) {
            auto candidate = materialize(attempt, vertices, edges, components);
            const auto old_masks = best.components.masks();
            const auto new_masks = candidate.components.masks();
            // Direction decisions are finished. Select the representative and
            // relabel its masks without reversing the component stack.
            if (best.empty() || best.canonical < candidate.canonical
                || (best.canonical == candidate.canonical
                    && std::lexicographical_compare(old_masks.begin(), old_masks.end(),
                                                    new_masks.begin(), new_masks.end())))
                best = std::move(candidate);
        }
        return best;
    }

public:
    // An empty result means this construction has the noncanonical orientation.
    Entry standardize(int vertices, std::span<const Int> edges, const Components& components,
                      Orientation* orientation = nullptr) const {
        if (vertices < 3 || vertices > MaxVertices || edges.empty() || edges.size() % 2
            || edges.size() > 6 * MaxLoop)
            throw std::invalid_argument("invalid triangle graph dimensions");
        Attempt initial;
        std::array<int, MaxVertices> degrees{};
        for (Int v : edges) {
            if (v >= vertices) throw std::invalid_argument("invalid triangle endpoint");
            ++degrees[v];
        }
        Mask covered = 0;
        for (std::size_t c = 0; c < components.size(); ++c) {
            const auto mask = components.masks()[c];
            if (!mask || (mask >> vertices)) throw std::invalid_argument("invalid component mask");
            covered |= mask;
            initial.component_colours[c] = hash(std::popcount(mask));
        }
        if (covered != (Mask{1} << vertices) - 1) throw std::invalid_argument("uncovered vertex");
        // Compare the size sequence with its reversal, from the ends inward.
        // A single component has no reversal alternative.
        Orientation choice = components.size() == 1 ? Orientation::Keep : Orientation::Both;
        for (std::size_t c = 0; c < components.size() / 2; ++c) {
            const int first = std::popcount(components.masks()[c]);
            const int last = std::popcount(components.masks()[components.size() - 1 - c]);
            if (first == last) continue;
            choice = first > last ? Orientation::Keep : Orientation::Discard;
            break;
        }
        if (choice == Orientation::Discard) {
            if (orientation) *orientation = choice;
            return {};
        }
        for (int v = 0; v < vertices; ++v) {
            initial.colours[v] = hash(degrees[v]);
            initial.order[v] = v;
        }
        initial.groups = 1;
        initial.ends[0] = vertices;

        std::vector<Attempt> attempts;
        const int alternatives = choice == Orientation::Both ? 2 : 1;
        for (int reverse = 0; reverse < alternatives; ++reverse) {
            auto start = initial;
            start.reversed = reverse;
            std::array<std::vector<int>, MaxVertices> membership;
            // Only an undecided size sequence needs both directional partitions.
            for (std::size_t position = 0; position < components.size(); ++position) {
                const auto index = reverse ? components.size() - 1 - position : position;
                auto mask = components.masks()[index];
                while (mask) {
                    membership[std::countr_zero(mask)].push_back(position);
                    mask &= mask - 1;
                }
            }
            std::sort(start.order.begin(), start.order.begin() + vertices,
                      [&](int a, int b) { return membership[a] < membership[b]; });
            start.groups = 0;
            for (int i = 1; i < vertices; ++i)
                if (membership[start.order[i]] != membership[start.order[i - 1]])
                    start.ends[start.groups++] = i;
            start.ends[start.groups++] = vertices;
            split_groups(start);
            attempts.push_back(std::move(start));
        }
        auto result = search(std::move(attempts), vertices, edges, components, choice);
        if (orientation) *orientation = choice;
        return result;
    }
};

} // namespace GraphGeneration
