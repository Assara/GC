#include "GraphGeneration/TransientGraph2Standardizer.hpp"

#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <set>
#include <string>
#include <type_traits>
#include <tuple>
#include <vector>

namespace {
using Edge = std::pair<Int, Int>;

// Storage is a caller policy, independent of transient_graph2.
template <typename T>
struct bucket {
    template <std::size_t... K>
    static auto group_types(std::index_sequence<K...>)
        -> std::tuple<std::vector<typename T::template split_graph_type<K>>...>;
    using groups_type = decltype(group_types(
        std::make_index_sequence<T::MAX_SHARED_EDGES + 1>{}));
    groups_type groups;

    bucket() = default;
    bucket(const bucket&) = delete;
    bucket& operator=(const bucket&) = delete;
    bucket(bucket&&) = default;

    template <typename G>
    void add(G&& child) {
        constexpr auto k = std::remove_cvref_t<G>::N_EDGES_ - T::N_EDGES_ - 1;
        std::get<k>(groups).push_back(std::forward<G>(child));
    }
};

template <typename T>
auto collect_splits(const T& input, Int vertex, Int preserve_valence = 2, Int min_valence = 2,
                    Int max_loop_number = std::numeric_limits<Int>::max()) {
    bucket<T> result;
    input.split(vertex, preserve_valence, min_valence, max_loop_number, result);
    return result;
}

struct counting_bucket {
    std::size_t total = 0;
    template <typename G>
    void add(G&&) { ++total; }
};

std::string key(std::vector<Edge> edges, Int old_vertex,
                Int new_vertex, bool swap) {
    auto label = [&](Int v) {
        if (!swap) return v;
        return v == old_vertex ? new_vertex : v == new_vertex ? old_vertex : v;
    };
    for (auto& [a, b] : edges) {
        a = label(a);
        b = label(b);
        if (a > b) std::swap(a, b);
    }
    std::sort(edges.begin(), edges.end());
    std::string result;
    for (const auto& [a, b] : edges)
        result += std::to_string(a) + "," + std::to_string(b) + ";";
    return result;
}

std::string canonical_key(const std::vector<Edge>& edges, Int v, Int w) {
    return std::min(key(edges, v, w, false), key(edges, v, w, true));
}

// Independent exhaustive assignments, then quotient by exchanging split vertices.
template <typename T>
auto reference(const T& input, Int vertex, Int preserve_valence, Int min_valence, Int max_loop_number) {
    std::vector<std::set<std::string>> result(T::MAX_SHARED_EDGES + 1);
    std::vector<Edge> fixed;
    std::vector<Int> neighbors;
    for (Int e = 0; e < T::N_EDGES_; ++e) {
        auto [a, b] = input.graph().getEdge(e);
        if (a == vertex) neighbors.push_back(b);
        else if (b == vertex) neighbors.push_back(a);
        else fixed.emplace_back(a, b);
    }
    std::size_t assignments = 1;
    for (std::size_t i = 0; i < neighbors.size(); ++i) assignments *= 3;
    for (std::size_t code = 0; code < assignments; ++code) {
        auto edges = fixed;
        std::size_t remaining = code, shared = 0;
        for (Int neighbor : neighbors) {
            const auto choice = remaining % 3;
            remaining /= 3;
            if (choice != 1) edges.emplace_back(vertex, neighbor);
            if (choice != 0) edges.emplace_back(T::N_VERTICES_, neighbor);
            shared += choice == 2;
        }
        edges.emplace_back(vertex, T::N_VERTICES_);
        if (static_cast<int>(edges.size()) - T::N_VERTICES_ > max_loop_number) continue;
        int left = 0, right = 0;
        for (auto [a, b] : edges) {
            left += (a == vertex) + (b == vertex);
            right += (a == T::N_VERTICES_) + (b == T::N_VERTICES_);
        }
        const auto parent_valence = neighbors.size();
        if (neighbors.size() != T::N_VERTICES_ - 1
            && std::max(left, right) > parent_valence) continue;
        if (left >= min_valence && right >= min_valence && std::max(left, right) >= preserve_valence)
            result[shared].insert(canonical_key(edges, vertex, T::N_VERTICES_));
    }
    return result;
}

template <typename T>
void check_splits(const T& input, Int vertex, Int preserve_valence = 2, Int min_valence = 2,
                    Int max_loop_number = std::numeric_limits<Int>::max()) {
    const auto expected = reference(input, vertex, preserve_valence, min_valence, max_loop_number);
    const auto result = collect_splits(input, vertex, preserve_valence, min_valence, max_loop_number);
    std::size_t shared = 0;
    std::apply([&](const auto&... groups) {
        auto check_group = [&](const auto& children) {
            using G = typename std::decay_t<decltype(children)>::value_type;
            static_assert(G::N_HAIR == 0);
            assert(G::N_VERTICES_ == T::N_VERTICES_ + 1);
            assert(G::N_EDGES_ == T::N_EDGES_ + 1 + shared);
            std::set<std::string> actual;
            for (const auto& child : children) {
                std::vector<Edge> edges;
                std::set<Edge> unique_edges;
                for (Int e = 0; e < G::N_EDGES_; ++e) {
                    auto [a, b] = child.getEdge(e);
                    assert(a < G::N_VERTICES_ && b < G::N_VERTICES_ && a != b);
                    edges.emplace_back(a, b);
                    assert(unique_edges.emplace(std::min(a, b), std::max(a, b)).second);
                }
                assert(actual.insert(canonical_key(edges,
                                                   vertex, T::N_VERTICES_)).second);
            }
            assert(actual == expected[shared]);
            ++shared;
        };
        (check_group(groups), ...);
    }, result.groups);
}

template <typename Collection>
std::size_t count(const Collection& collection) {
    return std::apply([](const auto&... groups) { return (groups.size() + ... + 0); }, collection.groups);
}

template <typename T>
void check_all(const T& input) {
    for (Int v = 0; v < T::N_VERTICES_; ++v) {
        for (Int minimum = 1; minimum <= 3; ++minimum)
            for (Int threshold = 0; threshold <= T::N_VERTICES_ + 2; ++threshold)
                for (Int cap = 0; cap <= T::N_EDGES_ + 1; ++cap)
                    check_splits(input, v, threshold, minimum, cap);
        check_splits(input, v, std::numeric_limits<Int>::max());
    }
    bucket<T> appended;
    std::size_t expected = 0;
    for (Int v = 0; v < T::N_VERTICES_; ++v) {
        expected += count(collect_splits(input, v));
        input.split(v, 2, 2, std::numeric_limits<Int>::max(), appended);
    }
    assert(count(appended) == expected);
    counting_bucket counter;
    for (int repeat = 0; repeat < 2; ++repeat) {
        for (Int v = 0; v < T::N_VERTICES_; ++v) {
            input.split(v, 2, 2, std::numeric_limits<Int>::max(), counter);
        }
        assert(counter.total == (repeat + 1) * expected);
    }
    GraphGeneration::transient_graph2_standardizer<T::N_VERTICES_, T::N_EDGES_> standardizer;
    const auto standardized = standardizer.standardize_no_sign(input);
    check_splits(standardized, 0);
}
} // namespace

int main() {
    using Star = GraphGeneration::transient_graph2<5, 4>;
    Star::graph_type star;
    for (Int e = 0; e < 4; ++e) star.setEdge(e, 0, e + 1);
    check_all(Star(star));
    // All incident edges may be shared, including the first one.
    const auto all_shared = collect_splits(Star(star), 0);
    assert(std::get<4>(all_shared.groups).size() == 1);

    using Complete = GraphGeneration::transient_graph2<4, 6>;
    Complete::graph_type complete;
    Int e = 0;
    for (Int a = 0; a < 4; ++a)
        for (Int b = a + 1; b < 4; ++b) complete.setEdge(e++, a, b);
    check_all(Complete(complete));

    using Path = GraphGeneration::transient_graph2<3, 2>;
    Path::graph_type path;
    path.setEdge(0, 0, 1);
    path.setEdge(1, 2, 1); // Include incidences at either endpoint position.
    check_all(Path(path));
    assert(std::get<2>(collect_splits(Path(path), 1).groups).size() == 1);
    // Separating the two path edges makes both new vertices bivalent.
    const auto middle = collect_splits(Path(path), 1);
    assert(std::get<0>(middle.groups).size() == 1);
    const auto middle_valences = std::get<0>(middle.groups).front().valence_array();
    assert(middle_valences[1] == 2 && middle_valences[3] == 2);
    // A non-universal leaf cannot increase its valence by sharing its edge.
    assert(count(collect_splits(Path(path), 2)) == 0);
    assert(count(collect_splits(Path(path), 0)) == 0);

    // A zero-loop cap allows ordinary tree splits but no shared edges.
    const auto tree_splits = collect_splits(Path(path), 1, 2, 2, 0);
    assert(std::get<0>(tree_splits.groups).size() == 1);
    assert(std::get<1>(tree_splits.groups).empty());
    assert(std::get<2>(tree_splits.groups).empty());
    const auto one_loop = collect_splits(Path(path), 1, 2, 2, 1);
    assert(std::get<1>(one_loop.groups).size() == 2);
    assert(std::get<2>(one_loop.groups).empty());
    // K4 already has three loops, so a cap of two rejects all children.
    assert(count(collect_splits(Complete(complete), 0, 2, 2, 2)) == 0);

    // Threshold three keeps (2,3), (3,2), and (3,3), but rejects (2,2).
    const auto preserved = collect_splits(Path(path), 1, 3);
    assert(std::get<0>(preserved.groups).empty());
    assert(std::get<1>(preserved.groups).size() == 2);
    assert(std::get<2>(preserved.groups).size() == 1);
    assert(count(collect_splits(Path(path), 1, 4)) == 0);

    // With the same degree but no universal vertex, only valence-preserving splits survive.
    using LongerPath = GraphGeneration::transient_graph2<4, 3>;
    LongerPath::graph_type longer;
    longer.setEdge(0, 0, 1);
    longer.setEdge(1, 1, 2);
    longer.setEdge(2, 2, 3);
    check_all(LongerPath(longer));
    assert(count(collect_splits(LongerPath(longer), 1)) == 1);
    assert(count(collect_splits(LongerPath(longer), 1, 3)) == 0);

    using Empty = GraphGeneration::transient_graph2<1, 0>;
    Empty::graph_type empty;
    check_all(Empty(empty));
    std::cout << "transient_graph2 tests passed\n";
}
