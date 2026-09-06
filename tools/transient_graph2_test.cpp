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
auto collect_splits(const T& input, Int vertex, Int min_valence = 2,
                    Int min_loop_number = 0, Int max_loop_number = std::numeric_limits<Int>::max(),
                    Int preserve_valence = 0, bool require_root = false,
                    Int max_valence = std::numeric_limits<Int>::max()) {
    bucket<T> result;
    if (require_root)
        input.split_with_bivalent_neighbours(vertex, preserve_valence, min_valence, max_valence,
                                             min_loop_number, max_loop_number, result);
    else
        input.split(vertex, preserve_valence, min_valence, max_valence, min_loop_number, max_loop_number, result);
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
auto reference(const T& input, Int vertex, Int min_valence, Int min_loop_number, Int max_loop_number, Int preserve_valence, bool require_root, Int max_valence) {
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
        const int loops = static_cast<int>(edges.size()) - T::N_VERTICES_;
        if (loops < min_loop_number || loops > max_loop_number) continue;
        if (require_root) {
            // Independent oracle: build adjacency sets and intersect neighbours of degree-two vertices.
            std::array<std::set<Int>, T::N_VERTICES_ + 1> adjacent;
            for (auto [a, b] : edges) {
                adjacent[a].insert(b);
                adjacent[b].insert(a);
            }
            if (adjacent[vertex].size() != 2 && adjacent[T::N_VERTICES_].size() != 2
                && std::ranges::any_of(adjacent, [](const auto& n) { return n.size() == 2; }))
                continue;
            bool has_root = false;
            for (Int candidate = 0; candidate <= T::N_VERTICES_; ++candidate) {
                if (adjacent[candidate].size() == 2) continue;
                bool supports_all = true;
                for (const auto& neighbors : adjacent)
                    if (neighbors.size() == 2 && !neighbors.contains(candidate)) supports_all = false;
                has_root |= supports_all;
            }
            if (!has_root) continue;
        }
        int left = 0, right = 0;
        for (auto [a, b] : edges) {
            left += (a == vertex) + (b == vertex);
            right += (a == T::N_VERTICES_) + (b == T::N_VERTICES_);
        }
        if (left >= min_valence && right >= min_valence
            && left <= max_valence && right <= max_valence
            && std::max(left, right) >= preserve_valence)
            result[shared].insert(canonical_key(edges, vertex, T::N_VERTICES_));
    }
    return result;
}

template <typename T>
void check_splits(const T& input, Int vertex, Int min_valence = 2,
                    Int min_loop_number = 0, Int max_loop_number = std::numeric_limits<Int>::max(),
                    Int preserve_valence = 0, bool require_root = false,
                    Int max_valence = std::numeric_limits<Int>::max()) {
    const auto expected = reference(input, vertex, min_valence, min_loop_number, max_loop_number, preserve_valence, require_root, max_valence);
    const auto result = collect_splits(input, vertex, min_valence, min_loop_number, max_loop_number, preserve_valence, require_root, max_valence);
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
            for (Int cap = 0; cap <= T::N_EDGES_ + 1; ++cap)
                for (Int floor = 0; floor <= T::N_EDGES_ + 2; ++floor)
                    for (Int preserve = 0; preserve <= T::N_VERTICES_ + 2; ++preserve)
                        for (bool require_root : {false, true})
                            check_splits(input, v, minimum, floor, cap, preserve, require_root);
        check_splits(input, v, std::numeric_limits<Int>::max());
        for (Int maximum = 0; maximum <= T::N_VERTICES_ + 1; ++maximum)
            for (bool require_root : {false, true})
                check_splits(input, v, 2, 0, T::N_EDGES_ + 1, 3, require_root, maximum);
    }
    bucket<T> appended;
    std::size_t expected = 0;
    for (Int v = 0; v < T::N_VERTICES_; ++v) {
        expected += count(collect_splits(input, v));
        input.split(v, 0, 2, std::numeric_limits<Int>::max(), 0, std::numeric_limits<Int>::max(), appended);
    }
    assert(count(appended) == expected);
    counting_bucket counter;
    for (int repeat = 0; repeat < 2; ++repeat) {
        for (Int v = 0; v < T::N_VERTICES_; ++v) {
            input.split(v, 0, 2, std::numeric_limits<Int>::max(), 0, std::numeric_limits<Int>::max(), counter);
        }
        assert(counter.total == (repeat + 1) * expected);
    }
    GraphGeneration::transient_graph2_standardizer<T::N_VERTICES_, T::N_EDGES_> standardizer;
    const auto standardized = standardizer.standardize_no_sign(input);
    check_splits(standardized, 0);
}
} // namespace

int main() {
    using Triangle = GraphGeneration::transient_graph2<3, 3>;
    Triangle::graph_type triangle;
    triangle.setEdge(0, 0, 1); triangle.setEdge(1, 1, 2); triangle.setEdge(2, 0, 2);
    const auto rooted = collect_splits(Triangle(triangle), 0, 2, 0, 3, 2, true);
    assert(std::get<0>(rooted.groups).empty()); // C4 has no eligible root.
    assert(std::get<1>(rooted.groups).size() == 2); // Diamond retains two roots.
    assert(std::get<2>(rooted.groups).size() == 1); // K4 has no bivalents.
    check_all(Triangle(triangle));
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
    // Extended edges retain the whole parent and add one two-edge path per incidence.
    for (Int vertex = 0; vertex < 3; ++vertex) {
        bucket<Path> extensions;
        Path(path).create_bivalent_vertices(vertex, 0, std::numeric_limits<Int>::max(), 1, 1, extensions);
        const auto& children = std::get<1>(extensions.groups);
        assert(children.size() == path.adjacent(vertex).size());
        std::set<Int> extended_neighbours;
        for (const auto& child : children) {
            assert(std::equal(path.half_edges.begin(), path.half_edges.end(), child.half_edges.begin()));
            const auto [root, added] = child.getEdge(2);
            const auto [same_added, neighbour] = child.getEdge(3);
            assert(root == vertex && added == 3 && same_added == 3);
            assert(child.valence_array()[3] == 2);
            assert(extended_neighbours.insert(neighbour).second);
            bool retained_edge = false;
            for (Int edge = 0; edge < 2; ++edge) {
                const auto [a, b] = path.getEdge(edge);
                retained_edge |= std::minmax(a, b) == std::minmax(vertex, neighbour);
            }
            assert(retained_edge);
        }
        assert(count(extensions) == children.size());
    }
    counting_bucket rejected_extensions;
    Path(path).create_bivalent_vertices(1, 0, std::numeric_limits<Int>::max(), 0, 0, rejected_extensions); // Above maximum.
    Path(path).create_bivalent_vertices(1, 0, std::numeric_limits<Int>::max(), 2, 3, rejected_extensions); // Below minimum.
    Path(path).create_bivalent_vertices(1, 0, std::numeric_limits<Int>::max(), 2, 1, rejected_extensions); // Inverted range.
    Complete(complete).create_bivalent_vertices(0, 0, std::numeric_limits<Int>::max(), 0, 3, rejected_extensions); // L=4 required.
    assert(rejected_extensions.total == 0);
    // Inclusive valence bounds apply to the chosen vertex after the addition.
    bucket<Path> at_bound;
    Path(path).create_bivalent_vertices(1, 3, 3, 1, 1, at_bound);
    assert(count(at_bound) == 2);
    Path(path).create_bivalent_vertices(1, 4, 5, 1, 1, rejected_extensions);
    Path(path).create_bivalent_vertices(1, 0, 2, 1, 1, rejected_extensions);
    Path(path).create_bivalent_vertices(1, 4, 3, 1, 1, rejected_extensions);
    assert(rejected_extensions.total == 0);
    bucket<Path> chosen_only;
    Path(path).create_bivalent_vertices(0, 2, 2, 1, 1, chosen_only);
    assert(count(chosen_only) == 1);
    const auto degrees = std::get<1>(chosen_only.groups).front().valence_array();
    assert(degrees[0] == 2 && degrees[1] == 3); // Neighbour may exceed the chosen-vertex cap.
    check_all(Path(path));
    assert(std::get<2>(collect_splits(Path(path), 1).groups).size() == 1);
    // Creating (3,3) while leaving the original leaves bivalent is forbidden
    // by the specialized method, even though the diamond has eligible roots.
    const auto bivalent_rule = collect_splits(Path(path), 1, 2, 0, 2, 0, true);
    assert(std::get<2>(bivalent_rule.groups).empty());
    assert(std::get<1>(bivalent_rule.groups).size() == 2); // New bivalent endpoint.

    // Minimum split valence three still allows original leaves to become
    // bivalent by sharing. Only the all-shared split of K1,2 survives.
    const auto trivalent = collect_splits(Path(path), 1, 3);
    assert(std::get<0>(trivalent.groups).empty());
    assert(std::get<1>(trivalent.groups).empty());
    assert(std::get<2>(trivalent.groups).size() == 1);
    const auto trivalent_valences = std::get<2>(trivalent.groups).front().valence_array();
    assert(trivalent_valences[1] == 3 && trivalent_valences[3] == 3);
    assert(trivalent_valences[0] == 2 && trivalent_valences[2] == 2);
    // Separating the two path edges makes both new vertices bivalent.
    const auto middle = collect_splits(Path(path), 1);
    assert(std::get<0>(middle.groups).size() == 1);
    const auto middle_valences = std::get<0>(middle.groups).front().valence_array();
    assert(middle_valences[1] == 2 && middle_valences[3] == 2);
    // Sharing a leaf edge creates two bivalent split vertices.
    assert(count(collect_splits(Path(path), 2)) == 1);
    assert(count(collect_splits(Path(path), 0)) == 1);

    // A valence cap of two accepts (2,2), but rejects (2,3) and (3,3).
    const auto capped = collect_splits(Path(path), 1, 2, 0, 2, 0, false, 2);
    assert(std::get<0>(capped.groups).size() == 1);
    assert(std::get<1>(capped.groups).empty());
    assert(std::get<2>(capped.groups).empty());
    assert(count(collect_splits(Path(path), 1, 2, 0, 2, 3, false, 2)) == 0);

    // A zero-loop cap allows ordinary tree splits but no shared edges.
    const auto tree_splits = collect_splits(Path(path), 1, 2, 0, 0);
    assert(std::get<0>(tree_splits.groups).size() == 1);
    assert(std::get<1>(tree_splits.groups).empty());
    assert(std::get<2>(tree_splits.groups).empty());
    const auto one_loop = collect_splits(Path(path), 1, 2, 0, 1);
    assert(std::get<1>(one_loop.groups).size() == 2);
    assert(std::get<2>(one_loop.groups).empty());
    // K4 already has three loops, so a cap of two rejects all children.
    assert(count(collect_splits(Complete(complete), 0, 2, 0, 2)) == 0);

    // Exact loop order requires exactly that many duplicated edges for this tree.
    const auto exact_loop = collect_splits(Path(path), 1, 2, 1, 1);
    assert(std::get<0>(exact_loop.groups).empty());
    assert(std::get<1>(exact_loop.groups).size() == 2);
    assert(std::get<2>(exact_loop.groups).empty());
    assert(count(collect_splits(Path(path), 1, 2, 3, 4)) == 0); // Unreachable minimum.
    assert(count(collect_splits(Path(path), 1, 2, 2, 1)) == 0); // Empty interval.
    assert(count(collect_splits(Complete(complete), 0, 2, 2, 6))
           == count(collect_splits(Complete(complete), 0, 2, 0, 6)));
    // Preserve applies to the split pair, and equality is accepted.
    const auto preserved = collect_splits(Path(path), 1, 2, 0, 2, 3);
    assert(std::get<0>(preserved.groups).empty()); // (2,2) rejected.
    assert(std::get<1>(preserved.groups).size() == 2); // (2,3) accepted.
    assert(std::get<2>(preserved.groups).size() == 1); // (3,3) accepted.
    assert(count(collect_splits(Path(path), 1, 2, 0, 2, 4)) == 0);
    assert(count(collect_splits(Path(path), 1, 4)) == 0);

    // A non-universal vertex may also increase its valence.
    using LongerPath = GraphGeneration::transient_graph2<4, 3>;
    LongerPath::graph_type longer;
    longer.setEdge(0, 0, 1);
    longer.setEdge(1, 1, 2);
    longer.setEdge(2, 2, 3);
    check_all(LongerPath(longer));
    assert(count(collect_splits(LongerPath(longer), 1)) == 4);
    assert(count(collect_splits(LongerPath(longer), 1, 3)) == 1);

    using Empty = GraphGeneration::transient_graph2<1, 0>;
    Empty::graph_type empty;
    counting_bucket empty_extensions;
    Empty(empty).create_bivalent_vertices(0, 0, std::numeric_limits<Int>::max(), 0, 3, empty_extensions);
    assert(empty_extensions.total == 0);
    check_all(Empty(empty));
    std::cout << "transient_graph2 tests passed\n";
}
