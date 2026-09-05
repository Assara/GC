// Small independent exhaustive oracle for the new pipeline.
#include "GraphGeneration/GraphGenerationPipeline.hpp"
#include "GraphGeneration/TransientToGCPipeline.hpp"
#include <bit>
#include <fstream>
#include <functional>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <unordered_map>

namespace {
#ifndef GC_VERIFY_MAX_LOOP
#define GC_VERIFY_MAX_LOOP 5
#endif
constexpr int max_loop = GC_VERIFY_MAX_LOOP;
constexpr int max_vertices = 2 * max_loop - 1;
static_assert(max_loop >= 3 && max_loop <= 6, "exhaustive oracle supports loop caps 3 through 6");
struct Record { std::uint64_t mask; bool odd, even; };
using Key = std::pair<int, int>;
std::map<Key, std::vector<Record>> generated;
std::map<Key, std::size_t> transient_counts;
std::map<Key, std::vector<std::uint64_t>> transient_records;

int edge_index(int n, int a, int b) {
    if (a > b) std::swap(a, b);
    return a * (2 * n - a - 1) / 2 + b - a - 1;
}

template <typename G>
std::uint64_t mask_of(const G& graph) {
    std::uint64_t mask = 0;
    for (Int e = 0; e < G::N_EDGES_; ++e) {
        auto [a, b] = graph.getEdge(e);
        if (a >= G::N_VERTICES_ || b >= G::N_VERTICES_ || a == b)
            throw std::runtime_error("invalid output endpoint");
        const auto bit = std::uint64_t{1} << edge_index(G::N_VERTICES_, a, b);
        if (mask & bit) throw std::runtime_error("parallel output edges");
        mask |= bit;
    }
    return mask;
}

template <int V, int L>
void convert(const std::filesystem::path& root) {
    constexpr int E = V - 1 + L;
    if constexpr (E <= V * (V - 1) / 2) {
        using G = Graph<V, E, 0, 0, 0, 0, fieldType>;
        const auto input = GraphGeneration::GraphGenerationPipeline<max_loop, max_vertices>::file_path(root / "transient", V, E);
        const auto output = root / "gc" / ("gc_V" + std::to_string(V) + "_E" + std::to_string(E) + ".gcg");
        GraphGeneration::MappedTransientGraph2Reader<G> transient(input);
        transient_counts[{V,E}] = transient.size();
        for (std::size_t i = 0; i < transient.size(); ++i) {
            transient_records[{V,E}].push_back(mask_of(transient[i]));
            const auto valences = transient[i].valence_array();
            if (!std::is_sorted(valences.begin(), valences.end()))
                throw std::runtime_error("unordered transient valences");
        }
        if (!std::filesystem::exists(output))
            GraphGeneration::TransientToGCPipeline<V, E>{}.run(input, output);
        GraphGeneration::MappedGCGraphReader<G> reader(output);
        auto& records = generated[{V,E}];
        std::size_t odd = 0, even = 0;
        for (std::size_t i = 0; i < reader.size(); ++i) {
            const auto record = reader[i];
            for (Int valence : record.graph.valence_array())
                if (valence < 3) throw std::runtime_error("GC output not at least trivalent");
            records.push_back({mask_of(record.graph), record.odd_gc, record.even_gc});
            odd += record.odd_gc; even += record.even_gc;
        }
        if (odd != reader.counts().odd_gc || even != reader.counts().even_gc)
            throw std::runtime_error("stored GC counts disagree with record flags");
    }
}

template <std::size_t V, std::size_t... L>
void convert_vertex(const std::filesystem::path& root, std::index_sequence<L...>) {
    (convert<V, L>(root), ...);
}

bool odd_permutation(const std::vector<int>& p) {
    bool odd = false;
    for (std::size_t i = 0; i < p.size(); ++i)
        for (std::size_t j = i + 1; j < p.size(); ++j) odd ^= p[i] > p[j];
    return odd;
}

// Enumerate labeled graphs by assigning each vertex's higher-labeled neighbors.
// Prune only by minimum valence 3 and the loop cap. Every simple graph appears
// once. Isomorphism classes and signs use full permutation orbits, independently
// of the repository's canonicalization implementation.
struct Oracle {
    int n;
    std::vector<int> degree;
    std::vector<std::pair<int,int>> pairs;
    std::unordered_map<std::uint64_t, std::size_t> orbit;
    std::vector<Record> representatives;
    std::size_t labeled = 0;

    explicit Oracle(int vertices) : n(vertices), degree(n) {
        for (int a = 0; a < n; ++a)
            for (int b = a + 1; b < n; ++b) pairs.emplace_back(a,b);
    }
    bool connected(std::uint64_t mask) const {
        unsigned seen = 1;
        bool changed;
        do {
            const auto previous = seen;
            for (std::size_t i = 0; i < pairs.size(); ++i) if (mask & (std::uint64_t{1} << i)) {
                auto [a,b] = pairs[i];
                if (seen & ((1U << a) | (1U << b))) seen |= (1U << a) | (1U << b);
            }
            changed = seen != previous;
        } while (changed);
        return seen == (1U << n) - 1;
    }
    // Only equal-degree vertices can be exchanged after ordering by degree.
    // Keep one exhaustive canonical representative, rather than storing all
    // labeled images of every transient graph.
    std::uint64_t canonical(std::uint64_t mask) const {
        std::vector<int> valences(n), order(n), labels(n);
        std::vector<std::pair<int,int>> edges;
        for (std::size_t i=0; i<pairs.size(); ++i) if (mask & (std::uint64_t{1} << i)) {
            const auto [a,b] = pairs[i];
            ++valences[a]; ++valences[b]; edges.push_back(pairs[i]);
        }
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&](int a,int b) { return valences[a]<valences[b]; });
        std::uint64_t best = ~std::uint64_t{0};
        auto enumerate = [&](auto&& self, int first) -> void {
            if (first == n) {
                for (int i=0; i<n; ++i) labels[order[i]] = i;
                std::uint64_t image = 0;
                for (auto [a,b] : edges)
                    image |= std::uint64_t{1} << edge_index(n,labels[a],labels[b]);
                best = std::min(best,image);
                return;
            }
            int last = first+1;
            while (last<n && valences[order[last]]==valences[order[first]]) ++last;
            do { self(self,last); }
            while (std::next_permutation(order.begin()+first,order.begin()+last));
        };
        enumerate(enumerate,0);
        return best;
    }
    void accept(std::uint64_t mask) {
        if (!connected(mask)) return;
        ++labeled;
        if (orbit.contains(mask)) return;
        const auto id = representatives.size();
        Record record{mask, true, true};
        std::vector<int> p(n), edge_ids(pairs.size(), -1);
        std::vector<std::pair<int,int>> edges;
        for (std::size_t i = 0; i < pairs.size(); ++i) if (mask & (std::uint64_t{1} << i)) {
            edge_ids[i] = edges.size();
            edges.push_back(pairs[i]);
        }
        std::iota(p.begin(), p.end(), 0);
        do {
            std::uint64_t image = 0;
            for (auto [a,b] : edges) image |= std::uint64_t{1} << edge_index(n,p[a],p[b]);
            orbit.emplace(image, id);
            if (image == mask) {
                bool odd_sign = odd_permutation(p);
                for (auto [a,b] : edges) odd_sign ^= p[a] > p[b];
                record.odd &= !odd_sign;
                std::vector<int> edge_permutation;
                for (auto [a,b] : edges) edge_permutation.push_back(edge_ids[edge_index(n,p[a],p[b])]);
                record.even &= !odd_permutation(edge_permutation);
            }
        } while (std::next_permutation(p.begin(), p.end()));
        representatives.push_back(record);
    }
    void enumerate(int vertex = 0, int edges = 0, std::uint64_t mask = 0) {
        if (edges > n - 1 + max_loop) return;
        int deficit = 0;
        for (int i = vertex; i < n; ++i) {
            if (degree[i] + n - vertex - 1 < 3) return;
            deficit += std::max(0, 3 - degree[i]);
        }
        if (edges + (deficit + 1) / 2 > n - 1 + max_loop) return;
        if (vertex == n) { accept(mask); return; }
        const int remaining = n - vertex - 1;
        for (unsigned subset = 0; subset < (1U << remaining); ++subset) {
            const int added = std::popcount(subset);
            if (degree[vertex] + added < 3 || edges + added > n - 1 + max_loop) continue;
            auto next = mask;
            for (int j = 0; j < remaining; ++j) if (subset & (1U << j)) {
                const int b = vertex + 1 + j;
                ++degree[b];
                next |= std::uint64_t{1} << edge_index(n,vertex,b);
            }
            enumerate(vertex + 1, edges + added, next);
            for (int j = 0; j < remaining; ++j) if (subset & (1U << j)) --degree[vertex+1+j];
        }
    }
};

std::string edges_text(int n, std::uint64_t mask) {
    std::string text;
    for (int a=0;a<n;++a) for (int b=a+1;b<n;++b)
        if (mask & (std::uint64_t{1} << edge_index(n,a,b)))
            text += "(" + std::to_string(a) + "," + std::to_string(b) + ")";
    return text;
}
} // namespace

int main(int argc, char** argv) {
    if (argc != 2) { std::cerr << "Usage: " << argv[0] << " RUN_DIRECTORY\n"; return 2; }
    try {
        const std::filesystem::path root(argv[1]);
        [&]<std::size_t... I>(std::index_sequence<I...>) {
            (convert_vertex<I+3>(root,std::make_index_sequence<max_loop + 1>{}), ...);
        }(std::make_index_sequence<max_vertices - 2>{});
        std::ofstream table(root / "counts.tsv");
        std::ofstream missing(root / "missing_graphs.tsv");
        table << "loop\tvertices\tedges\ttransient\ttotal\toddGC\tevenGC\treference_total\treference_oddGC\treference_evenGC\tmissing\tunexpected\tsign_mismatches\n";
        missing << "loop\tvertices\tedges\toddGC\tevenGC\tedge_list\n";
        std::size_t missing_total = 0, unexpected_total = 0, sign_errors = 0;
        for (int v=3;v<=max_vertices;++v) {
            Oracle oracle(v);
            if (v >= 4 && v <= 2 * (max_loop - 1)) oracle.enumerate();
            std::cout << "Reference V=" << v << " labeled=" << oracle.labeled
                      << " classes=" << oracle.representatives.size() << std::endl;
            for (int l=0;l<=max_loop;++l) {
                const int e=v-1+l;
                if (e>v*(v-1)/2) continue;
                Oracle transient_oracle(v);
                std::set<std::uint64_t> unique_transients;
                for (auto mask : transient_records[{v,e}]) {
                    if (!transient_oracle.connected(mask))
                        throw std::runtime_error("disconnected transient graph");
                    if (!unique_transients.insert(transient_oracle.canonical(mask)).second)
                        throw std::runtime_error("isomorphic duplicate transient graphs");
                }
                const auto& actual=generated[{v,e}];
                std::set<std::size_t> found;
                std::size_t odd=0,even=0,ref_total=0,ref_odd=0,ref_even=0,unexpected=0,errors=0,miss=0;
                for (const auto& record : actual) {
                    odd += record.odd; even += record.even;
                    auto it=oracle.orbit.find(record.mask);
                    if (it == oracle.orbit.end() || !found.insert(it->second).second) { ++unexpected; continue; }
                    const auto& ref=oracle.representatives[it->second];
                    errors += record.odd != ref.odd || record.even != ref.even;
                }
                for (std::size_t i=0;i<oracle.representatives.size();++i) {
                    const auto& ref=oracle.representatives[i];
                    if (std::popcount(ref.mask)!=e) continue;
                    ++ref_total; ref_odd += ref.odd; ref_even += ref.even;
                    if (!found.contains(i)) {
                        ++miss;
                        missing << l << '\t' << v << '\t' << e << '\t' << ref.odd << '\t' << ref.even
                                << '\t' << edges_text(v,ref.mask) << '\n';
                    }
                }
                table << l << '\t' << v << '\t' << e << '\t' << transient_counts[{v,e}] << '\t'
                      << actual.size() << '\t' << odd << '\t' << even << '\t' << ref_total << '\t'
                      << ref_odd << '\t' << ref_even << '\t' << miss << '\t' << unexpected << '\t' << errors << '\n';
                missing_total+=miss; unexpected_total+=unexpected; sign_errors+=errors;
            }
        }
        std::cout << "missing=" << missing_total << " unexpected_or_duplicate=" << unexpected_total
                  << " sign_mismatches=" << sign_errors << '\n';
        return missing_total || unexpected_total || sign_errors ? 1 : 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 2; }
}
