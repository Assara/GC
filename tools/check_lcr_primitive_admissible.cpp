#include <iostream>
#include <optional>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "coalgebra_utils.hpp"
#include "examplegraphs.hpp"

namespace {

using FullGC = GC<13, 24, 0, 0, 0, 1>;
using ContGC = GC<12, 23, 0, 0, 0, 1>;
using FullGraph = typename FullGC::GraphType;

std::vector<std::vector<bool>> all_sequences(Int length) {
	std::vector<std::vector<bool>> result;
	const std::size_t count = std::size_t{1} << length;
	result.reserve(count);

	for (std::size_t mask = 0; mask < count; ++mask) {
		std::vector<bool> sequence;
		sequence.reserve(length);
		for (Int i = 0; i < length; ++i) {
			sequence.push_back(((mask >> i) & 1U) != 0U);
		}
		result.push_back(std::move(sequence));
	}

	return result;
}

std::string sequence_string(const std::vector<bool>& sequence) {
	std::string result;
	for (bool bit : sequence) {
		result.push_back(bit ? 'R' : 'L');
	}
	return result.empty() ? "empty" : result;
}

FullGC build_sum(Int first_length, Int second_length) {
	FullGC result;
	const auto first_sequences = all_sequences(first_length);
	const auto second_sequences = all_sequences(second_length);

	for (const auto& first : first_sequences) {
		for (const auto& second : second_sequences) {
			std::vector<std::vector<bool>> blocks{first, second};
			result += FullGC(iterated_V_graph_lcr_extension<12>(blocks));
		}
	}

	result.standardize_all();
	result.sort_elements();
	return result;
}

std::string signature_from_half_edges(Int n_vertices, const auto& half_edges) {
	std::string result = std::to_string(static_cast<int>(n_vertices)) + ":";
	for (std::size_t i = 0; i < half_edges.size(); ++i) {
		if (i > 0) {
			result.push_back(',');
		}
		result += std::to_string(static_cast<int>(half_edges[i]));
	}
	return result;
}

template <Int N_VERTICES, Int N_EDGES>
std::string build_canonical_signature_from_selected_edges(
	const FullGraph& graph,
	const std::vector<Int>& selected_edges,
	const std::vector<Int>& relabel
) {
	using SubGraph = Graph<N_VERTICES, N_EDGES, 0, 0, 0, 1, fieldType>;
	SubGraph subgraph;

	for (Int i = 0; i < N_EDGES; ++i) {
		const auto [u, v] = graph.getEdge(selected_edges[i]);
		subgraph.setEdge(i, relabel[u], relabel[v]);
	}

	typename SubGraph::Basis be(subgraph);
	const auto canonical = SubGraph::canonized(be).getValue();
	return signature_from_half_edges(N_VERTICES, canonical.half_edges);
}

std::optional<std::string> admissible_signature(
	const FullGraph& graph,
	const std::vector<Int>& selected_edges
) {
	const Int n_edges = static_cast<Int>(selected_edges.size());
	if (n_edges == 0 || n_edges == FullGraph::N_EDGES_) {
		return std::nullopt;
	}

	std::array<bool, FullGraph::N_VERTICES_> used{};
	std::array<Int, FullGraph::N_VERTICES_> degrees{};
	std::array<std::vector<Int>, FullGraph::N_VERTICES_> adjacency{};

	for (Int edge_index : selected_edges) {
		const auto [u, v] = graph.getEdge(edge_index);
		used[u] = true;
		used[v] = true;
		++degrees[u];
		++degrees[v];
		adjacency[u].push_back(v);
		adjacency[v].push_back(u);
	}

	Int n_vertices = 0;
	Int start = FullGraph::N_VERTICES_;
	for (Int v = 0; v < FullGraph::N_VERTICES_; ++v) {
		if (!used[v]) {
			continue;
		}
		++n_vertices;
		if (degrees[v] < 3) {
			return std::nullopt;
		}
		if (start == FullGraph::N_VERTICES_) {
			start = v;
		}
	}

	if (n_edges - 2 * n_vertices + 2 != 0) {
		return std::nullopt;
	}

	std::array<bool, FullGraph::N_VERTICES_> visited{};
	std::vector<Int> stack{start};
	visited[start] = true;
	Int reached = 0;
	while (!stack.empty()) {
		const Int v = stack.back();
		stack.pop_back();
		++reached;
		for (Int w : adjacency[v]) {
			if (!visited[w]) {
				visited[w] = true;
				stack.push_back(w);
			}
		}
	}

	if (reached != n_vertices) {
		return std::nullopt;
	}

	std::vector<Int> relabel(FullGraph::N_VERTICES_, FullGraph::N_VERTICES_);
	Int next = 0;
	for (Int v = 0; v < FullGraph::N_VERTICES_; ++v) {
		if (used[v]) {
			relabel[v] = next++;
		}
	}

	switch (n_vertices) {
		case 2: return build_canonical_signature_from_selected_edges<2, 2>(graph, selected_edges, relabel);
		case 3: return build_canonical_signature_from_selected_edges<3, 4>(graph, selected_edges, relabel);
		case 4: return build_canonical_signature_from_selected_edges<4, 6>(graph, selected_edges, relabel);
		case 5: return build_canonical_signature_from_selected_edges<5, 8>(graph, selected_edges, relabel);
		case 6: return build_canonical_signature_from_selected_edges<6, 10>(graph, selected_edges, relabel);
		case 7: return build_canonical_signature_from_selected_edges<7, 12>(graph, selected_edges, relabel);
		case 8: return build_canonical_signature_from_selected_edges<8, 14>(graph, selected_edges, relabel);
		case 9: return build_canonical_signature_from_selected_edges<9, 16>(graph, selected_edges, relabel);
		case 10: return build_canonical_signature_from_selected_edges<10, 18>(graph, selected_edges, relabel);
		case 11: return build_canonical_signature_from_selected_edges<11, 20>(graph, selected_edges, relabel);
		case 12: return build_canonical_signature_from_selected_edges<12, 22>(graph, selected_edges, relabel);
		default: return std::nullopt;
	}
}

template <Int N>
void add_v_catalog_entries(std::unordered_map<std::string, std::string>& catalog) {
	for (auto sequence : all_sequences((N - 5) / 2)) {
		auto graph = V_graph<N>(sequence);
		typename OddGraphdegZero<N + 1>::Basis be(graph);
		const auto canonical = OddGraphdegZero<N + 1>::canonized(be).getValue();
		catalog.emplace(
			signature_from_half_edges(OddGraphdegZero<N + 1>::N_VERTICES_, canonical.half_edges),
			"V_" + std::to_string(N) + "(" + sequence_string(sequence) + ")"
		);
	}
}

std::unordered_map<std::string, std::string> build_v_catalog() {
	std::unordered_map<std::string, std::string> catalog;
	add_v_catalog_entries<5>(catalog);
	add_v_catalog_entries<7>(catalog);
	add_v_catalog_entries<9>(catalog);
	add_v_catalog_entries<11>(catalog);
	return catalog;
}

std::unordered_set<std::string> admissible_subgraph_signatures(const FullGraph& graph) {
	const auto banks = coalgebra_utils::connected_triangle_banks(graph);
	std::unordered_set<std::string> signatures;
	std::vector<Int> selected_edges;

	const auto recurse = [&](auto&& self, std::size_t bank_index) -> void {
		if (bank_index == banks.size()) {
			const auto signature = admissible_signature(graph, selected_edges);
			if (signature.has_value()) {
				signatures.insert(*signature);
			}
			return;
		}

		self(self, bank_index + 1);

		const std::size_t old_size = selected_edges.size();
		selected_edges.insert(selected_edges.end(), banks[bank_index].begin(), banks[bank_index].end());
		self(self, bank_index + 1);
		selected_edges.resize(old_size);
	};

	recurse(recurse, 0);
	return signatures;
}

} // namespace

int main() {
	const auto v_catalog = build_v_catalog();

	const FullGC sum_1_then_2 = build_sum(1, 2);
	const ContGC d_1_then_2 = [&]() {
		auto tmp = sum_1_then_2;
		return tmp.d_contraction();
	}();

	auto primitive_opt = d_1_then_2.try_find_cont_primitive();
	if (!primitive_opt.has_value()) {
		std::cout << "primitive found: no\n";
		return 0;
	}

	auto primitive = *primitive_opt;
	primitive.standardize_all();
	primitive.sort_elements();

	std::size_t zero_admissible = 0;
	std::size_t only_v_admissible = 0;
	std::size_t has_non_v_admissible = 0;
	std::unordered_map<std::string, std::size_t> non_v_frequency;

	std::size_t i = 0;
	for (const auto& be : primitive.data()) {
		const auto& graph = be.getValue();
		const auto admissible = admissible_subgraph_signatures(graph);

		std::vector<std::string> matched_v;
		std::vector<std::string> non_v;
		for (const auto& signature : admissible) {
			const auto it = v_catalog.find(signature);
			if (it == v_catalog.end()) {
				non_v.push_back(signature);
				++non_v_frequency[signature];
			} else {
				matched_v.push_back(it->second);
			}
		}

		std::sort(matched_v.begin(), matched_v.end());
		matched_v.erase(std::unique(matched_v.begin(), matched_v.end()), matched_v.end());

		std::cout << "term " << i << ": admissible subgraphs = " << admissible.size()
			  << ", matched V-subgraphs = " << matched_v.size()
			  << ", non-V admissible subgraphs = " << non_v.size();
		if (!matched_v.empty()) {
			std::cout << " [";
			for (std::size_t j = 0; j < matched_v.size(); ++j) {
				if (j > 0) {
					std::cout << ", ";
				}
				std::cout << matched_v[j];
			}
			std::cout << "]";
		}
		std::cout << '\n';

		if (admissible.empty()) {
			++zero_admissible;
		} else if (non_v.empty()) {
			++only_v_admissible;
		} else {
			++has_non_v_admissible;
		}
		++i;
	}

	std::cout << "summary: primitive terms = " << primitive.data().size()
		  << ", zero admissible = " << zero_admissible
		  << ", only V-type admissible = " << only_v_admissible
		  << ", with non-V admissible subgraphs = " << has_non_v_admissible
		  << '\n';
	std::cout << "distinct non-V admissible signatures = " << non_v_frequency.size() << '\n';
	for (const auto& [signature, frequency] : non_v_frequency) {
		std::cout << "non-V signature seen in " << frequency << " primitive terms: " << signature << '\n';
	}

	return 0;
}
