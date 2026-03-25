#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "DegreeZeroGraphWedge.hpp"
#include "examplegraphs.hpp"

namespace {

using FullGC = GC<13, 24, 0, 0, 0, 1>;
using ContGC = GC<12, 23, 0, 0, 0, 1>;
using OneWedge = DegreeZeroGraphWedge<48, 1>;
using TwoWedge = DegreeZeroGraphWedge<48, 2>;
using TwoWedgeLinComb = VectorSpace::LinComb<TwoWedge, fieldType>;

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

std::string sequence_string(const std::vector<bool>& sequence) {
	std::string result;
	for (bool bit : sequence) {
		result.push_back(bit ? 'R' : 'L');
	}
	return result.empty() ? "empty" : result;
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

template <Int N_VERTICES, Int N_EDGES>
std::string factor_signature_from_wedge(const TwoWedge& wedge, std::size_t factor_index) {
	using GraphType = Graph<N_VERTICES, N_EDGES, 0, 0, 0, 1, fieldType>;
	GraphType graph;
	const std::size_t offset = wedge.factor_offset(factor_index);
	for (std::size_t i = 0; i < static_cast<std::size_t>(2 * N_EDGES); ++i) {
		graph.half_edges[i] = wedge.data()[offset + i];
	}
	typename GraphType::Basis be(graph);
	const auto canonical = GraphType::canonized(be).getValue();
	return signature_from_half_edges(N_VERTICES, canonical.half_edges);
}

std::pair<std::string, std::string> classify_factor(const TwoWedge& wedge, std::size_t factor_index,
	const std::unordered_map<std::string, std::string>& v_catalog) {
	const Int count = wedge.factor_half_edge_count(factor_index);
	std::string signature;
	switch (count) {
		case 20: signature = factor_signature_from_wedge<6, 10>(wedge, factor_index); break;
		case 28: signature = factor_signature_from_wedge<8, 14>(wedge, factor_index); break;
		case 36: signature = factor_signature_from_wedge<10, 18>(wedge, factor_index); break;
		case 44: signature = factor_signature_from_wedge<12, 22>(wedge, factor_index); break;
		default: return {"?", ""};
	}

	const auto it = v_catalog.find(signature);
	return it == v_catalog.end() ? std::pair{"?", signature} : std::pair{it->second, signature};
}

TwoWedgeLinComb cobracket_of_gc(const FullGC& gc) {
	TwoWedgeLinComb result;

	for (const auto& be : gc.data()) {
		const auto wedge = OneWedge::from_graph(be.getValue());
		auto term = wedge.cobracket();
		term.scalar_multiply(be.getCoefficient());
		result += term;
	}

	result.standardize_and_sort();
	return result;
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

	const auto cobracket = cobracket_of_gc(primitive);

	std::map<std::pair<int, int>, std::size_t> split_shape_counts;
	for (const auto& be : cobracket) {
		const auto& wedge = be.getValue();
		split_shape_counts[{wedge.factor_half_edge_count(0), wedge.factor_half_edge_count(1)}]++;
	}

	std::cout << "primitive terms: " << primitive.data().size() << '\n';
	std::cout << "cobracket support size: " << cobracket.size() << '\n';
	for (const auto& [shape, count] : split_shape_counts) {
		std::cout << "split shape (" << shape.first << ", " << shape.second << "): " << count << " terms\n";
	}
	for (const auto& be : cobracket) {
		const auto& wedge = be.getValue();
		const auto left = classify_factor(wedge, 0, v_catalog);
		const auto right = classify_factor(wedge, 1, v_catalog);
		std::cout << "term coeff = " << be.getCoefficient()
			  << ", factors = "
			  << left.first
			  << " ^ "
			  << right.first
			  << '\n';
		if (left.first == "?") {
			std::cout << "  left signature: " << left.second << '\n';
		}
		if (right.first == "?") {
			std::cout << "  right signature: " << right.second << '\n';
		}
	}

	return 0;
}
