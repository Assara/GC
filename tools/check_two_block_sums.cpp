#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "examplegraphs.hpp"

namespace {

using FullGraph = Graph<15, 28, 0, 0, 0, 1, fieldType>;
using FullGC = GC<15, 28, 0, 0, 0, 1>;
using ContGC = GC<14, 27, 0, 0, 0, 1>;
using DiffLinComb = typename ContGC::L;
using FullLinComb = typename FullGC::L;

void initialize_seed(FullGraph& graph, Int& edge_it) {
	graph.setEdge(edge_it++, 0, 1);
	graph.setEdge(edge_it++, 0, 2);
	graph.setEdge(edge_it++, 0, 3);
	graph.setEdge(edge_it++, 1, 2);
	graph.setEdge(edge_it++, 2, 3);
}

void append_arc_vertices(
	FullGraph& graph,
	Int& edge_it,
	Int& next_vertex,
	Int center,
	Int right,
	Int left,
	Int arc_vertex_count
) {
	Int arc_pos = left;
	for (Int i = 0; i < arc_vertex_count; ++i) {
		graph.setEdge(edge_it++, arc_pos, next_vertex);
		graph.setEdge(edge_it++, center, next_vertex);
		arc_pos = next_vertex;
		++next_vertex;
	}
	graph.setEdge(edge_it++, arc_pos, right);
}

FullGraph build_two_block_graph(
	const std::vector<bool>& first_sequence,
	const std::vector<bool>& second_sequence,
	const std::array<Int, 4>& attachment_pattern
) {
	FullGraph graph;
	Int edge_it = 0;
	initialize_seed(graph, edge_it);

	Int left = 1;
	Int center = 2;
	Int right = 3;
	Int next_vertex = 4;

	append_V_sequence_block(graph, edge_it, next_vertex, left, center, right, first_sequence);
	graph.setEdge(edge_it++, left, right);

	extend_V_seed_block_general(
		graph,
		edge_it,
		next_vertex,
		attachment_pattern[0],
		attachment_pattern[1],
		attachment_pattern[2],
		attachment_pattern[3],
		left,
		center,
		right
	);

	append_V_sequence_block(graph, edge_it, next_vertex, left, center, right, second_sequence);
	append_arc_vertices(graph, edge_it, next_vertex, center, right, left, 2);

	if (next_vertex != Int{15} || edge_it != Int{28}) {
		std::cerr << "unexpected graph size while building a two-block graph\n";
		std::exit(EXIT_FAILURE);
	}

	return graph;
}

std::vector<std::vector<bool>> all_sequences(const Int length, const bool require_left_start) {
	std::vector<std::vector<bool>> result;
	if (length == 0 || (length == 1 && require_left_start)) {
		result.push_back({});
		if (length == 1 && require_left_start) {
			result.front().push_back(false);
		}
		return result;
	}

	if (require_left_start) {
		const std::size_t count = std::size_t{1} << (length - 1);
		result.reserve(count);

		for (std::size_t mask = 0; mask < count; ++mask) {
			std::vector<bool> sequence;
			sequence.reserve(length);
			sequence.push_back(false);
			for (Int i = 1; i < length; ++i) {
				sequence.push_back(((mask >> (i - 1)) & 1U) != 0U);
			}
			result.push_back(std::move(sequence));
		}
		return result;
	}

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

struct SumData {
	FullGC sum;
	std::size_t raw_graph_count = 0;
};

SumData build_sum(const Int first_length, const Int second_length, const bool require_left_start) {
	SumData result;
	const auto first_sequences = all_sequences(first_length, require_left_start);
	const auto second_sequences = all_sequences(second_length, require_left_start);

	for (const auto& first : first_sequences) {
		FullGraph first_block_graph;
		Int edge_it = 0;
		initialize_seed(first_block_graph, edge_it);

		Int left = 1;
		Int center = 2;
		Int right = 3;
		Int next_vertex = 4;
		append_V_sequence_block(first_block_graph, edge_it, next_vertex, left, center, right, first);
		first_block_graph.setEdge(edge_it++, left, right);

		const auto attachment_patterns = all_V_seed_attachment_patterns(0, left, center, right);

		for (const auto& second : second_sequences) {
			for (const auto& pattern : attachment_patterns) {
				result.sum += FullGC(build_two_block_graph(first, second, pattern));
				++result.raw_graph_count;
			}
		}
	}

	result.sum.standardize_all();
	return result;
}

FullLinComb difference(const FullGC& left, const FullGC& right) {
	FullLinComb diff = left.data();
	diff = diff.add_scaled(right.data(), fieldType{-1});
	diff.standardize_all();
	diff.sort_elements();
	return diff;
}

FullLinComb doubled(const FullGC& gc) {
	FullLinComb result = gc.data();
	result = result.add_scaled(gc.data(), fieldType{1});
	result.standardize_all();
	result.sort_elements();
	return result;
}

FullLinComb scaled_copy(const FullGC& gc, const fieldType& scalar) {
	FullLinComb result;
	result = result.add_scaled(gc.data(), scalar);
	result.standardize_all();
	result.sort_elements();
	return result;
}

DiffLinComb difference(const ContGC& left, const ContGC& right) {
	DiffLinComb diff = left.data();
	diff = diff.add_scaled(right.data(), fieldType{-1});
	diff.standardize_all();
	diff.sort_elements();
	return diff;
}

std::size_t overlap_count(const DiffLinComb& left, const DiffLinComb& right) {
	const auto& a = left.raw_elements();
	const auto& b = right.raw_elements();
	std::size_t i = 0;
	std::size_t j = 0;
	std::size_t overlap = 0;

	while (i < a.size() && j < b.size()) {
		const auto cmp = a[i].getValue().compare(b[j].getValue());
		if (cmp < 0) {
			++i;
		} else if (cmp > 0) {
			++j;
		} else {
			++overlap;
			++i;
			++j;
		}
	}

	return overlap;
}

void print_summary(
	const char* label,
	const SumData& sum_data,
	const ContGC& differential
) {
	std::cout
		<< label
		<< ": raw graphs = " << sum_data.raw_graph_count
		<< ", sum has " << sum_data.sum.data().size() << " term(s), "
		<< "d_contraction has " << differential.data().size() << " term(s)\n";
}

} // namespace

int main() {
	const SumData sum_1_then_2_left = build_sum(1, 2, true);
	const SumData sum_2_then_1_left = build_sum(2, 1, true);
	const SumData sum_1_then_2_all = build_sum(1, 2, false);
	const SumData sum_2_then_1_all = build_sum(2, 1, false);

	const ContGC d_1_then_2 = [&]() {
		FullGC tmp = sum_1_then_2_left.sum;
		return tmp.d_contraction();
	}();

	const ContGC d_2_then_1 = [&]() {
		FullGC tmp = sum_2_then_1_left.sum;
		return tmp.d_contraction();
	}();

	print_summary("length-1 first, length-2 second", sum_1_then_2_left, d_1_then_2);
	print_summary("length-2 first, length-1 second", sum_2_then_1_left, d_2_then_1);
	std::cout << "input sums equal: " << ((sum_1_then_2_left.sum.data() == sum_2_then_1_left.sum.data()) ? "yes" : "no") << '\n';

	const bool equal = (d_1_then_2.data() == d_2_then_1.data());
	std::cout << "differentials equal: " << (equal ? "yes" : "no") << '\n';
	std::cout << "common differential terms: " << overlap_count(d_1_then_2.data(), d_2_then_1.data()) << '\n';

	const DiffLinComb diff = difference(d_1_then_2, d_2_then_1);
	std::cout << "difference has " << diff.size() << " term(s)\n";

	if (diff.size() != 0U && diff.size() <= 10U) {
		diff.print();
	}

	const FullLinComb doubled_1_then_2_left = doubled(sum_1_then_2_left.sum);
	const FullLinComb doubled_2_then_1_left = doubled(sum_2_then_1_left.sum);
	const FullLinComb quadrupled_1_then_2_left = scaled_copy(sum_1_then_2_left.sum, fieldType{4});
	const FullLinComb quadrupled_2_then_1_left = scaled_copy(sum_2_then_1_left.sum, fieldType{4});
	const FullLinComb diff_1_then_2_double = [&]() {
		FullLinComb diff = sum_1_then_2_all.sum.data();
		diff = diff.add_scaled(doubled_1_then_2_left, fieldType{-1});
		diff.standardize_all();
		diff.sort_elements();
		return diff;
	}();
	const FullLinComb diff_2_then_1_double = [&]() {
		FullLinComb diff = sum_2_then_1_all.sum.data();
		diff = diff.add_scaled(doubled_2_then_1_left, fieldType{-1});
		diff.standardize_all();
		diff.sort_elements();
		return diff;
	}();
	const FullLinComb diff_1_then_2_quadruple = [&]() {
		FullLinComb diff = sum_1_then_2_all.sum.data();
		diff = diff.add_scaled(quadrupled_1_then_2_left, fieldType{-1});
		diff.standardize_all();
		diff.sort_elements();
		return diff;
	}();
	const FullLinComb diff_2_then_1_quadruple = [&]() {
		FullLinComb diff = sum_2_then_1_all.sum.data();
		diff = diff.add_scaled(quadrupled_2_then_1_left, fieldType{-1});
		diff.standardize_all();
		diff.sort_elements();
		return diff;
	}();

	std::cout
		<< "all-sequence sum, length-1 then 2: raw graphs = " << sum_1_then_2_all.raw_graph_count
		<< ", standardised terms = " << sum_1_then_2_all.sum.data().size()
		<< ", equals 0: " << (sum_1_then_2_all.sum.data().size() == 0U ? "yes" : "no")
		<< ", equals 2x left-starting: " << (diff_1_then_2_double.size() == 0U ? "yes" : "no")
		<< ", equals 4x left-starting: " << (diff_1_then_2_quadruple.size() == 0U ? "yes" : "no")
		<< '\n';
	std::cout
		<< "all-sequence sum, length-2 then 1: raw graphs = " << sum_2_then_1_all.raw_graph_count
		<< ", standardised terms = " << sum_2_then_1_all.sum.data().size()
		<< ", equals 0: " << (sum_2_then_1_all.sum.data().size() == 0U ? "yes" : "no")
		<< ", equals 2x left-starting: " << (diff_2_then_1_double.size() == 0U ? "yes" : "no")
		<< ", equals 4x left-starting: " << (diff_2_then_1_quadruple.size() == 0U ? "yes" : "no")
		<< '\n';

	return EXIT_SUCCESS;
}
