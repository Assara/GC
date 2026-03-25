#include <iostream>
#include <vector>

#include "examplegraphs.hpp"

namespace {

using FullGC = GC<13, 24, 0, 0, 0, 1>;
using ContGC = GC<12, 23, 0, 0, 0, 1>;

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

struct SumData {
	FullGC sum;
	std::size_t raw_graph_count = 0;
};

SumData build_sum(Int first_length, Int second_length) {
	SumData result;
	const auto first_sequences = all_sequences(first_length);
	const auto second_sequences = all_sequences(second_length);

	for (const auto& first : first_sequences) {
		for (const auto& second : second_sequences) {
			std::vector<std::vector<bool>> blocks{first, second};
			result.sum += FullGC(iterated_V_graph_lcr_extension<12>(blocks));
			++result.raw_graph_count;
		}
	}

	result.sum.standardize_all();
	result.sum.sort_elements();
	return result;
}

template <typename LinCombType>
std::size_t overlap_count(const LinCombType& left, const LinCombType& right) {
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

void print_summary(const char* label, const SumData& sum_data, const ContGC& differential) {
	std::cout
		<< label
		<< ": raw graphs = " << sum_data.raw_graph_count
		<< ", sum terms = " << sum_data.sum.data().size()
		<< ", d_contraction terms = " << differential.data().size()
		<< '\n';
}

void print_valence_summary(const char* label, const ContGC& differential) {
	std::size_t has_5_valent = 0;
	std::size_t has_6_valent = 0;
	std::size_t has_5_or_6_valent = 0;
	Int max_valence_seen = 0;

	for (const auto& be : differential.data()) {
		const auto valences = be.getValue().valence_array();
		bool term_has_5 = false;
		bool term_has_6 = false;
		Int term_max_valence = 0;

		for (const Int valence : valences) {
			if (valence == 5) {
				term_has_5 = true;
			}
			if (valence == 6) {
				term_has_6 = true;
			}
			if (valence > term_max_valence) {
				term_max_valence = valence;
			}
		}

		if (term_has_5) {
			++has_5_valent;
		}
		if (term_has_6) {
			++has_6_valent;
		}
		if (term_has_5 || term_has_6) {
			++has_5_or_6_valent;
		}
		if (term_max_valence > max_valence_seen) {
			max_valence_seen = term_max_valence;
		}
	}

	std::cout
		<< label
		<< ": terms with a 5-valent vertex = " << has_5_valent
		<< ", terms with a 6-valent vertex = " << has_6_valent
		<< ", terms with a 5- or 6-valent vertex = " << has_5_or_6_valent
		<< ", max valence seen = " << static_cast<int>(max_valence_seen)
		<< '\n';
}

} // namespace

int main() {
	const SumData sum_1_then_2 = build_sum(1, 2);
	const SumData sum_2_then_1 = build_sum(2, 1);

	const auto raw_d_1_then_2 = sum_1_then_2.sum.d_contraction_without_sort();
	const auto raw_d_2_then_1 = sum_2_then_1.sum.d_contraction_without_sort();

	const ContGC d_1_then_2 = [&]() {
		auto tmp = sum_1_then_2.sum;
		return tmp.d_contraction();
	}();

	const ContGC d_2_then_1 = [&]() {
		auto tmp = sum_2_then_1.sum;
		return tmp.d_contraction();
	}();

	print_summary("length-1 then length-2", sum_1_then_2, d_1_then_2);
	print_summary("length-2 then length-1", sum_2_then_1, d_2_then_1);

	std::cout << "raw d(length-1 then 2) terms: " << raw_d_1_then_2.size() << '\n';
	std::cout << "raw d(length-2 then 1) terms: " << raw_d_2_then_1.size() << '\n';
	print_valence_summary("d(length-1 then 2)", d_1_then_2);
	print_valence_summary("d(length-2 then 1)", d_2_then_1);
	std::cout << "input sums equal: " << ((sum_1_then_2.sum.data() == sum_2_then_1.sum.data()) ? "yes" : "no") << '\n';
	std::cout << "differentials equal: " << ((d_1_then_2.data() == d_2_then_1.data()) ? "yes" : "no") << '\n';
	std::cout << "common differential terms: " << overlap_count(d_1_then_2.data(), d_2_then_1.data()) << '\n';

	return 0;
}
