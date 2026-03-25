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

void check_case(const char* label, const ContGC& target) {
	std::cout << "==== " << label << " ====\n";
	std::cout << "target terms: " << target.data().size() << '\n';

	auto primitive_opt = target.try_find_cont_primitive();
	if (!primitive_opt.has_value()) {
		std::cout << "primitive found: no\n";
		return;
	}

	auto primitive = *primitive_opt;
	primitive.standardize_all();
	primitive.sort_elements();

	auto d_primitive = primitive.d_contraction();
	d_primitive.standardize_all();
	d_primitive.sort_elements();

	std::cout << "primitive found: yes\n";
	std::cout << "primitive terms: " << primitive.data().size() << '\n';
	std::cout << "d(primitive) terms: " << d_primitive.data().size() << '\n';
	std::cout << "d(primitive) == target: " << ((d_primitive.data() == target.data()) ? "yes" : "no") << '\n';
}

} // namespace

int main() {
	const FullGC sum_1_then_2 = build_sum(1, 2);
	const FullGC sum_2_then_1 = build_sum(2, 1);

	const ContGC d_1_then_2 = [&]() {
		auto tmp = sum_1_then_2;
		return tmp.d_contraction();
	}();

	const ContGC d_2_then_1 = [&]() {
		auto tmp = sum_2_then_1;
		return tmp.d_contraction();
	}();

	check_case("length-1 then length-2", d_1_then_2);
	check_case("length-2 then length-1", d_2_then_1);

	return 0;
}
