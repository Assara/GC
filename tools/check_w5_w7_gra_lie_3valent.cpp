#include <iostream>
#include <vector>

#include "coalgebra_utils.hpp"
#include "examplegraphs.hpp"

namespace {

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

} // namespace

int main() {
	std::vector<bool> L{false};
	std::vector<bool> LL{false, false};
	std::vector<bool> LR{false, true};

	auto w5 = OddGCdegZero<6>(V_graph<5>(L));

	auto w7 = OddGCdegZero<8>(V_graph<7>(LL));
	w7 += OddGCdegZero<8>(V_graph<7>(LR));

	const auto pre_5_7 = coalgebra_utils::gra_pre_lie_3valent_only(w5, w7);
	const auto pre_7_5 = coalgebra_utils::gra_pre_lie_3valent_only(w7, w5);
	const auto lie = coalgebra_utils::gra_lie_3valent_only(w5, w7);
	const auto raw_d_pre_5_7 = pre_5_7.d_contraction_without_sort();
	const auto raw_d_pre_7_5 = pre_7_5.d_contraction_without_sort();
	const auto d_pre_5_7 = [&]() {
		auto tmp = pre_5_7;
		return tmp.d_contraction();
	}();
	const auto d_pre_7_5 = [&]() {
		auto tmp = pre_7_5;
		return tmp.d_contraction();
	}();
	const auto raw_d_lie = lie.d_contraction_without_sort();
	const auto d_lie = [&]() {
		auto tmp = lie;
		return tmp.d_contraction();
	}();

	std::cout << "3-valent-only pre_5_7 terms: " << pre_5_7.data().size() << '\n';
	std::cout << "3-valent-only pre_7_5 terms: " << pre_7_5.data().size() << '\n';
	std::cout << "3-valent-only lie terms: " << lie.data().size() << '\n';
	std::cout << "3-valent-only raw d(pre_5_7) terms: " << raw_d_pre_5_7.size() << '\n';
	std::cout << "3-valent-only raw d(pre_7_5) terms: " << raw_d_pre_7_5.size() << '\n';
	std::cout << "3-valent-only standardised d(pre_5_7) terms: " << d_pre_5_7.data().size() << '\n';
	std::cout << "3-valent-only standardised d(pre_7_5) terms: " << d_pre_7_5.data().size() << '\n';
	std::cout << "3-valent-only common standardised terms in the two pre-Lie differentials: "
		  << overlap_count(d_pre_5_7.data(), d_pre_7_5.data()) << '\n';
	std::cout << "3-valent-only raw d(lie) terms: " << raw_d_lie.size() << '\n';
	std::cout << "3-valent-only standardised d(lie) terms: " << d_lie.data().size() << '\n';
	std::cout << "3-valent-only d(lie) is zero: " << (d_lie.data().size() == 0 ? "yes" : "no") << '\n';

	return 0;
}
