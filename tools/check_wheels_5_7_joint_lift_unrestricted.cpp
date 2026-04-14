#include <iostream>
#include <unordered_set>

#include "coalgebra_search.hpp"
#include "coalgebra_utils.hpp"
#include "examplegraphs.hpp"

namespace {

using W5GC = OddGCdegZero<6>;
using W7GC = OddGCdegZero<8>;
using SearchGC = coalgebra_utils::GraInsertLowValentGC<W5GC, W7GC>;
using GraphType = SearchGC::GraphType;

W5GC build_w5() {
	return W5GC(wheel_graph<5>());
}

W7GC build_w7() {
	return W7GC(wheel_graph<7>());
}

std::unordered_set<GraphType> build_unrestricted_lie_support(const W5GC& w5, const W7GC& w7) {
	const auto lie = coalgebra_utils::gra_lie_unrestricted(w5, w7);

	std::unordered_set<GraphType> support;
	support.reserve(lie.data().size());
	for (const auto& be : lie.data()) {
		support.insert(be.getValue());
	}
	return support;
}

} // namespace

int main() {
	const auto w5 = build_w5();
	const auto w7 = build_w7();
	const auto unrestricted_lie = coalgebra_utils::gra_lie_unrestricted(w5, w7);
	const auto support = build_unrestricted_lie_support(w5, w7);

	std::cout << "unrestricted lie size: " << unrestricted_lie.data().size() << std::endl;
	std::cout << "unrestricted lie support size: " << support.size() << std::endl;

	auto cycle_opt = coalgebra_search::find_cobracket_lift_in_support<SearchGC>(
		support,
		w5,
		w7,
		true
	);

	if (cycle_opt.has_value() == false) {
		std::cout << "no cycle with target cobracket W5 ^ W7 found in unrestricted wheel lie support" << std::endl;
		return 0;
	}

	const auto& cycle = *cycle_opt;
	std::cout << "found cycle with target cobracket W5 ^ W7" << std::endl;
	std::cout << "cycle support size: " << cycle.data().size() << std::endl;

	return 0;
}
