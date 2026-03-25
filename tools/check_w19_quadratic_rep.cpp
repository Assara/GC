#include <iostream>
#include <vector>

#include "examplegraphs.hpp"

namespace {

std::vector<std::vector<bool>> all_sequences(int length) {
	std::vector<std::vector<bool>> result;
	const int count = 1 << length;
	result.reserve(count);

	for (int mask = 0; mask < count; ++mask) {
		std::vector<bool> seq;
		seq.reserve(length);
		for (int bit = 0; bit < length; ++bit) {
			seq.push_back(((mask >> bit) & 1) != 0);
		}
		result.push_back(seq);
	}

	return result;
}

template <typename GCType>
void normalize(GCType& gc) {
	gc.standardize_all();
	gc.sort_elements();
}

template <Int N>
typename OddGCdegZero<N + 1>::SplitGC build_u_total_all_sequences() {
	typename OddGCdegZero<N + 1>::SplitGC total;
	for (int k = 0; k <= (static_cast<int>(N) - 5) / 2; ++k) {
		for (auto seq : all_sequences(k)) {
			total += typename OddGCdegZero<N + 1>::SplitGC(U_graph<N>(seq));
		}
	}
	normalize(total);
	return total;
}

template <Int N>
OddGCdegZero<N + 1> build_v_empty() {
	std::vector<bool> empty;
	auto gc = OddGCdegZero<N + 1>(V_graph<N>(empty));
	normalize(gc);
	return gc;
}

template <typename GraphType>
bool is_three_four_valent(const GraphType& graph) {
	for (Int v : graph.valence_array()) {
		if (v != 3 && v != 4) {
			return false;
		}
	}
	return true;
}

} // namespace

int main() {
	using GC19 = OddGCdegZero<20>;

	auto u19 = build_u_total_all_sequences<19>();
	auto du19 = u19.d_contraction();
	normalize(du19);

	auto v_empty = build_v_empty<19>();
	auto lhs = v_empty;
	lhs += du19;
	normalize(lhs);

	std::size_t non_quadratic = 0;
	std::size_t quadratic = 0;
	std::size_t wheel_matches = 0;

	for (const auto& be : lhs.data()) {
		if (be.getValue() == v_empty.data().front().getValue()) {
			++wheel_matches;
		}
		if (is_three_four_valent(be.getValue())) {
			++quadratic;
		} else {
			++non_quadratic;
		}
	}

	std::cout << "d(U_19) support size = " << du19.size() << '\n';
	std::cout << "V_19(empty) + d(U_19) support size = " << lhs.size() << '\n';
	std::cout << "quadratic (3/4-valent) terms = " << quadratic << '\n';
	std::cout << "non-quadratic terms = " << non_quadratic << '\n';
	std::cout << "wheel matches remaining = " << wheel_matches << '\n';

	return 0;
}
