#include <iostream>
#include <vector>

#include "examplegraphs.hpp"

namespace {

std::vector<std::vector<bool>> sequences_starting_left(int length) {
	if (length == 0) {
		return {{}};
	}

	std::vector<std::vector<bool>> result;
	const int free_bits = length - 1;
	const int count = 1 << free_bits;
	result.reserve(count);

	for (int mask = 0; mask < count; ++mask) {
		std::vector<bool> seq;
		seq.reserve(length);
		seq.push_back(false);
		for (int bit = 0; bit < free_bits; ++bit) {
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
typename OddGCdegZero<N + 1>::SplitGC build_u_total() {
	typename OddGCdegZero<N + 1>::SplitGC total;
	for (int k = 0; k <= (static_cast<int>(N) - 5) / 2; ++k) {
		for (auto seq : sequences_starting_left(k)) {
			total += typename OddGCdegZero<N + 1>::SplitGC(U_graph<N>(seq));
		}
	}
	normalize(total);
	return total;
}

template <typename GCType>
bool same_normalized(GCType a, GCType b) {
	normalize(a);
	normalize(b);
	return a.data() == b.data();
}

} // namespace

int main() {
	auto u13 = build_u_total<13>();
	auto du13 = u13.d_contraction();
	normalize(du13);
	std::vector<bool> L{false};
	std::vector<bool> R{true};
	auto v13L = OddGCdegZero<14>(V_graph<13>(L));
	auto v13R = OddGCdegZero<14>(V_graph<13>(R));
	normalize(v13L);
	normalize(v13R);
	std::cout << "size " << du13.size() << '\n';
	int shown = 0;
	for (const auto& be : du13.data()) {
		OddGCdegZero<14> term(be);
		normalize(term);
		std::cout << "coeff " << be.getCoefficient() << '\n';
		for (Int e = 0; e < decltype(du13)::GraphType::N_EDGES_; ++e) {
			auto edge = be.getValue().getEdge(e);
			std::cout << "(" << static_cast<int>(edge.first) << "," << static_cast<int>(edge.second) << ")";
			if (e + 1 < decltype(du13)::GraphType::N_EDGES_) {
				std::cout << ' ';
			}
		}
		std::cout << '\n';
		std::cout << "matches V13(L): " << (term.data() == v13L.data()) << '\n';
		std::cout << "matches V13(R): " << (term.data() == v13R.data()) << '\n';
		++shown;
		if (shown == 2) {
			break;
		}
	}
	return 0;
}
