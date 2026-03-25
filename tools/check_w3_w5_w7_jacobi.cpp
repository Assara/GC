#include <iostream>
#include <vector>

#include "coalgebra_utils.hpp"
#include "examplegraphs.hpp"

int main() {
	std::vector<bool> L{false};
	std::vector<bool> LL{false, false};
	std::vector<bool> LR{false, true};

	// W_3 is already 3-valent, so we use the wheel graph itself.
	auto w3 = OddGCdegZero<4>(wheel_graph<3>());

	// Use the same low-valence representatives as in the current paper/examples.
	auto w5 = OddGCdegZero<6>(V_graph<5>(L));

	auto w7 = OddGCdegZero<8>(V_graph<7>(LL));
	w7 += OddGCdegZero<8>(V_graph<7>(LR));

	std::cout << "computing first-level brackets..." << std::endl;
	const auto bracket_5_7 = coalgebra_utils::gra_lie(w5, w7);
	const auto bracket_7_3 = coalgebra_utils::gra_lie(w7, w3);
	const auto bracket_3_5 = coalgebra_utils::gra_lie(w3, w5);
	std::cout << "[W5,W7] terms: " << bracket_5_7.data().size() << std::endl;
	std::cout << "[W7,W3] terms: " << bracket_7_3.data().size() << std::endl;
	std::cout << "[W3,W5] terms: " << bracket_3_5.data().size() << std::endl;

	std::cout << "computing nested bracket 1..." << std::endl;
	const auto jacobi_1 = coalgebra_utils::gra_lie(w3, bracket_5_7);
	std::cout << "[W3,[W5,W7]] terms: " << jacobi_1.data().size() << std::endl;
	std::cout << "computing nested bracket 2..." << std::endl;
	const auto jacobi_2 = coalgebra_utils::gra_lie(w5, bracket_7_3);
	std::cout << "[W5,[W7,W3]] terms: " << jacobi_2.data().size() << std::endl;
	std::cout << "computing nested bracket 3..." << std::endl;
	const auto jacobi_3 = coalgebra_utils::gra_lie(w7, bracket_3_5);
	std::cout << "[W7,[W3,W5]] terms: " << jacobi_3.data().size() << std::endl;

	std::cout << "summing Jacobi terms..." << std::endl;
	auto jacobi_sum = jacobi_1;
	jacobi_sum += jacobi_2;
	jacobi_sum += jacobi_3;
	jacobi_sum.standardize_all();
	jacobi_sum.sort_elements();

	std::cout << "Jacobi sum terms: " << jacobi_sum.data().size() << std::endl;
	std::cout << "Jacobi identity holds: " << (jacobi_sum.data().size() == 0 ? "yes" : "no") << std::endl;

	return 0;
}
