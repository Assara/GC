#include <iostream>

#include "examplegraphs.hpp"

int main() {
	using GraphType = OddGraphdegZero<16>;
	using Standardizer = GraphStandardizer<
		GraphType::N_VERTICES_,
		GraphType::N_EDGES_,
		GraphType::N_OUT_HAIR_,
		GraphType::N_IN_HAIR_,
		GraphType::C_,
		GraphType::D_,
		fieldType
	>;
	using Basis = BasisElement<GraphType, fieldType>;

	Standardizer standardizer;

	GraphType original_graph = wheel_graph<15>();
	Basis original(original_graph, fieldType{1});

	GraphType standard_graph = original_graph;
	Basis regular = standardizer.standardize(standard_graph, fieldType{1});

	Basis standard_copy(regular.getValue(), regular.getCoefficient());
	Basis via_original = standardizer.standardize2(original);
	Basis via_standard_copy = standardizer.standardize2(standard_copy);

	const bool same_graphs = via_original.getValue() == via_standard_copy.getValue();
	const bool same_coefficients = via_original.getCoefficient() == via_standard_copy.getCoefficient();

	std::cout << "standardize2(original) == standardize2(standardized copy): "
	          << (same_graphs && same_coefficients ? "yes" : "no") << '\n';
	std::cout << "same graphs: " << (same_graphs ? "yes" : "no") << '\n';
	std::cout << "same coefficients: " << (same_coefficients ? "yes" : "no") << '\n';
	std::cout << "regular coefficient: " << regular.getCoefficient() << '\n';
	std::cout << "standardize2(original) coefficient: " << via_original.getCoefficient() << '\n';
	std::cout << "standardize2(standardized copy) coefficient: "
	          << via_standard_copy.getCoefficient() << '\n';

	return same_graphs && same_coefficients ? 0 : 1;
}
