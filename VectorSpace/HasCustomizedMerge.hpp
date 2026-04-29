#pragma once

#include <concepts>

template<typename T, typename k>
class BasisElement;

namespace VectorSpace {

	template<typename T, typename k>
		concept HasCustomizedMerge =
		requires (BasisElement<T, k> a, BasisElement<T, k> b) {
			{ T::merge(a, b) } -> std::same_as<BasisElement<T, k>>;
		};

}
