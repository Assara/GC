#pragma once

#include "BasisElement.hpp"
#include "HasCompare.hpp"
#include <concepts>
#include <vector>

namespace VectorSpace {

template<typename T, typename k>
concept ValidBasisElement =
    //HasCompare<T> &&
    
    
    requires(BasisElement<T, k>& b) {
        { T::standardized(b) } -> std::same_as<BasisElement<T, k>>;
    };

    
} // namespace VectorSpace