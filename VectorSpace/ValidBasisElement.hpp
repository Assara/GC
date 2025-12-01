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
        { T::std(b) } -> std::same_as<void>;
    };

    
} // namespace VectorSpace