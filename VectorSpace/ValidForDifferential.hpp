#pragma once

#include "BasisElement.hpp"
#include "HasCompare.hpp"
#include <concepts>
#include <vector>

namespace VectorSpace {

template<typename T>
concept ValidBasisElement =
    HasCompare<T> &&
    requires(BasisElement<T>& b) {
        { T::std(b) } -> std::same_as<void>;
    };
} // namespace VectorSpace