#pragma once

#include <concepts>
#include "types.hpp"

namespace VectorSpace {

// Syntactic field interface; the field laws are the implementation's contract.
template <class F>
concept Field = std::regular<F> && std::constructible_from<F, int>
    && requires(F a, const F b, SmallSignedInt small) {
        { -b } -> std::same_as<F>;
        { a + b } -> std::same_as<F>;
        { a - b } -> std::same_as<F>;
        { a * b } -> std::same_as<F>;
        { a / b } -> std::same_as<F>;
        { a += b } -> std::same_as<F&>;
        { a -= b } -> std::same_as<F&>;
        { a *= b } -> std::same_as<F&>;
        { a /= b } -> std::same_as<F&>;
        { b.inv() } -> std::same_as<F>;
        { b * small } -> std::same_as<F>;
        { small * b } -> std::same_as<F>;
        { a *= small } -> std::same_as<F&>;
    };

} // namespace VectorSpace
