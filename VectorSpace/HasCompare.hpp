
#pragma once
#include "../Types.hpp"
#include <concepts>

// ✅ Any type with a compare(const T&) -> signedInt

template<typename T>
concept HasCompare = requires(const T& a, const T& b) {
    { a.compare(b) } -> std::same_as<signedInt>;
    { a < b } -> std::convertible_to<bool>;
    { a == b } -> std::convertible_to<bool>;
};