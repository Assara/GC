#pragma once
#include <string>
#include "Z32783.hpp"

// Change this alias to change the integer type project-wide
using Int = std::uint8_t;
using signedInt = int;

using bigInt = std::vector<Int>::size_type;
using fieldType = Z32783;


template<typename T>
struct TypeName {
    static std::string name() { return "unknown"; }
};

template<>
struct TypeName<double> {
    static std::string name() { return "double"; }
};

#include "Z32783.hpp"

template<>
struct TypeName<Z32783> {
    static std::string name() { return "Z32783"; }
};
