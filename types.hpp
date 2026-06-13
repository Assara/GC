#pragma once
#include <cstdint>
#include <vector>
#include "Q.hpp"
#include "Z2179564669.hpp"
#include "Z34821139123.hpp"
#include "Z32783.hpp"

// Change this alias to change the integer type project-wide
using Int = std::uint8_t;
using signedInt = int;
using bigInt = std::vector<Int>::size_type;
//using fieldType = Z2179564669;
using fieldType = Z34821139123;
//using fieldType = Z32783;
//using fieldType = Q;
