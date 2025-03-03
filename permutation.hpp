#pragma once
#include <array>
#include "types.hpp"


template<Int n>
class Permutation {
public:
    std::array<Int, n> p;

    // Default constructor: creates the identity permutation.
    Permutation() {
        for (Int i = 0; i < n; ++i) {
            p[i] = i;
        }
    }

    // Constructor that initializes from an existing array.
    explicit Permutation(const std::array<Int, n>& arr)
        : p(arr)
    {}

    const Int& operator[](Int i) const {
        return p[i];
    }

    // Check whether the permutation is valid (i.e., a rearrangement of 0, 1, ..., n-1).
    bool is_valid() const {
        std::array<bool, n> seen{};
        for (Int i = 0; i < n; ++i) {
            if (p[i] < 0 || p[i] >= n)
                return false;
            if (seen[p[i]])
                return false;
            seen[p[i]] = true;
        }
        return true;
    }

    // Compute the inverse permutation.
    Permutation inverse() const {
        Permutation inv;
        for (Int i = 0; i < n; ++i) {
            inv.p[p[i]] = i;
        }
        return inv;
    }

    // Compose this permutation with another.
    // Returns a permutation q such that: q(i) = p(other.p(i)).
    Permutation compose(const Permutation& other) const {
        Permutation result;
        for (Int i = 0; i < n; ++i) {
            result.p[i] = p[other.p[i]];
        }
        return result;
    }


    signedInt sign() const {
        std::array<bool, n> visited{};  // n is our template parameter for size.
        signedInt s = 1;
        for (Int i = 0; i < n; ++i) {
            if (!visited[i]) {
                Int cycle_length = 0;
                Int j = i;
                while (!visited[j]) {
                    visited[j] = true;
                    ++cycle_length;
                    j = p[j];
                }
                // A cycle of length L has sign (-1)^(L-1)
                if ((cycle_length - 1) % 2 != 0) {
                    s = -s;
                }
            }
        }
        return s;
    }

};