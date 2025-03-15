#pragma once
#include <vector>
#include <utility>
#include "Types.hpp"
#include <concepts>

// Define a concept that requires a type T to have a standardize(G, fieldType)
// method that returns something convertible to std::pair<G, fieldType>.
template<typename Standardizer, typename G, typename FieldType>
concept HasStandardize = requires(Standardizer s, G g, FieldType k) {
    { s.standardize(g, k) } -> std::convertible_to<std::pair<G, FieldType>>;
};

template<typename G, typename Standardizer, typename Coefficient = fieldType>
requires HasStandardize<Standardizer, G, fieldType>
class VectorSpace {
public:
    // Each term is a pair: (object of type G, coefficient).
    std::vector<std::pair<G, Coefficient>> terms;

    VectorSpace() = default;

    // Adds a term to the linear combination.
    // If g is already present (as determined by operator== on G),
    // its coefficient is updated. If the coefficient becomes zero, the term is removed.
    void addTerm(const G &g, Coefficient c) {
        for (auto it = terms.begin(); it != terms.end(); ++it) {
            if (it->first == g) { // Requires that G has an operator==
                it->second += c;
                if (it->second == 0)
                    terms.erase(it);
                return;
            }
        }
        if (c != 0)
            terms.push_back({g, c});
    }

    // Vector addition.
    VectorSpace operator+(const VectorSpace &other) const {
        VectorSpace result = *this;
        for (const auto &term : other.terms)
            result.addTerm(term.first, term.second);
        return result;
    }

    // Scalar multiplication.
    VectorSpace operator*(Coefficient scalar) const {
        VectorSpace result;
        for (const auto &term : terms) {
            Coefficient newCoeff = term.second * scalar;
            if (newCoeff != 0)
                result.terms.push_back({term.first, newCoeff});
        }
        return result;
    }

    // Standardizes an object using the provided Standardizer.
    // This method creates a Standardizer instance and calls its standardize() method.
    std::pair<G, fieldType> standardize(const G &g, fieldType k) const {
        Standardizer s;
        return s.standardize(g, k);
    }
};