#pragma once
#include "../Types.hpp"
#include "HasCompare.hpp"

// BasisElement class template storing a value of type T and a coefficient
template<typename T>
class BasisElement {
private:
    T value;
    fieldType coefficient;
public:
    // Constructor to initialize value and coefficient (default coefficient = 1.0)
    BasisElement(const T& val, fieldType coeff = 1.0f)
        : value(val), coefficient(coeff) {}

    // Accessors for the value and coefficient
    const T& getValue() const { return value; }
    fieldType getCoefficient() const { return coefficient; }
    fieldType& getCoefficientRef() { return coefficient; }

    // Compare this BasisElement with another by comparing their values using T.compare
    signedInt compare(const BasisElement& other) const {
        return value.compare(other.value);
    }

};
