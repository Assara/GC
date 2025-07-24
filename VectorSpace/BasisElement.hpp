#pragma once

#include "../Types.hpp"
#include "HasCompare.hpp"
#include <memory>

// BasisElement class template storing a pointer to a value of type T and a coefficient
template<typename T>
class BasisElement {
private:
    std::unique_ptr<T> value;
    fieldType coefficient;

public:
    // Constructor from value
    BasisElement(const T& val, fieldType coeff = 1.0f)
        : value(std::make_unique<T>(val)), coefficient(coeff) {}

    // Constructor from unique_ptr
    BasisElement(std::unique_ptr<T>&& val, fieldType coeff)
        : value(std::move(val)), coefficient(coeff) {}

    // ✅ Deep copy constructor
    BasisElement(const BasisElement& other)
        : value(std::make_unique<T>(*other.value)), coefficient(other.coefficient) {}

    // ✅ Copy assignment operator
    BasisElement& operator=(const BasisElement& other) {
        if (this != &other) {
            value = std::make_unique<T>(*other.value);
            coefficient = other.coefficient;
        }
        return *this;
    }

    // ✅ Move constructor
    BasisElement(BasisElement&&) noexcept = default;
    BasisElement& operator=(BasisElement&&) noexcept = default;

    const T& getValue() const { return *value; }
    fieldType getCoefficient() const { return coefficient; }
    fieldType& getCoefficientRef() { return coefficient; }

    signedInt compare(const BasisElement& other) const {
        return value->compare(*other.value);
    }
};

