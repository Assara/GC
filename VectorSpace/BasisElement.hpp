#pragma once

#include "../Types.hpp"
#include "HasCompare.hpp"
#include <memory>

template<typename T, typename k>
class BasisElement {
private:
    std::unique_ptr<T> value;
    k coefficient;

public:
    // Construct from value / unique_ptr
    explicit BasisElement(const T& val, k coeff = k{1})
        : value(std::make_unique<T>(val)), coefficient(coeff) {}

    BasisElement(std::unique_ptr<T>&& val, k coeff)
        : value(std::move(val)), coefficient(coeff) {}

    // Deep copy
    BasisElement(const BasisElement& other)
        : value(std::make_unique<T>(*other.value)), coefficient(other.coefficient) {}

    BasisElement& operator=(const BasisElement& other) {
        if (this != &other) {
            value = std::make_unique<T>(*other.value);
            coefficient = other.coefficient;
        }
        return *this;
    }

    // Moves OK
    BasisElement(BasisElement&&) noexcept = default;
    BasisElement& operator=(BasisElement&&) noexcept = default;

    // Mutators
    void multiplyCoefficient(k factor) noexcept { coefficient *= factor; }

    // Accessors (const-correct; note order: `const` then `noexcept`)
    const T* borrowValue() const noexcept { return value.get(); }
    T*       borrowValue()       noexcept { return value.get(); }

    const T& getValue() const noexcept { return *value; }
    T&       getValue()       noexcept { return *value; }

    k  getCoefficient() const noexcept { return coefficient; }
    k& getCoefficientRef()     noexcept { return coefficient; }

    // If you really need this, expose raw ptrs, not the unique_ptr itself:
    const T* getInternalPtr() const noexcept { return value.get(); }

    // Compare (assumes T::compare is const)
    signedInt compare(const BasisElement& other) const {
        return value->compare(*other.value);
    }
};