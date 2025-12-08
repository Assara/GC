#pragma once

#include "HasCompare.hpp"
// #include <memory>   // no longer needed

template<typename T, typename k>
class BasisElement {
private:
    T value;         // stored by value, not via pointer
    k coefficient;

public:
    // --- ctors ---

    // Default-constructible if T and k are
    BasisElement() = default;

    // Construct from value (copy)
    explicit BasisElement(const T& val, k coeff = k{1})
        : value(val), coefficient(coeff) {}

    // Construct from value (move)
    explicit BasisElement(T&& val, k coeff = k{1})
        : value(std::move(val)), coefficient(coeff) {}

    // Copy / move are now trivial / defaulted
    BasisElement(const BasisElement&)            = default;
    BasisElement& operator=(const BasisElement&) = default;
    BasisElement(BasisElement&&) noexcept        = default;
    BasisElement& operator=(BasisElement&&) noexcept = default;


    BasisElement(std::unique_ptr<T>&&, k) = delete;

    // --- mutators ---

    void multiplyCoefficient(k factor) noexcept { coefficient *= factor; }

    // --- accessors ---

    const T* borrowValue() const noexcept { return &value; }
    T*       borrowValue()       noexcept { return &value; }

    const T& getValue() const noexcept { return value; }
    T&       getValue()       noexcept { return value; }

    k  getCoefficient() const noexcept { return coefficient; }
    k& getCoefficientRef()     noexcept { return coefficient; }

    // Optional "dereference" operators if you like pointer-like syntax:
    const T& operator*() const noexcept { return value; }
    T&       operator*()       noexcept { return value; }

    const T* operator->() const noexcept { return &value; }
    T*       operator->()       noexcept { return &value; }

    // --- comparisons (value-based) ---
	signedInt compare(const BasisElement& other) const noexcept {
		const T& a = value;
		const T& b = other.value;

		if constexpr (HasCompare<T>) {
			// Underlying value type has its own compare()
			return a.compare(b);
		} else {
			// Fallback: use relational operators (works for size_t, ints, etc.)
			// (a > b) - (a < b) → -1, 0, +1
			return static_cast<signedInt>((a > b) - (a < b));
		}
	}

    bool operator<(const BasisElement& rhs) const noexcept {
		return value < rhs.getValue();
    }

    bool operator==(const BasisElement& rhs) const noexcept {
        return value == rhs.getValue(); // equality by value (coeff ignored)
    }

    void set_coefficient(k c) noexcept { coefficient = c; }

    void add_if_same(BasisElement& other) {
        if (*this == other) {
            coefficient += other.getCoefficient();
            other.set_coefficient(k{0});
        }
    
    }
};
