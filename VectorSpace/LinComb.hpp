#pragma once

#include "BasisElement.hpp"
#include "ValidForDifferential.hpp"
#include <vector>
#include <memory>

namespace VectorSpace {

template<ValidBasisElement T>
class LinComb {
public:
    using Element = BasisElement<T>;
    using ElementPtr = std::unique_ptr<Element>;
    
    void standardize_all() {
        for (auto& elem : elements) {
            elem-> getValue().std(*elem);
        }
    }

    explicit LinComb(const BasisElement<T>& elem) {
        elements.push_back(std::make_unique<Element>(elem));
    }

    LinComb& operator+=(LinComb const& other) {
        std::vector<ElementPtr> result;
        result.reserve(this->elements.size() + other.elements.size());

        std::size_t i = 0, j = 0;
        auto& A = this->elements;
        auto& B = other.elements;

        // Merge‐scan loop:
        while (i < A.size() && j < B.size()) {
            T const& valA = A[i]->getValue();
            T const& valB = B[j]->getValue();
            signedInt cmp = valA.compare(valB);

            if (cmp < 0) {
                // A[i] < B[j], so copy A[i] into result
                result.push_back(
                    std::make_unique<Element>(valA, A[i]->getCoefficient())
                );
                ++i;
            }
            else if (cmp > 0) {
                // A[i] > B[j], so copy B[j] into result
                result.push_back(
                    std::make_unique<Element>(valB, B[j]->getCoefficient())
                );
                ++j;
            }
            else {
                // valA == valB: combine coefficients
                fieldType sumCoeff = A[i]->getCoefficient() + B[j]->getCoefficient();
                if (sumCoeff != static_cast<fieldType>(0)) {
                    // Only insert if the new coefficient is nonzero
                    result.push_back(
                        std::make_unique<Element>(valA, sumCoeff)
                    );
                }
                ++i;
                ++j;
            }
        }

        // Copy any remaining from A
        while (i < A.size()) {
            T const& v = A[i]->getValue();
            result.push_back(
                std::make_unique<Element>(v, A[i]->getCoefficient())
            );
            ++i;
        }

        // Copy any remaining from B
        while (j < B.size()) {
            T const& v = B[j]->getValue();
            result.push_back(
                std::make_unique<Element>(v, B[j]->getCoefficient())
            );
            ++j;
        }

        // Swap in the merged vector as our new elements
        this->elements.swap(result);
        return *this;
    }

    //--------------------------------------------------------------------------
    // Out‐of‐place addition: returns A + B as a new LinComb
    LinComb operator+(LinComb const& other) const {
        LinComb copy = *this;   // make a shallow copy of our own elements
        copy += other;          // merge‐add in 'other'
        return copy;
    }

    //--------------------------------------------------------------------------
    // Optional: for convenience, you can expose a method 'add' that is exactly
    // the same as operator+=.
    void add(LinComb const& other) {
        *this += other;
    }

    auto begin() const {
        return elements.begin();
    }

    auto end() const {
        return elements.end();
    }

private:
    std::vector<ElementPtr> elements;
};

} // namespace VectorSpace
