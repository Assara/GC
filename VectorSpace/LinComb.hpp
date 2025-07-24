#pragma once

#include "BasisElement.hpp"
#include "ValidForDifferential.hpp"
#include <vector>

namespace VectorSpace {

template<ValidBasisElement T>
class LinComb {
public:
    using Element = BasisElement<T>;

    LinComb() = default;

    void standardize_all() {
        for (auto& elem : elements) {
            elem.getValue().std(elem);
        }
    }

    explicit LinComb(const BasisElement<T>& elem) {
        elements.push_back(elem);
    }

    LinComb& operator+=(const LinComb& other) {
        std::vector<Element> result;
        result.reserve(this->elements.size() + other.elements.size());

        std::size_t i = 0, j = 0;
        const auto& A = this->elements;
        const auto& B = other.elements;

        // Merge-scan loop:
        while (i < A.size() && j < B.size()) {
            const T& valA = A[i].getValue();
            const T& valB = B[j].getValue();
            signedInt cmp = valA.compare(valB);

            if (cmp < 0) {
                result.push_back(A[i]);
                ++i;
            } else if (cmp > 0) {
                result.push_back(B[j]);
                ++j;
            } else {
                fieldType sumCoeff = A[i].getCoefficient() + B[j].getCoefficient();
                if (sumCoeff != static_cast<fieldType>(0)) {
                    result.emplace_back(valA, sumCoeff);
                }
                ++i;
                ++j;
            }
        }

        while (i < A.size()) {
            result.push_back(A[i]);
            ++i;
        }
        while (j < B.size()) {
            result.push_back(B[j]);
            ++j;
        }

        this->elements.swap(result);
        return *this;
    }

    LinComb operator+(const LinComb& other) const {
        LinComb copy = *this;
        copy += other;
        return copy;
    }

    void add(const LinComb& other) {
        *this += other;
    }

    auto begin() const { return elements.begin(); }
    auto end() const { return elements.end(); }

    void sort_elements() {
        if (elements.size() <= 1) return;

        std::vector<LinComb> wrapped;
        wrapped.reserve(elements.size());
        for (auto& e : elements) {
            if (e.getCoefficient() != static_cast<fieldType>(0)) {
                LinComb single;
                single.elements.push_back(std::move(e));
                wrapped.push_back(std::move(single));
            }
        }
        elements.clear();

        while (wrapped.size() > 1) {
            std::vector<LinComb> merged;
            for (size_t i = 0; i + 1 < wrapped.size(); i += 2) {
                wrapped[i] += wrapped[i + 1];
                merged.push_back(std::move(wrapped[i]));
            }
            if (wrapped.size() % 2 != 0) {
                merged.push_back(std::move(wrapped.back()));
            }
            wrapped = std::move(merged);
        }

        if (!wrapped.empty()) {
            *this = std::move(wrapped.front());
        }
    }

    void standardize_and_sort() {
        standardize_all();
        sort_elements();
    }

    const std::vector<Element>& raw_elements() const {
        return elements;
    }

private:
    std::vector<Element> elements;
};

} // namespace VectorSpace
