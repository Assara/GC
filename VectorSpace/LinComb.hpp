#pragma once

#include "BasisElement.hpp"
#include "ValidBasisElement.hpp"
#include <vector>
#include "tags.hpp"

namespace VectorSpace {

template<typename T, typename k>
class LinComb {
public:
    using Element = BasisElement<T, k>;

private:
    std::vector<Element> elements;

public:
    LinComb() = default;

    void standardize_all() {
        for (auto& elem : elements) {
            elem.getValue().std(elem);
        }

    }

    explicit LinComb(const BasisElement<T, k>& elem) {
        elements.push_back(elem);
        standardize_all();
    }

    
    explicit LinComb(const BasisElement<T, k>& elem, AssumeBasisOrderTag) {
        elements.push_back(elem);
    }

    explicit LinComb(std::vector<Element>&& elems)
            : elements(std::move(elems)) {
            standardize_and_sort();
    }

    explicit LinComb(std::vector<std::unique_ptr<T>>&& graphs, const k& coeff) {
        elements.reserve(graphs.size());
        for (auto& p : graphs) {
            if (p && coeff != k{}) elements.emplace_back(std::move(p), coeff);
        }
        standardize_and_sort();
    }

    explicit LinComb(const T& val, k coeff = k{1}) {
            elements.emplace_back(val, coeff);
            standardize_and_sort();
    }

     explicit LinComb(std::vector<Element>&& elems, AssumeBasisOrderTag) noexcept
        : elements(std::move(elems)) {}


    explicit LinComb(typename std::vector<Element>::const_iterator first,
            typename std::vector<Element>::const_iterator last)
        : elements(first, last) {}

    
    explicit LinComb(const T& val, AssumeBasisOrderTag, const k& coeff = k{1}) {
            elements.emplace_back(val, coeff);
    }

    LinComb& operator+=(const LinComb& other) {
        std::vector<Element> result;
        result.reserve(this->elements.size() + other.elements.size());

        std::size_t i = 0, j = 0;
        const auto& A = this->elements;
        const auto& B = other.elements;

        // Merge-scan loop (reversed priority)
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
                k sumCoeff = A[i].getCoefficient() + B[j].getCoefficient();
                if (sumCoeff != static_cast<k>(0)) {
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

    bigInt size() const {
        return elements.size();
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


    BasisElement<T, k>& back() {
            return elements.back();
    }

    BasisElement<T, k>& front() {
            return elements.front();
    }

    inline void reserve(size_t size) {
            elements.reserve(size);
    }

    inline void sort_without_deduplicate() {
        std::sort(elements.begin(), elements.end());
    }

    void remove_zeros() {
        elements.erase(
            std::remove_if(
                elements.begin(),
                elements.end(),
                [](const Element& e) {
                    return e.getCoefficient() == static_cast<k>(0);
                }
            ),
            elements.end()
        );
    }

    void sort_elements() {
        sort_without_deduplicate();
        deduplicate_sorted();
    }
        
    void deduplicate_sorted() {
        if (elements.size() == 0) return;
        for (std::size_t i = 0; i < elements.size()-1; ++i) {
           elements[i+1].add_if_same(elements[i]);
        }
        remove_zeros();
    }

    void standardize_and_sort() {
        standardize_all();
        cout << "after standardizing:" << endl;
        print();
        sort_elements();
    }

    LinComb& scalar_multiply(k scalar) {
            if (scalar == 0) {
                    elements.clear();
                    return *this;
            }

            for (auto& be : elements) {
                    be.getCoefficientRef() *= scalar;
            }
            return *this;
    }

    LinComb operator*(k scalar) const {
        LinComb result = *this;
        result.scalar_multiply(scalar);
        return result;
    }

    const std::vector<Element>& raw_elements() const {
        return elements;
    }

    std::vector<Element>& raw_elements_nonconst() {
        return elements;
    }


    void append_in_basis_order(const T& val, k coeff) {
        if (coeff == k{}) return;
        elements.emplace_back(val, coeff);
    }


    void append_in_basis_order(Element& be) {
        elements.emplace_back(be);
    }


public:
    void print() const {
            cout << "LinCombPrint: "<< endl;

            for (const auto& elem : elements) {
                    elem.getValue().print();
                    std::cout << "Coefficient: " << elem.getCoefficient() << "\n\n";
            }
    }

};

} // namespace VectorSpace
