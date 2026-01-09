#pragma once

#include "BasisElement.hpp"
#include "HasCanonized.hpp"
#include <algorithm>  
#include <vector>
#include <iostream>
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
		if constexpr (!HasCanonized<T,k>) {
			return;
		} else {	
			std::vector<Element> standardized;
			for (auto& elem : elements) {

				Element canon = elem.getValue().canonized(elem);
				if (canon.getCoefficient() != 0) {  
					standardized.emplace_back(std::move(canon));
				}
			}

			elements = std::move(standardized);
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
    
    void append_in_basis_order(const T& val, k coeff) {
        if (coeff == k{}) return;
        elements.emplace_back(val, coeff);
    }

    void append_in_basis_order(const Element& be) {
        elements.emplace_back(be);
    }
    
    
	 void append_scaled_without_sort(const LinComb& other, const k& scalar) {
		for (const auto& term : other) {
			k c = term.getCoefficient() * scalar;
			elements.emplace_back(term.getValue(), c);
		}
	}

    
    
	LinComb operator+(const LinComb& other) const {
		LinComb result;
		const auto& A = this->elements;
		const auto& B = other.elements;

		result.reserve(A.size() + B.size());

		std::size_t i = 0, j = 0;

		while (i < A.size() && j < B.size()) {
			const T& valA = A[i].getValue();
			const T& valB = B[j].getValue();

			signedInt cmp;
			if constexpr (HasCompare<T>) {
				cmp = valA.compare(valB);
			} else {
				cmp = static_cast<signedInt>((valA > valB) - (valA < valB));
			}

			if (cmp < 0) {
				result.append_in_basis_order(A[i]);
				++i;
			} else if (cmp > 0) {
				result.append_in_basis_order(B[j]);
				++j;
			} else {
				k sumCoeff = A[i].getCoefficient() + B[j].getCoefficient();
				if (sumCoeff != k{}) {
					result.append_in_basis_order(valA, sumCoeff);
				}
				++i;
				++j;
			}
		}

		while (i < A.size()) {
			result.append_in_basis_order(A[i++]);
		}
		while (j < B.size()) {
			result.append_in_basis_order(B[j++]);
		}

		return result;
	}

	LinComb& operator+=(const LinComb& other) {
		*this = *this + other;
		return *this;
	}
	
	LinComb add_scaled(const LinComb& other, k scalar) {
		LinComb result;
		if (scalar == k{} || other.size() == 0) {
			return *this; 
		}

		const auto& A = this->elements;
		const auto& B = other.elements;
		
		result.reserve(A.size() + B.size());

		std::size_t i = 0, j = 0;

		while (i < A.size() && j < B.size()) {
			const T& valA = A[i].getValue();
			const T& valB = B[j].getValue();

			signedInt cmp;
			if constexpr (HasCompare<T>) {
				cmp = valA.compare(valB);
			} else {
				cmp = static_cast<signedInt>((valA > valB) - (valA < valB));
			}

			if (cmp < 0) {
				result.append_in_basis_order(A[i]);
				++i;
			} else if (cmp > 0) {
				result.append_in_basis_order(B[j].getValue(), B[j].getCoefficient() * scalar);
				++j;
			} else {
				k sumCoeff = A[i].getCoefficient() + B[j].getCoefficient() * scalar;
				if (sumCoeff != k{}) {
					result.append_in_basis_order(valA, sumCoeff);
				}
				++i;
				++j;
			}
		}

		while (i < A.size()) {
			result.append_in_basis_order(A[i++]);
		}
		while (j < B.size()) {
			result.append_in_basis_order(B[j].getValue(), B[j].getCoefficient() * scalar);
			++j;
		}
		
		return result;
	}
	

    bigInt size() const {
        return elements.size();
    }
    
    bool empty() const {
        return elements.empty();
    }
    
    bool is_zero() const {
        return elements.empty();
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
        sort_elements();
    }
    
    void shrink_to_fit() {
		elements.shrink_to_fit();
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
	
    
    LinComb& operator*=(k scalar) {
			return scalar_multiply(scalar);
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


public:
    void print() const {
            std::cout << "LinCombPrint: "<< std::endl;

            for (const auto& elem : elements) {
                    elem.getValue().print();
                    std::cout << "Coefficient: " << elem.getCoefficient() << "\n\n";
            }
    }

};

} // namespace VectorSpace
