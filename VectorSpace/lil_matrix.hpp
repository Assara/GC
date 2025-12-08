#pragma once

#include <vector>
#include <memory>
#include <cstddef>
#include <optional>

#include "VectorSpace/LinComb.hpp"  // adjust include path

template<typename k>
class lil_matrix {
public:
    using DomainVec = VectorSpace::LinComb<std::size_t, k>; // we are expecting indices to be in [0, cols.size())
    using ImageVec = VectorSpace::LinComb<std::size_t, k>;  // we are expecting indices to be in [0, image_dim)

private:
    std::vector<ImageVec> cols_;   // dynamic number of rows
    size_t image_dim_ = 0;

public:
    lil_matrix() = default;

    explicit lil_matrix(std::size_t reserve_rows) {
        cols_.reserve(reserve_rows);
    }

    // --- dimensions ---

    std::size_t image_dim() const noexcept {
		return image_dim_;
	}

    // Compute number of columns lazily: max column index + 1
    std::size_t domain_dim() const {
        return cols_.size();
    }

    // --- adding rows dynamically ---

    // Add an existing col (move)
	void add_col(ImageVec&& col) {
		if (col.size() > 0) {
			const auto& last = col.raw_elements().back();
			std::size_t row_index = last.getValue();  // image index
			image_dim_ = std::max(image_dim_, row_index + 1);
		}
		cols_.emplace_back(std::move(col));
	}

    inline std::size_t size() const {
        std::size_t result = 0;
        for (const ImageVec& col : cols_) {
            result += col.size();
        }
        return result;
    }

        
   // used for testing. could be smarter
    const ImageVec evaluate(DomainVec& input) const {
        ImageVec result;
		
		//could divide and conquer for better complexity
		for (const auto& be : input) {
			result += cols_[be.getValue()]*be.getCoefficient();
		
		}
        return result;
    }
    
    
    
    
    std::optional<DomainVec> solve_MX_equals_y(const ImageVec& y) const {
		const std::size_t m = image_dim_;      // number of rows (Image dimension)
		const std::size_t n = cols_.size();    // number of columns (Domain dimension)

		// Trivial cases
		if (m == 0 || n == 0) {
			// 0 * X = y  ⇒ solvable iff y == 0
			for (const auto& term : y) {
				if (term.getCoefficient() != k{}) {
					return std::nullopt;
				}
			}
			return DomainVec{}; // X = 0
		}

		// --- Build row-wise representation from columns ---

		// rows[i]: equation for row i, as a LinComb of domain indices (0..n-1)
		std::vector<DomainVec> rows(m);

		// For each column j (domain index), distribute its entries into rows[i]
		for (std::size_t j = 0; j < n; ++j) {
			const ImageVec& col = cols_[j];
			for (const auto& term : col) {
				std::size_t i = term.getValue();  // image index (row)
				if (i >= m) {
					// If this ever triggers, image_dim_ is inconsistent with col indices.
					continue;                     // or assert(false);
				}
				k coeff = term.getCoefficient();
				if (coeff != k{}) {
					// variable j appears in row i with coefficient 'coeff'
					rows[i].append_in_basis_order(j, coeff);
				}
			}
		}

		// --- Build RHS vector b from y (dense length m) ---

		std::vector<k> b(m, k{});
		for (const auto& term : y) {
			std::size_t i = term.getValue();      // row index
			if (i < m) {
				b[i] = term.getCoefficient();
			}
		}

		// Helper: get coefficient of variable 'col' in a row
		auto get_coeff = [](const DomainVec& row, std::size_t col) -> k {
			for (const auto& term : row) {
				std::size_t j = term.getValue();
				if (j == col) return term.getCoefficient();
				if (j > col) break;               // rows are kept sorted
			}
			return k{}; // zero
		};

		// pivot_row_for_col[j] = index of pivot row for variable j, or -1 if free
		std::vector<int> pivot_row_for_col(n, -1);

		std::size_t pivot_row = 0;

		// --- Sparse Gaussian elimination on rows ---

		for (std::size_t col = 0; col < n && pivot_row < m; ++col) {
			// find pivot row with nonzero coeff in column 'col'
			std::size_t pivot = pivot_row;
			while (pivot < m && get_coeff(rows[pivot], col) == k{}) {
				++pivot;
			}

			if (pivot == m) {
				// no pivot in this column ⇒ free variable
				continue;
			}

			// swap pivot row up
			if (pivot != pivot_row) {
				std::swap(rows[pivot], rows[pivot_row]);
				std::swap(b[pivot],    b[pivot_row]);
			}

			// normalize pivot row so that coeff in column 'col' is 1
			k pivot_val = get_coeff(rows[pivot_row], col);
			// Use whatever inverse you use for k:
			//   k inv_pivot = pivot_val.inv();
			// or, if appropriate:
			k inv_pivot = k{1} / pivot_val;

			rows[pivot_row].scalar_multiply(inv_pivot);
			b[pivot_row] *= inv_pivot;

			// eliminate this variable from all other rows
			for (std::size_t i = 0; i < m; ++i) {
				if (i == pivot_row) continue;

				k factor = get_coeff(rows[i], col);
				if (factor == k{}) continue;

				// row_i ← row_i - factor * row_pivot
				DomainVec tmp = rows[pivot_row];   // copy pivot row
				tmp.scalar_multiply(-factor);
				rows[i] += tmp;                    // uses LinComb::operator+=

				// *** THIS WAS MISSING BEFORE: update RHS as well ***
				b[i] -= factor * b[pivot_row];
			}

			pivot_row_for_col[col] = static_cast<int>(pivot_row);
			++pivot_row;
		}

		// --- Check for inconsistency: 0 ... 0 = nonzero ---

		auto is_zero_row = [](const DomainVec& row) {
			for (const auto& term : row) {
				if (term.getCoefficient() != k{}) return false;
			}
			return true;
		};

		for (std::size_t i = 0; i < m; ++i) {
			if (is_zero_row(rows[i]) && b[i] != k{}) {
				// 0 = nonzero ⇒ no solution
				return std::nullopt;
			}
		}

		// --- Read off one particular solution: free vars = 0 (no support minimization) ---

		DomainVec solution;
		solution.reserve(n);

		for (std::size_t col = 0; col < n; ++col) {
			int r = pivot_row_for_col[col];
			if (r < 0) {
				// free variable ⇒ choose 0
				continue;
			}	
			k x_col = b[static_cast<std::size_t>(r)];
			if (x_col != k{}) {
				solution.append_in_basis_order(col, x_col); // domain index
			}
		}

		solution.shrink_to_fit();
		return solution;
	}

	

    
   
};
