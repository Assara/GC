#pragma once

#include <vector>
#include <memory>
#include <cstddef>
#include <optional>


#include "LinComb.hpp"
#include "compressed_sparse_matrix.hpp"



template<typename k>
class lil_matrix {
public:
    using DomainVec = VectorSpace::LinComb<std::size_t, k>; // we are expecting indices to be in [0, cols.size())
    using ImageVec = VectorSpace::LinComb<std::size_t, k>;  // we are expecting indices to be in [0, image_dim)
    
    using DenseDomainVec = std::unique_ptr<k[]>; //assumed to hold cols.size() elements
    using DenseImageVec = std::unique_ptr<k[]>; //assumed to hold  image_dim elements

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
	
	void add_element(std::size_t image_index,
                     std::size_t domain_index,
                     const k& coefficient) {
        if (coefficient == k{}) return;

        if (domain_index >= cols_.size()) {
            cols_.resize(domain_index + 1);  // extend with empty rows
        }

        cols_[domain_index].append_in_basis_order(image_index, coefficient);
        image_dim_ = std::max(image_dim_, image_index + 1);
    }
	

    inline std::size_t size() const {
        std::size_t result = 0;
        for (const ImageVec& col : cols_) {
            result += col.size();
        }
        return result;
    }
    
    DenseImageVec reserve_dense_image_vec() const {
			return std::make_unique<k[]>(image_dim());
	}
	
	DenseDomainVec reserve_dense_domain_vec() const {
			return std::make_unique<k[]>(cols_.size());
	}
	
	DenseImageVec make_dense_image_vec_zero() const {
		DenseImageVec v = reserve_dense_image_vec();
		for (std::size_t i = 0; i < image_dim(); ++i) {
			v[i] = k{0};
		}
		return v;
	}
	
	DenseDomainVec make_dense_domain_vec_zero() const {
		DenseDomainVec v = reserve_dense_domain_vec();
		for (std::size_t i = 0; i < domain_dim(); ++i) {
			v[i] = k{0};
		}
		return v;
	}
	
	DenseImageVec make_dense_image_vec_random(std::uint64_t seed = 0xC0FFEEULL) const {		     
        DenseImageVec rv = reserve_dense_image_vec();

        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<std::uint64_t> dist(0, k::modulus() - 1);

        for (std::size_t i = 0; i < image_dim(); ++i) {
            rv[i] = k(dist(rng));
        }
        return rv;

	}
	
	void sort_cols() {
        for (auto& col : cols_) {
            col.standardize_and_sort();
        }
    }
	
	
	DenseDomainVec make_dense_domain_vec_random(std::uint64_t seed = 0xC0FFEEULL) const {
        DenseDomainVec rv = reserve_dense_domain_vec();

        std::mt19937_64 rng(seed);
        std::uniform_int_distribution<std::uint64_t> dist(0, k::modulus() - 1);

        for (std::size_t i = 0; i < domain_dim(); ++i) {
            rv[i] = k(dist(rng));
        }
        return rv;

	}



	
    
	// y = M * x
	DenseImageVec evaluate_from_dense(const DenseDomainVec& input) const {
			DenseImageVec result = make_dense_image_vec_zero();  // size = image_dim(), zero-initialized

			for (std::size_t i = 0; i < domain_dim(); ++i) {     // loop columns
					const auto& col = cols_[i];
					for (const auto& be : col) {
							std::size_t row = be.getValue();             // row index
							result[row] += input[i] * be.getCoefficient();
					}
			}

			return result;   // moved
	}

	// x = M^T * y
	DenseDomainVec evaluate_transpose_dense(const DenseImageVec& input) const {
		DenseDomainVec result = reserve_dense_domain_vec(); // sized to domain_dim(); no need to be zeroed

		const std::size_t n = domain_dim();

		#pragma omp parallel for schedule(static)
		for (std::int64_t ii = 0; ii < static_cast<std::int64_t>(n); ++ii) {
			const std::size_t i = static_cast<std::size_t>(ii);
			const auto& col = cols_[i];

			k acc = k{0};
			for (const auto& be : col) {
				const std::size_t row = be.getValue(); // row index
				acc += input[row] * be.getCoefficient();
			}
			result[i] = acc;
		}

		return result; // moved
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
    
    // Return the transposed matrix Mt = M^T.
	// M : image_dim_ x domain_dim()
	// Mt: domain_dim() x image_dim_
	lil_matrix transpose() const {
		const std::size_t m = image_dim_;      // rows in M
		const std::size_t n = domain_dim();    // cols in M

		lil_matrix Mt(m); // reserve rows (= domain indices) up to m

		for (std::size_t j = 0; j < n; ++j) {
			const auto& col = cols_[j];
			for (const auto& term : col) {
				const std::size_t i = term.getValue(); // row index in M
				Mt.add_element(
					/*image_index=*/ j,                 // row in Mt
					/*domain_index=*/ i,                // col in Mt
					/*coefficient=*/ term.getCoefficient()
				);
			}
		}

		return Mt;
	}
    
    
    compressed_sparse_matrix<k> to_compressed_sparse_matrix() const {
		using CSM = compressed_sparse_matrix<k>;

		CSM M(static_cast<typename CSM::indexType>(image_dim_));
		M.rows_and_coeffs_.reserve(size());          // total nnz
		M.col_ptr_.reserve(domain_dim() + 1);         // columns + sentinel

		for (const auto& col : cols_) {
			M.add_col_size_t(col.raw_elements());
		}

		return M;
	}


private: //solver helpers
    using RowSystem   = std::vector<DomainVec>;
    using DenseVector = std::vector<k>;

    // Coefficient of variable `col` in a row
    static k get_coeff_in_col(const DomainVec& row, std::size_t col) {
        for (const auto& term : row) {
            std::size_t j = term.getValue();
            if (j == col) return term.getCoefficient();
            if (j > col) break; // basis-order invariant
        }
        return k{};
    }

    // Build row-wise representation rows[i] from cols_
    void build_rows_from_columns(RowSystem& rows) const {
        const std::size_t m = image_dim_;
        const std::size_t n = cols_.size();

        rows.assign(m, DomainVec{});

        for (std::size_t j = 0; j < n; ++j) {
            const ImageVec& col = cols_[j];
            for (const auto& term : col) {
                std::size_t i = term.getValue();  // row index in image
                if (i >= m) continue;             // or assert
                k coeff = term.getCoefficient();
                if (coeff != k{}) {
                    rows[i].append_in_basis_order(j, coeff); // var j
                }
            }
        }
    }

    // Build dense RHS b from sparse y
    void build_rhs_from_y(const ImageVec& y, DenseVector& b) const {
        const std::size_t m = image_dim_;
        b.assign(m, k{});
        for (const auto& term : y) {
            std::size_t i = term.getValue(); // row index
            if (i < m) {
                b[i] = term.getCoefficient();
            }
        }
    }

    // Forward elimination: put (rows, b) into row-echelon form
    // Returns the rank and fills pivot_row_for_col, pivot_col_for_row.
    std::size_t forward_elimination(
        RowSystem& rows,
        DenseVector& b,
        std::vector<int>& pivot_row_for_col,
        std::vector<int>& pivot_col_for_row
    ) const {
        const std::size_t m = rows.size();
        const std::size_t n = cols_.size();

        pivot_row_for_col.assign(n, -1);
        pivot_col_for_row.assign(m, -1);

        std::size_t pivot_row = 0;
        const std::size_t max_rank = std::min(m, n);

        for (std::size_t col = 0; col < n && pivot_row < max_rank; ++col) {
            // Find pivot row with nonzero coeff in column `col`
            std::size_t pivot = pivot_row;
            while (pivot < m && get_coeff_in_col(rows[pivot], col) == k{}) {
                ++pivot;
            }
            if (pivot == m) {
                // No pivot in this column ⇒ free variable
                continue;
            }

            // Swap pivot row up
            if (pivot != pivot_row) {
                std::swap(rows[pivot], rows[pivot_row]);
                std::swap(b[pivot],    b[pivot_row]);
            }

            // Normalize pivot row so that coeff in column `col` is 1
            k pivot_val = get_coeff_in_col(rows[pivot_row], col);
            // Adjust depending on your field type:
            //   k inv_pivot = pivot_val.inv();
            k inv_pivot = k{1} / pivot_val;

            rows[pivot_row].scalar_multiply(inv_pivot);
            b[pivot_row] *= inv_pivot;

            // Eliminate this variable from rows BELOW the pivot row
            for (std::size_t i = pivot_row + 1; i < m; ++i) {
                k factor = get_coeff_in_col(rows[i], col);
                if (factor == k{}) continue;

                // row_i ← row_i - factor * row_pivot
                rows[i] = rows[i].add_scaled(rows[pivot_row], -factor);
                b[i]   -= factor * b[pivot_row];
            }

            pivot_row_for_col[col]       = static_cast<int>(pivot_row);
            pivot_col_for_row[pivot_row] = static_cast<int>(col);
            ++pivot_row;
        }

        return pivot_row; // rank
    }

    // Check for inconsistency: row = 0 but RHS != 0
    static bool has_inconsistency(const RowSystem& rows, const DenseVector& b) {
        const std::size_t m = rows.size();

        auto is_zero_row = [](const DomainVec& row) {
            for (const auto& term : row) {
                if (term.getCoefficient() != k{}) return false;
            }
            return true;
        };

        for (std::size_t i = 0; i < m; ++i) {
            if (is_zero_row(rows[i]) && b[i] != k{}) {
                return true; // 0 = nonzero ⇒ no solution
            }
        }
        return false;
    }

    // Back-substitution given row-echelon form
    DomainVec back_substitute(
        const RowSystem& rows,
        const DenseVector& b,
        const std::vector<int>& pivot_col_for_row,
        std::size_t rank
    ) const {
        const std::size_t n = cols_.size();
        std::vector<k> x(n, k{}); // dense domain vector

        if (rank > 0) {
            // Work from last pivot row up to first
            for (int i = static_cast<int>(rank) - 1; i >= 0; --i) {
                int pc = pivot_col_for_row[static_cast<std::size_t>(i)];
                if (pc < 0) continue; // shouldn't happen for i < rank

                std::size_t col = static_cast<std::size_t>(pc);

                // sum_{j != col} a_ij * x_j
                k sum = k{};
                const DomainVec& row = rows[static_cast<std::size_t>(i)];
                for (const auto& term : row) {
                    std::size_t j = term.getValue();
                    if (j == col) continue;
                    if (x[j] != k{}) {
                        sum += term.getCoefficient() * x[j];
                    }
                }

                // pivot coefficient is 1 after normalization
                x[col] = b[static_cast<std::size_t>(i)] - sum;
            }
        }

        // Convert dense x into sparse DomainVec
        DomainVec solution;
        solution.reserve(n);
        for (std::size_t col = 0; col < n; ++col) {
            if (x[col] != k{}) {
                solution.append_in_basis_order(col, x[col]);
            }
        }
        solution.shrink_to_fit();
        return solution;
    }
    
public:
    std::optional<DomainVec> solve_MX_equals_y(const ImageVec& y) const {
        const std::size_t m = image_dim_;
        const std::size_t n = cols_.size();

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

        // 1. Build rows and RHS
        RowSystem rows;
        DenseVector b;
        build_rows_from_columns(rows);
        build_rhs_from_y(y, b);

        // 2. Forward elimination
        std::vector<int> pivot_row_for_col;
        std::vector<int> pivot_col_for_row;
        std::size_t rank = forward_elimination(
            rows, b, pivot_row_for_col, pivot_col_for_row
        );

        // 3. Inconsistency check
        if (has_inconsistency(rows, b)) {
            return std::nullopt;
        }

        // 4. Back-substitution (free vars = 0)
        DomainVec solution = back_substitute(rows, b, pivot_col_for_row, rank);
        return solution;
    }
};
