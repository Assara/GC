#pragma once

#include <vector>
#include <cstddef>
#include <optional>

#include "VectorSpace/LinComb.hpp"

template<typename k>
class row_based_lil_matrix {
public:
    using DomainVec = VectorSpace::LinComb<std::size_t, k>; // indices in [0, domain_dim_)
    using ImageVec  = VectorSpace::LinComb<std::size_t, k>; // indices in [0, rows_.size())

private:
    std::vector<DomainVec> rows_;   // row i is equation for image index i
    std::size_t domain_dim_ = 0;    // max column index + 1

public:
    row_based_lil_matrix() = default;

    explicit row_based_lil_matrix(std::size_t reserve_rows) {
        rows_.reserve(reserve_rows);
    }

    // --- dimensions ---

    // number of rows (image dimension)
    std::size_t image_dim() const noexcept {
        return rows_.size();
    }

    // number of columns (domain dimension)
    std::size_t domain_dim() const noexcept {
        return domain_dim_;
    }

    // --- adding entries dynamically ---

    // Add entry M[image_index, domain_index] += coefficient
    void add_element(std::size_t image_index,
                     std::size_t domain_index,
                     const k& coefficient)
    {
        if (coefficient == k{}) return;

        if (image_index >= rows_.size()) {
            rows_.resize(image_index + 1);  // extend with empty rows
        }

        rows_[image_index].append_in_basis_order(domain_index, coefficient);
        domain_dim_ = std::max(domain_dim_, domain_index + 1);
    }

    // Optionally standardize/sort each row
    void sort_rows() {
        for (auto& row : rows_) {
            row.standardize_and_sort();
        }
    }

    // total number of stored terms (nnz)
    std::size_t size() const {
        std::size_t result = 0;
        for (const DomainVec& row : rows_) {
            result += static_cast<std::size_t>(row.size());
        }
        return result;
    }

private:
    using DenseVector = std::vector<k>;

    static k get_coeff_in_col(const DomainVec& row, std::size_t col) {
        for (const auto& term : row) {
            std::size_t j = term.getValue();
            if (j == col) return term.getCoefficient();
            if (j > col) break; // basis-order invariant
        }
        return k{};
    }

    static bool is_zero_row(const DomainVec& row) {
        // Assuming LinComb invariant: no explicit zero coefficients
        return row.size() == 0;
    }


public:
    // Solve M * X = y, where:
    //  - M : Domain -> Image, stored row-wise in rows_
    //  - y : ImageVec, indices are row/image indices
    // Returns:
    //  - DomainVec X if solution exists
    //  - std::nullopt if inconsistent
    std::optional<DomainVec> solve_MX_equals_y(const ImageVec& y) {
        const std::size_t m = rows_.size();   // image dimension (#rows)
        const std::size_t n = domain_dim_;    // domain dimension (#cols)

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

        // --- Build dense RHS b from y (length m) ---

        DenseVector b(m, k{});
        for (const auto& term : y) {
            std::size_t i = term.getValue(); // row index
            if (i < m) {
                b[i] = term.getCoefficient();
            }
        }

        // --- Forward elimination on rows_ (destructive) ---

        // pivot_row_for_col[j] = row index of pivot in column j, or -1 if free
        std::vector<int> pivot_row_for_col(n, -1);
        // pivot_col_for_row[i] = pivot column in row i, or -1 if no pivot
        std::vector<int> pivot_col_for_row(m, -1);

        std::size_t pivot_row = 0;
        const std::size_t max_rank = std::min(m, n);

        for (std::size_t col = 0; col < n && pivot_row < max_rank; ++col) {
            // Find pivot row with nonzero coeff in column `col`
            std::size_t pivot = pivot_row;
            while (pivot < m && get_coeff_in_col(rows_[pivot], col) == k{}) {
                ++pivot;
            }
            if (pivot == m) {
                // No pivot in this column ⇒ free variable
                continue;
            }

            // Swap pivot row up
            if (pivot != pivot_row) {
                std::swap(rows_[pivot], rows_[pivot_row]);
                std::swap(b[pivot],    b[pivot_row]);
            }

            // Normalize pivot row so that coeff in column `col` is 1
            k pv = get_coeff_in_col(rows_[pivot_row], col);
            // Adjust depending on your field type:
            //   k inv_pivot = pv.inv();
            k inv_pivot = k{1} / pv;

            rows_[pivot_row].scalar_multiply(inv_pivot);
            b[pivot_row] *= inv_pivot;

            // Eliminate this variable from rows BELOW the pivot row
            for (std::size_t i = pivot_row + 1; i < m; ++i) {
                k factor = get_coeff_in_col(rows_[i], col);
                if (factor == k{}) continue;

                // row_i ← row_i - factor * row_pivot
                rows_[i] = rows_[i].add_scaled(rows_[pivot_row], -factor);
                b[i]    -= factor * b[pivot_row];
            }

            pivot_row_for_col[col]        = static_cast<int>(pivot_row);
            pivot_col_for_row[pivot_row]  = static_cast<int>(col);
            ++pivot_row;
        }

        const std::size_t rank = pivot_row;

        // --- Check for inconsistency: 0 ... 0 = nonzero ---

        for (std::size_t i = 0; i < m; ++i) {
            if (is_zero_row(rows_[i]) && b[i] != k{}) {
                // 0 = nonzero ⇒ no solution
                return std::nullopt;
            }
        }

        // --- Back substitution: solve for pivot variables, free vars = 0 ---

        std::vector<k> x(n, k{});  // dense domain vector (all zeros initially)

        if (rank > 0) {
            for (int i = static_cast<int>(rank) - 1; i >= 0; --i) {
                int pc = pivot_col_for_row[static_cast<std::size_t>(i)];
                if (pc < 0) continue; // shouldn't happen for i < rank

                std::size_t col = static_cast<std::size_t>(pc);

                // sum_{j != col} a_ij * x_j
                k sum = k{};
                const DomainVec& row = rows_[static_cast<std::size_t>(i)];
                for (const auto& term : row) {
                    std::size_t j = term.getValue();
                    if (j == col) continue;
                    if (x[j] != k{}) {
                        sum += term.getCoefficient() * x[j];
                    }
                }

                // pivot coeff is 1 after normalization
                x[col] = b[static_cast<std::size_t>(i)] - sum;
            }
        }

        // --- Convert dense x into sparse DomainVec ---

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

};
