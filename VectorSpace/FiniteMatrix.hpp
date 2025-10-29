#pragma once
#include <memory>
#include <cstddef> // size_t
#include <span>
#include <iomanip>
#include <cassert>
#include <optional>
#include <vector>

namespace VectorSpace {

template<typename k>
class FiniteMatrix {
public:
    using Row    = std::unique_ptr<k[]>;   // length = n_cols
    using Col    = std::unique_ptr<k[]>;   // length = n_rows
    using Matrix = std::unique_ptr<k[]>;   // single contiguous block (column-major)

private:
    // ---- private data at the top ----
    size_t n_rows;
    size_t n_cols;
    Matrix M; // contiguous column-major storage

    static inline size_t idx(size_t r, size_t c, size_t n_rows_) { return c * n_rows_ + r; }

public:
    // Column-major contiguous storage: idx(r,c) = c*n_rows + r
    FiniteMatrix() noexcept : n_rows(0), n_cols(0), M(nullptr) {}
    FiniteMatrix(size_t rows, size_t cols)
        : n_rows(rows), n_cols(cols), M(std::make_unique<k[]>(rows * cols)) {}

    k&       operator()(size_t r, size_t c)       { return M[idx(r, c, n_rows)]; }
    const k& operator()(size_t r, size_t c) const { return M[idx(r, c, n_rows)]; }

    // pointer to column c (contiguous block of length n_rows)
    k*       col_ptr(size_t c)       { return M.get() + c * n_rows; }
    const k* col_ptr(size_t c) const { return M.get() + c * n_rows; }

    size_t rows() const noexcept { return n_rows; }
    size_t cols() const noexcept { return n_cols; }


    inline void set(size_t row, size_t col, k value) {
        if (row >= rows()) {
            cout << "invalid row";
        }

        if (col >= cols()) {
            cout << "invalid col";
        }
        M[row + n_rows * col] = value;
    }

    // ---------- y = M^T * x ----------
    // x has length n_rows, y has length n_cols.
    Row mul_T(const Col& x) const {
        Row y = std::make_unique<k[]>(n_cols);
        for (size_t j = 0; j < n_cols; ++j) {
            const k* cj = col_ptr(j); // contiguous column j
            k sum = k{};
            for (size_t r = 0; r < n_rows; ++r) {
                sum = sum + cj[r] * x[r];
            }
            y[j] = sum;
        }
        return y;
    }

    // Greedy solve for M^T * x ≈ y, LEFT→RIGHT with prefix preservation.
    // For each i, try to make y_hat[i] == y[i] without changing y_hat[0..i-1].
    // PRE: y points to n_cols values (not modified).
    Col solve(Row& y) const {
        Col x     = std::make_unique<k[]>(n_rows);
        Row ycur  = std::make_unique<k[]>(n_cols);
        auto used = std::make_unique<unsigned char[]>(n_rows);

        // init
        for (size_t r = 0; r < n_rows; ++r) { x[r] = k{}; used[r] = 0; }
        for (size_t j = 0; j < n_cols; ++j)  ycur[j] = k{};

        // columns: left → right
        for (size_t i = 0; i < n_cols; ++i) {
            const k residual = y[i] - ycur[i];
            if (residual == k{}) continue;                 // already matched

            const k* col_i = col_ptr(i);                   // column i (len n_rows)
            assert(col_i);

            // find an unused pivot row p with:
            //   M[p,i] != 0  AND  M[p,j] == 0 for all j < i   (won't disturb earlier matches)
            size_t p = n_rows;
            for (size_t r = 0; r < n_rows; ++r) {
                if (used[r] || col_i[r] == k{}) continue;

                bool ok = true;
                for (size_t j = 0; j < i; ++j) {
                    if (col_ptr(j)[r] != k{}) { ok = false; break; }
                }
                if (ok) { p = r; break; }                  // pick first/top-most; change order if desired
            }
            if (p == n_rows) continue;                     // can't fix this column without breaking prefix

            const k a     = col_i[p];
            const k delta = residual / a;
            if (delta == k{}) continue;                    // nothing to add

            x[p]    += delta;
            used[p]  = 1;

            // ycur += delta * (row p of M); j < i are guaranteed zero, so start at i
            for (size_t j = i; j < n_cols; ++j) {
                ycur[j] += col_ptr(j)[p] * delta;          // M[p,j] = col_ptr(j)[p]
            }
        }

        return x;  // move out
    }

    // ---------- Compose: y_hat = M^T * solve(y) ----------
    Row solve_then_mul_T(Row& y) const {
        Col x = solve(y);
        return mul_T(x);
    }

    // Solve M^T * x = y exactly via Gauss–Jordan elimination.
    // Returns std::nullopt if inconsistent; otherwise a Col (unique_ptr<k[]>)
    // containing one solution (free vars set to 0).
    std::optional<Col> solve_exact_MT(const Row& y) const {
        const size_t m = n_cols; // rows of A = M^T (number of equations)
        const size_t n = n_rows; // cols of A (number of unknowns)

        // Quick sanity: y must have length m (user responsibility in this design).
        // We proceed assuming that's true.

        // Build A = M^T (row-major here for simpler row ops)
        std::vector<k> A(m * n);
        for (size_t i = 0; i < m; ++i) {          // row i of A
            for (size_t j = 0; j < n; ++j) {      // col j of A
                // A[i,j] = (M^T)[i,j] = M[j,i]
                A[i * n + j] = M[idx(j, i, n_rows)];
            }
        }

        // Copy RHS b = y (length m)
        std::vector<k> b(m);
        for (size_t i = 0; i < m; ++i) b[i] = y[i];

        auto Aat = [&](size_t i, size_t j) -> k& { return A[i * n + j]; };

        // Gauss–Jordan elimination to RREF
        size_t row = 0;                   // current pivot row
        std::vector<size_t> pivot_col;    // pivot column per used row
        pivot_col.reserve(std::min(m, n));

        for (size_t col = 0; col < n && row < m; ++col) {
            // Find a pivot in column 'col' at or below 'row'
            size_t piv = row;
            while (piv < m && Aat(piv, col) == k{}) ++piv;
            if (piv == m) continue; // no pivot in this column → free variable

            // Swap the pivot row up if needed
            if (piv != row) {
                for (size_t j = 0; j < n; ++j) std::swap(Aat(piv, j), Aat(row, j));
                std::swap(b[piv], b[row]);
            }

            // Normalize pivot row so pivot becomes 1
            const k pivot = Aat(row, col);
            assert(pivot != k{});
            const k inv = k{1} / pivot;
            for (size_t j = 0; j < n; ++j) Aat(row, j) = Aat(row, j) * inv;
            b[row] = b[row] * inv;

            // Eliminate column 'col' in all other rows
            for (size_t i = 0; i < m; ++i) {
                if (i == row) continue;
                const k f = Aat(i, col);
                if (f == k{}) continue;
                for (size_t j = 0; j < n; ++j) Aat(i, j) = Aat(i, j) - f * Aat(row, j);
                b[i] = b[i] - f * b[row];
            }

            pivot_col.push_back(col);
            ++row;
        }

        // Check for inconsistency: 0 == row of A but b != 0
        for (size_t i = 0; i < m; ++i) {
            bool all_zero = true;
            for (size_t j = 0; j < n; ++j) {
                if (Aat(i, j) != k{}) { all_zero = false; break; }
            }
            if (all_zero && b[i] != k{}) {
                return std::nullopt; // no solution
            }
        }

        // Construct one solution: free variables = 0; basic = b at pivot rows.
        Col x = std::make_unique<k[]>(n);
        for (size_t j = 0; j < n; ++j) x[j] = k{};
        for (size_t r = 0; r < pivot_col.size(); ++r) {
            const size_t c = pivot_col[r];
            // After Gauss–Jordan, row r has 1 at column c and zeros elsewhere.
            x[c] = b[r];
        }

        return std::optional<Col>(std::move(x));
    }


    std::optional<Row> solve_exact_then_mul_T(Row& y) const {
            if (auto x = solve_exact_MT(y); x) {              // x is std::optional<Col>
                    return mul_T(*x);                    // mul_T takes const Col&; *x is a Col&
            }
            return std::nullopt;                     // no solution → empty
    }


    public:
        void print(std::ostream& os = std::cout,
                int width = 8,
                int precision = -1) const
        {
            if (precision >= 0) os << std::fixed << std::setprecision(precision);
            for (size_t r = 0; r < n_rows; ++r) {
                for (size_t c = 0; c < n_cols; ++c) {
                    os << std::setw(width) << M[idx(r, c, n_rows)];
                    if (c + 1 < n_cols) os << ' ';
                }
                os << '\n';
            }
        }



        void print_transpose(std::ostream& os = std::cout,
                int width = 8,
                int precision = -1) const
        {
            if (precision >= 0) os << std::fixed << std::setprecision(precision);
            
            for (size_t c = 0; c < n_cols; ++c) {
                for (size_t r = 0; r < n_rows; ++r) {
                    os << std::setw(width) << M[idx(r, c, n_cols)];
                    if (c + 1 < n_cols) os << ' ';
                }
                os << '\n';
            }
        }


};

} // namespace VectorSpace