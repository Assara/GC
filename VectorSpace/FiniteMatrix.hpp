#pragma once
#include <memory>
#include <cstddef> // size_t

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

    // ---------- Solve M^T * x ≈ y (greedy, right-to-left) ----------
    // y has length n_cols; returns x of length n_rows.
    Col solve(Row y) const {
        Col x    = std::make_unique<k[]>(n_rows);
        Row ycur = std::make_unique<k[]>(n_cols);
        auto used = std::make_unique<bool[]>(n_rows);

        for (size_t r = 0; r < n_rows; ++r) { x[r] = k{}; used[r] = false; }
        for (size_t j = 0; j < n_cols; ++j)  ycur[j] = k{};

        for (size_t i = n_cols; i-- > 0; ) {
            const k* col_i = col_ptr(i);     // column i contiguous
            k residual = y[i] - ycur[i];

            // choose bottom-most unused nonzero pivot in column i
            size_t p = n_rows;
            for (size_t r = n_rows; r-- > 0; ) {
                if (!used[r] && col_i[r] != k{}) { p = r; break; }
            }
            if (p == n_rows) continue;       // can't satisfy y[i]

            const k a = col_i[p];
            const k delta = residual / a;
            if (delta != k{}) {
                x[p]   = x[p] + delta;
                used[p] = true;

                // ycur += delta * (row p of M, transposed): ycur[j] += M[p][j] * delta
                for (size_t j = 0; j < n_cols; ++j) {
                    ycur[j] = ycur[j] + M[idx(p, j, n_rows)] * delta;
                }
            } else {
                used[p] = true; // keep pivot uniqueness
            }
        }
        return x;
    }

    // ---------- Compose: y_hat = M^T * solve(y) ----------
    Row solve_then_mul_T(Row y) const {
        Col x = solve(std::move(y)); // x ∈ k^{n_rows}
        return mul_T(x);             // y_hat ∈ k^{n_cols}
    }
};

} // namespace VectorSpace