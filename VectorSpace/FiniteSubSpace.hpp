#pragma once
#include "LinComb.hpp"
#include "FiniteMatrix.hpp"   // VectorSpace::FiniteMatrix<k> (single-block, column-major)
#include <algorithm>
#include <vector>

namespace VectorSpace {

template<ValidBasisElement T, typename k>
class FiniteSubSpace {
public:
    using FM  = FiniteMatrix<k>;
    using Row = typename FM::Row;  // length = cols()
    using Col = typename FM::Col;  // length = rows()

    // OWNED canonical basis (copied), sorted & deduped by value
    std::vector<T> elements;

    // Dense matrix (value, column-major single block)
    FM M;

    // Build from a list of LinComb<T>; no sources kept afterward
    explicit FiniteSubSpace(const std::vector<LinComb<T>>& vecs)
        : elements(build_elements(vecs))
        , M(vecs.size(), elements.size())
    {
        // Fill each matrix row directly (no temporary Row)
        // Column-major: address of (row=j, col=0) is M.col_ptr(0) + j; row stride = M.rows()
        for (size_t j = 0; j < M.rows(); ++j) {
            k* row_base = M.col_ptr(0) + j;
            projectVec_into(vecs[j], row_base);
        }
    }

    // sizes / access
    size_t rows() const noexcept { return M.rows(); }
    size_t cols() const noexcept { return M.cols(); }

    k&       at(size_t j, size_t i)       { return M(j, i); }
    const k& at(size_t j, size_t i) const { return M(j, i); }

    const std::vector<T>& basis() const noexcept { return elements; }

private:
    // Build & dedupe basis elements (copied by value)
    static std::vector<T> build_elements(const std::vector<LinComb<T>>& vecs) {
        std::vector<T> elems;
        size_t cap = 0;
        for (const auto& lc : vecs) cap += lc.raw_elements().size();
        elems.reserve(cap);

        for (const auto& lc : vecs) {
            for (const auto& be : lc.raw_elements()) {
                elems.push_back(*be.borrowValue()); // COPY value
            }
        }
        std::sort(elems.begin(), elems.end()); // needs T::operator<
        elems.erase(std::unique(elems.begin(), elems.end(),
                                [](const T& a, const T& b){ return a == b; }), // needs T::operator==
                    elems.end());
        return elems;
    }

    // Project a LinComb<T> onto the owned 'elements' basis and write directly into a strided matrix row.
    // row_out points to (row=j, col=0). Stride = M.rows(), ncols = M.cols().
    void projectVec_into(const LinComb<T>& lc, k* row_out) const {
    // Column-major FiniteMatrix: consecutive columns in a fixed row are spaced by rows()
        const size_t stride = M.rows();
        const size_t ncols  = M.cols();

        // Zero the destination row (length = ncols) with the correct stride
        for (size_t i = 0; i < ncols; ++i) row_out[i * stride] = k{};

        const auto& ve = lc.raw_elements(); // assumed sorted by value; duplicates allowed
        if (ve.empty() || elements.empty() || ncols == 0) return;

        auto eq_val = [](const T& a, const T& b){ return a == b; };

        struct Frame { size_t v_lo, v_hi, e_lo, e_hi; };
        std::vector<Frame> st;
        st.push_back({0, ve.size(), 0, ncols});

        while (!st.empty()) {
            Frame f = st.back(); st.pop_back();
            if (f.v_lo >= f.v_hi || f.e_lo >= f.e_hi) continue;

            // Pick a pivot term from lc’s slice and locate it in elements[f.e_lo, f.e_hi)
            const size_t mid = (f.v_lo + f.v_hi) >> 1;
            const T& key = *ve[mid].borrowValue();

            auto it  = std::lower_bound(elements.begin() + f.e_lo,
                                        elements.begin() + f.e_hi,
                                        key);                  // uses T::operator<
            const size_t pos = static_cast<size_t>(it - elements.begin());

            if (it != elements.begin() + f.e_hi && eq_val(*it, key)) {
                // Key exists in basis -> fold duplicates around mid and write sum
                k sum = static_cast<k>(ve[mid].getCoefficient());

                size_t L = mid;
                while (L > f.v_lo && eq_val(*ve[L - 1].borrowValue(), key)) {
                    --L; sum = sum + static_cast<k>(ve[L].getCoefficient());
                }
                size_t R = mid + 1;
                while (R < f.v_hi && eq_val(*ve[R].borrowValue(), key)) {
                    sum = sum + static_cast<k>(ve[R].getCoefficient()); ++R;
                }

                // Write directly into row_out at column `pos`
                row_out[pos * stride] = row_out[pos * stride] + sum;

                // Recurse on the two remaining rectangles that could still match
                if (R < f.v_hi && pos + 1 < f.e_hi) st.push_back({R, f.v_hi, pos + 1, f.e_hi});
                if (f.v_lo < L && f.e_lo < pos)     st.push_back({f.v_lo, L, f.e_lo, pos});
            } else {
                // Key not in basis slice -> split search rectangles and continue (effectively ignore this key)
                if (f.v_lo < mid && f.e_lo < pos)           st.push_back({f.v_lo, mid, f.e_lo, pos});
                if (mid + 1 < f.v_hi && pos < f.e_hi)       st.push_back({mid + 1, f.v_hi, pos, f.e_hi});
            }
        }
    }


    LinComb<T> project_solve_then_map(const LinComb<T>& lc) const {
        // 1) Project lc → dense contiguous row y (length = M.cols())
        Row y = std::make_unique<k[]>(M.cols());
        projectVec_into(lc, y.get());

        // 2) Compute ŷ = Mᵀ * solve(y)
        Row y_hat = M.solve_then_mul_T(std::move(y));

        // 3) Build LinComb directly (no external vector), assuming basis order
        LinComb<T> out(LinComb<T>::AssumeBasisOrder);          // tag-only ctor: no sort/standardize
        // If available, you can reserve to avoid reallocations:
        // out.reserve(M.cols());

        for (size_t i = 0; i < M.cols(); ++i) {
            const k coeff = y_hat[i];
            if (coeff != k{}) {
                // Append in basis order; no canonicalization

                //fix the casting later
                out.append_in_basis_order(elements[i], static_cast<fieldType>(coeff));
            
        }
        return out;
    }

};

} // namespace VectorSpace
