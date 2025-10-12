#pragma once
#include "LinComb.hpp"
#include "FiniteMatrix.hpp"   // VectorSpace::FiniteMatrix<k> (single-block, column-major)
#include <algorithm>
#include <vector>

namespace VectorSpace {

template<typename T, typename k>
requires ValidBasisElement<T, k>
class FiniteSubSpace {
public:
    using FM  = FiniteMatrix<k>;
    using Row = typename FM::Row;  // length = cols()
    using Col = typename FM::Col;  // length = rows()

    // OWNED canonical basis (copied), sorted & deduped by value
    std::vector<T> elements;

    // Dense matrix (value, column-major single block)
    FM M;

    // Build from a list of LinComb<T, k>; no sources kept afterward
    explicit FiniteSubSpace(const std::vector<LinComb<T, k>>& vecs)
        : elements(build_elements(vecs))
        , M(vecs.size(), elements.size())
    {
        // Fill each matrix row directly (no temporary Row)
        // Column-major: address of (row=j, col=0) is M.col_ptr(0) + j; row stride = M.rows()
        for (size_t j = 0; j < M.rows(); ++j) {
            k* row_base = M.col_ptr(0) + j;
            projectVec_into(vecs[j], row_base, M.rows());
        }
    }

    // sizes / access
    size_t rows() const noexcept { return M.rows(); }
    size_t cols() const noexcept { return M.cols(); }

    k&       at(size_t j, size_t i)       { return M(j, i); }
    const k& at(size_t j, size_t i) const { return M(j, i); }

    const std::vector<T>& basis() const noexcept { return elements; }
    
    LinComb<T, k> project_solve_then_map(const LinComb<T, k>& lc) const {
            // 1) Project lc → dense contiguous row y (length = M.cols())
            
            Row y = std::make_unique<k[]>(M.cols());
        
            projectVec_into(lc, y.get(), 1);


            // 2) Compute ŷ = Mᵀ * solve(y)
            Row y_hat =  M.solve_then_mul_T(y);

            LinComb<T, k> out{};         

            for (size_t i = 0; i < cols(); i++) {
                   // cout << "y[" << i << "] = " << y[i]; 
                    //cout << "                y_hat[" << i << "] = " << y_hat[i] << endl; 
                    // elements[i].print();

                    out.append_in_basis_order(elements[i], -y_hat[i]);
            }
            

            return out;
    }


private:
    // Build & dedupe basis elements (copied by value)
    static std::vector<T> build_elements(const std::vector<LinComb<T, k>>& vecs) {
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

   void projectVec_into(const LinComb<T,k>& lc, k* out, size_t stride) const {
        const size_t ncols = M.cols();
        // zero the row
        for (size_t j = 0; j < ncols; ++j)
            out[j * stride] = k{};

        const auto& ve = lc.raw_elements();
        if (ve.empty() || elements.empty() || ncols == 0) return;

        auto eq_val = [](const T& a, const T& b){ return a == b; };

        struct Frame { size_t v_lo, v_hi, e_lo, e_hi; };
        std::vector<Frame> st;
        st.push_back({0, ve.size(), 0, elements.size()});

        while (!st.empty()) {
            Frame f = st.back(); st.pop_back();
            if (f.v_lo >= f.v_hi || f.e_lo >= f.e_hi) continue;

            const size_t mid = (f.v_lo + f.v_hi) >> 1;
            const T& key = *ve[mid].borrowValue();

            auto it = std::lower_bound(elements.begin() + f.e_lo,
                                    elements.begin() + f.e_hi,
                                    key);
            const size_t pos = static_cast<size_t>(it - elements.begin());

            if (it != elements.begin() + f.e_hi && eq_val(*it, key)) {
                k sum = static_cast<k>(ve[mid].getCoefficient());

                size_t L = mid;
                while (L > f.v_lo && eq_val(*ve[L-1].borrowValue(), key)) {
                    --L; sum += static_cast<k>(ve[L].getCoefficient());
                }
                size_t R = mid + 1;
                while (R < f.v_hi && eq_val(*ve[R].borrowValue(), key)) {
                    sum += static_cast<k>(ve[R].getCoefficient()); ++R;
                }

                // stride-aware write
                assert(pos < ncols);
                out[pos * stride] += sum;

                if (R < f.v_hi && pos + 1 < f.e_hi)
                    st.push_back({R, f.v_hi, pos + 1, f.e_hi});
                if (f.v_lo < L && f.e_lo < pos)
                    st.push_back({f.v_lo, L, f.e_lo, pos});
            } else {
                if (f.v_lo < mid && f.e_lo < pos)
                    st.push_back({f.v_lo, mid, f.e_lo, pos});
                if (mid + 1 < f.v_hi && pos < f.e_hi)
                    st.push_back({mid + 1, f.v_hi, pos, f.e_hi});
            }
        }
    }


    public: //debug

        void printMatrix() {
            M.print();
        }


}; // namespace VectorSpace
}