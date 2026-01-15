template<typename fieldType>
class compressed_sparse_matrix {
public:
    using indexType   = std::size_t;   // row indices, dimensions
    using offset_type = std::size_t;   // offsets into storage
    using k = fieldType;

    using Basis = BasisElement<indexType, k>;

    using DenseDomainVec = std::unique_ptr<k[]>;
    using DenseImageVec  = std::unique_ptr<k[]>;

    std::vector<Basis> rows_and_coeffs_;   // concatenated columns
    std::vector<offset_type> col_ptr_;     // size = ncols + 1
    indexType image_dim_{0};

public:
    // --- constructor ---
    explicit compressed_sparse_matrix(indexType image_dim)
        : image_dim_(image_dim)
    {
        col_ptr_.push_back(offset_type{0});  // sentinel for column 0
    }

    // --- dimensions ---
    indexType image_dim() const noexcept {
        return image_dim_;
    }

    indexType domain_dim() const noexcept {
        return static_cast<indexType>(col_ptr_.size() - 1);
    }

    // --- building ---
    void add_col(const std::vector<Basis>& col) {
        rows_and_coeffs_.insert(rows_and_coeffs_.end(), col.begin(), col.end());
        col_ptr_.push_back(static_cast<offset_type>(rows_and_coeffs_.size()));
    }

    // --- column view ---
    struct column_view {
        const Basis* b;
        const Basis* e;
        const Basis* begin() const { return b; }
        const Basis* end()   const { return e; }
    };

    column_view get_column(indexType i) const {
        const offset_type b = col_ptr_[i];
        const offset_type e = col_ptr_[i + 1];
        return { rows_and_coeffs_.data() + b,
                 rows_and_coeffs_.data() + e };
    }

    // --- dense helpers ---
    DenseDomainVec reserve_dense_domain_vec() const {
        return std::make_unique<k[]>(domain_dim());
    }

    DenseImageVec reserve_dense_image_vec() const {
        return std::make_unique<k[]>(image_dim_);
    }

    // --- x = M^T * y ---
    DenseDomainVec evaluate_transpose(const DenseImageVec& input) const {
        DenseDomainVec result = reserve_dense_domain_vec();
        const indexType ncols = domain_dim();

        #pragma omp parallel for schedule(static)
        for (std::int64_t ci = 0; ci < static_cast<std::int64_t>(ncols); ++ci) {
            const indexType col = static_cast<indexType>(ci);

            k acc = k{0};
            for (const auto& be : get_column(col)) {
                acc += input[be.getValue()] * be.getCoefficient();
            }
            result[col] = acc;
        }

        return result;
    }
};
