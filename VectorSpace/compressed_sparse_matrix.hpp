#pragma once

#include <cstdint>
#include <cstddef>
#include <memory>
#include <vector>

#include "BasisElement.hpp"
#include "Field.hpp"

template<VectorSpace::Field fieldType>
class compressed_sparse_matrix {
	public:
		using indexType   = std::uint32_t; // 32-bit row indices, dimensions
		using offset_type = std::size_t;   // offsets into storage (nnz addressing)
		using k = fieldType;

		using Basis = BasisElement<indexType, k>;

		using DenseDomainVec = std::unique_ptr<k[]>;
		using DenseImageVec  = std::unique_ptr<k[]>;

		std::vector<Basis> rows_and_coeffs_;   // concatenated columns (CSC-like)
		std::vector<offset_type> col_ptr_;     // size = ncols + 1
		indexType image_dim_{0};

	public:
		explicit compressed_sparse_matrix(indexType image_dim)
			: image_dim_(image_dim)
		{
			col_ptr_.push_back(offset_type{0});  // sentinel for column 0
		}

		indexType image_dim() const noexcept {
			return image_dim_;
		}

		indexType domain_dim() const noexcept {
			return static_cast<indexType>(col_ptr_.size() - 1);
		}


		Basis indexType_conversion(const BasisElement<std::size_t,fieldType>& be) const {
			return Basis(static_cast<indexType>(be.getValue()), be.getCoefficient());
		}

		void add_col_size_t(const std::vector<BasisElement<std::size_t,fieldType>>& col) {
			for (const auto &be : col) {
				rows_and_coeffs_.emplace_back(indexType_conversion(be));
			}

			col_ptr_.push_back(static_cast<offset_type>(rows_and_coeffs_.size()));
		}


		void add_col(const std::vector<Basis>& col) {
			rows_and_coeffs_.insert(rows_and_coeffs_.end(), col.begin(), col.end());
			col_ptr_.push_back(static_cast<offset_type>(rows_and_coeffs_.size()));
		}

		struct column_view {
			const Basis* b;
			const Basis* e;
			const Basis* begin() const { return b; }
			const Basis* end()   const { return e; }
		};

		column_view get_column(indexType i) const {
			const std::size_t is = static_cast<std::size_t>(i);
			const offset_type b = col_ptr_[is];
			const offset_type e = col_ptr_[is + 1];
			if (rows_and_coeffs_.empty()) {
				return {nullptr, nullptr};
			}

			return { rows_and_coeffs_.data() + b,
				rows_and_coeffs_.data() + e };
		}

        compressed_sparse_matrix transpose() const {
            compressed_sparse_matrix result(domain_dim());
            result.col_ptr_.assign(std::size_t(image_dim())+1,0);
            for(const auto& term:rows_and_coeffs_)++result.col_ptr_[std::size_t(term.getValue())+1];
            for(std::size_t i=1;i<result.col_ptr_.size();++i)result.col_ptr_[i]+=result.col_ptr_[i-1];
            result.rows_and_coeffs_.resize(rows_and_coeffs_.size());
            auto next=result.col_ptr_;
            for(std::size_t c=0;c<domain_dim();++c)
                for(const auto& term:get_column(c))
                    result.rows_and_coeffs_[next[term.getValue()]++]=Basis(c,term.getCoefficient());
            return result;
        }

		DenseDomainVec reserve_dense_domain_vec() const {
			return std::make_unique<k[]>(static_cast<std::size_t>(domain_dim()));
		}

		DenseImageVec reserve_dense_image_vec() const {
			return std::make_unique<k[]>(static_cast<std::size_t>(image_dim_));
		}


		// x = M^T * y
		DenseDomainVec evaluate_transpose(const DenseImageVec& input) const {
			DenseDomainVec result = reserve_dense_domain_vec();
			evaluate_transpose(input, result);
			return result;
		}


		void evaluate_transpose(const DenseImageVec& input, DenseDomainVec& result) const {
			const indexType ncols = domain_dim();

#pragma omp parallel for schedule(static)
			for (std::int64_t ci = 0; ci < static_cast<std::int64_t>(ncols); ++ci) {
				const indexType col = static_cast<indexType>(ci);

				k acc = k{0};
				for (const auto& be : get_column(col)) {
					acc += input[static_cast<std::size_t>(be.getValue())] * be.getCoefficient();
				}
				result[static_cast<std::size_t>(col)] = acc;
			}
		}


};
