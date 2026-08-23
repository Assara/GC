#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "compressed_sparse_matrix.hpp"

namespace VectorSpace {

// Exact rank and kernel dimension of a column-compressed sparse matrix.
//
// compressed_sparse_matrix stores the columns of an image_dim x domain_dim
// matrix.  We incrementally reduce those columns against a sparse echelon
// basis, always using the least row index as the pivot.  The input need not
// have sorted columns: duplicate row entries are combined and explicit zeros
// are discarded before elimination.
//
// This is intended as a deterministic correctness baseline.  It never makes
// the whole matrix dense, but Gaussian elimination can still create fill-in;
// very large graph-complex matrices will eventually need pivot heuristics or a
// black-box rank algorithm on top of this interface.
template<typename Field>
class sparse_rank {
public:
	using field_type = Field;
	using matrix_type = compressed_sparse_matrix<Field>;
	using index_type = typename matrix_type::indexType;

	struct result_type {
		std::size_t rank{0};
		std::size_t nullity{0};

		friend bool operator==(const result_type&, const result_type&) = default;
	};

	explicit sparse_rank(const matrix_type& matrix)
		: result_(compute(matrix)) {}

	std::size_t rank() const noexcept {
		return result_.rank;
	}

	// Dimension of ker(M), hence domain_dim(M) - rank(M).
	std::size_t nullity() const noexcept {
		return result_.nullity;
	}

	result_type result() const noexcept {
		return result_;
	}

	static result_type compute(const matrix_type& matrix) {
		using pivot_slot_type = index_type;
		constexpr pivot_slot_type no_pivot =
			std::numeric_limits<pivot_slot_type>::max();

		const std::size_t image_dim =
			static_cast<std::size_t>(matrix.image_dim());
		const std::size_t domain_dim = validated_domain_dim(matrix);

		// pivot_slot[row] identifies the normalized echelon column whose
		// first entry is in that row.  A 32-bit slot is sufficient because
		// compressed_sparse_matrix also uses 32-bit matrix dimensions.
		std::vector<pivot_slot_type> pivot_slot(image_dim, no_pivot);
		std::vector<sparse_column> echelon_columns;
		echelon_columns.reserve(std::min(image_dim, domain_dim));

		for (std::size_t column_index = 0;
			column_index < domain_dim;
			++column_index) {
			sparse_column column = normalized_column(
				matrix,
				static_cast<index_type>(column_index)
			);

			while (!column.empty()) {
				const index_type pivot_row = column.front().row;
				const pivot_slot_type slot =
					pivot_slot[static_cast<std::size_t>(pivot_row)];

				if (slot == no_pivot) {
					normalize_pivot(column);
					pivot_slot[static_cast<std::size_t>(pivot_row)] =
						static_cast<pivot_slot_type>(echelon_columns.size());
					echelon_columns.push_back(std::move(column));
					break;
				}

				// Stored pivot columns have leading coefficient one.
				const Field factor = -column.front().coefficient;
				column = add_scaled(
					column,
					echelon_columns[static_cast<std::size_t>(slot)],
					factor
				);
			}
		}

		const std::size_t rank = echelon_columns.size();
		return result_type{rank, domain_dim - rank};
	}

private:
	struct entry {
		index_type row{};
		Field coefficient{};
	};

	using sparse_column = std::vector<entry>;

	static std::size_t validated_domain_dim(const matrix_type& matrix) {
		if (matrix.col_ptr_.empty()) {
			throw std::invalid_argument(
				"compressed sparse matrix is missing its initial column offset"
			);
		}
		if (matrix.col_ptr_.front() != typename matrix_type::offset_type{0}) {
			throw std::invalid_argument(
				"compressed sparse matrix must begin at offset zero"
			);
		}

		const std::size_t domain_dim = matrix.col_ptr_.size() - 1;
		if (domain_dim
			> static_cast<std::size_t>(std::numeric_limits<index_type>::max())) {
			throw std::length_error(
				"compressed sparse matrix domain exceeds its index type"
			);
		}

		typename matrix_type::offset_type previous = 0;
		for (const typename matrix_type::offset_type offset : matrix.col_ptr_) {
			if (offset < previous || offset > matrix.rows_and_coeffs_.size()) {
				throw std::invalid_argument(
					"compressed sparse matrix has invalid column offsets"
				);
			}
			previous = offset;
		}
		if (matrix.col_ptr_.back() != matrix.rows_and_coeffs_.size()) {
			throw std::invalid_argument(
				"compressed sparse matrix has unassigned stored entries"
			);
		}
		return domain_dim;
	}

	static sparse_column normalized_column(
		const matrix_type& matrix,
		index_type column_index
	) {
		sparse_column entries;
		const auto input = matrix.get_column(column_index);
		if (input.begin() != input.end()) {
			entries.reserve(
				static_cast<std::size_t>(input.end() - input.begin())
			);
		}

		for (const auto& term : input) {
			const index_type row = term.getValue();
			if (row >= matrix.image_dim()) {
				throw std::out_of_range(
					"compressed sparse matrix contains an out-of-range row"
				);
			}
			const Field coefficient = term.getCoefficient();
			if (coefficient != Field{}) {
				entries.push_back(entry{row, coefficient});
			}
		}

		std::sort(
			entries.begin(),
			entries.end(),
			[](const entry& lhs, const entry& rhs) {
				return lhs.row < rhs.row;
			}
		);

		sparse_column result;
		result.reserve(entries.size());
		std::size_t begin = 0;
		while (begin < entries.size()) {
			const index_type row = entries[begin].row;
			Field coefficient{};
			do {
				coefficient += entries[begin].coefficient;
				++begin;
			} while (begin < entries.size() && entries[begin].row == row);

			if (coefficient != Field{}) {
				result.push_back(entry{row, coefficient});
			}
		}
		return result;
	}

	static void normalize_pivot(sparse_column& column) {
		const Field inverse = Field{1} / column.front().coefficient;
		for (entry& term : column) {
			term.coefficient *= inverse;
		}
	}

	// lhs + scalar * rhs.  Both operands are sorted and contain no zero or
	// duplicate entries; the result preserves those invariants.
	static sparse_column add_scaled(
		const sparse_column& lhs,
		const sparse_column& rhs,
		const Field& scalar
	) {
		if (scalar == Field{}) {
			return lhs;
		}

		sparse_column result;
		result.reserve(lhs.size() + rhs.size());

		std::size_t left = 0;
		std::size_t right = 0;
		while (left < lhs.size() && right < rhs.size()) {
			if (lhs[left].row < rhs[right].row) {
				result.push_back(lhs[left++]);
				continue;
			}
			if (rhs[right].row < lhs[left].row) {
				const Field coefficient = rhs[right].coefficient * scalar;
				if (coefficient != Field{}) {
					result.push_back(entry{rhs[right].row, coefficient});
				}
				++right;
				continue;
			}

			const Field coefficient =
				lhs[left].coefficient + rhs[right].coefficient * scalar;
			if (coefficient != Field{}) {
				result.push_back(entry{lhs[left].row, coefficient});
			}
			++left;
			++right;
		}

		result.insert(result.end(), lhs.begin() + left, lhs.end());
		for (; right < rhs.size(); ++right) {
			const Field coefficient = rhs[right].coefficient * scalar;
			if (coefficient != Field{}) {
				result.push_back(entry{rhs[right].row, coefficient});
			}
		}
		return result;
	}

	result_type result_{};
};

} // namespace VectorSpace
