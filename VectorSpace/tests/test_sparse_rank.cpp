#include "VectorSpace/sparse_rank.hpp"
#include "types.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace {

using Field = fieldType;
using Matrix = compressed_sparse_matrix<Field>;
using Rank = VectorSpace::sparse_rank<Field>;
using DenseMatrix = std::vector<std::vector<Field>>;

bool expect(bool condition, std::string_view message) {
	if (!condition) {
		std::cerr << "FAIL: " << message << '\n';
	}
	return condition;
}

Matrix compressed_from_dense(const DenseMatrix& values, std::size_t columns) {
	const std::size_t rows = values.size();
	Matrix matrix(static_cast<Matrix::indexType>(rows));
	for (std::size_t column = 0; column < columns; ++column) {
		std::vector<Matrix::Basis> sparse_column;
		for (std::size_t row = 0; row < rows; ++row) {
			if (values[row][column] != Field{}) {
				sparse_column.emplace_back(
					static_cast<Matrix::indexType>(row),
					values[row][column]
				);
			}
		}
		matrix.add_col(sparse_column);
	}
	return matrix;
}

std::size_t dense_rank(DenseMatrix matrix, std::size_t columns) {
	const std::size_t rows = matrix.size();
	std::size_t pivot_row = 0;
	for (std::size_t column = 0;
		column < columns && pivot_row < rows;
		++column) {
		std::size_t pivot = pivot_row;
		while (pivot < rows && matrix[pivot][column] == Field{}) {
			++pivot;
		}
		if (pivot == rows) {
			continue;
		}
		std::swap(matrix[pivot], matrix[pivot_row]);

		const Field inverse = Field{1} / matrix[pivot_row][column];
		for (std::size_t j = column; j < columns; ++j) {
			matrix[pivot_row][j] *= inverse;
		}
		for (std::size_t row = pivot_row + 1; row < rows; ++row) {
			const Field factor = matrix[row][column];
			if (factor == Field{}) {
				continue;
			}
			for (std::size_t j = column; j < columns; ++j) {
				matrix[row][j] -= factor * matrix[pivot_row][j];
			}
		}
		++pivot_row;
	}
	return pivot_row;
}

bool test_empty_and_zero_matrices() {
	Matrix empty(0);
	const Rank empty_rank(empty);
	bool ok = expect(empty_rank.result() == Rank::result_type{0, 0},
		"0 x 0 matrix has rank and nullity zero");

	Matrix zero(5);
	for (int column = 0; column < 3; ++column) {
		zero.add_col({});
	}
	const Rank zero_rank(zero);
	ok &= expect(zero_rank.result() == Rank::result_type{0, 3},
		"five-by-three zero matrix has nullity three");
	return ok;
}

bool test_rectangular_matrix() {
	// c2 = 3*c0 - 2*c1, c3 is zero, and c4 adds the third pivot.
	const DenseMatrix dense{
		{Field{1}, Field{0}, Field{3}, Field{0}, Field{0}},
		{Field{0}, Field{1}, Field{-2}, Field{0}, Field{0}},
		{Field{2}, Field{0}, Field{6}, Field{0}, Field{1}},
	};
	const Matrix matrix = compressed_from_dense(dense, 5);
	const Rank rank(matrix);
	return expect(rank.rank() == 3, "rectangular matrix has rank three")
		&& expect(rank.nullity() == 2, "rectangular matrix has nullity two");
}

bool test_unsorted_duplicates_and_explicit_zeros() {
	Matrix matrix(3);
	using Basis = Matrix::Basis;

	// The first column normalizes to (1, 0, 2)^T.
	matrix.add_col({
		Basis{2, Field{5}},
		Basis{0, Field{1}},
		Basis{2, Field{-3}},
		Basis{1, Field{0}},
	});
	// Duplicate entries cancel, so this is the zero column.
	matrix.add_col({
		Basis{1, Field{7}},
		Basis{0, Field{-1}},
		Basis{1, Field{-7}},
		Basis{0, Field{1}},
	});
	// Independent of the first column.
	matrix.add_col({Basis{1, Field{4}}});

	const Rank rank(matrix);
	return expect(rank.result() == Rank::result_type{2, 1},
		"duplicates are combined and explicit zeros are ignored");
}

bool test_invalid_row_is_rejected() {
	Matrix matrix(2);
	matrix.add_col({Matrix::Basis{2, Field{1}}});
	try {
		[[maybe_unused]] const Rank rank(matrix);
	} catch (const std::out_of_range&) {
		return true;
	}
	return expect(false, "out-of-range row throws std::out_of_range");
}

bool test_invalid_compressed_structure_is_rejected() {
	auto throws_invalid_argument = [](const Matrix& matrix) {
		try {
			[[maybe_unused]] const Rank rank(matrix);
		} catch (const std::invalid_argument&) {
			return true;
		}
		return false;
	};

	Matrix missing_sentinel(2);
	missing_sentinel.col_ptr_.clear();
	bool ok = expect(throws_invalid_argument(missing_sentinel),
		"missing CSC sentinel is rejected");

	Matrix descending_offsets(2);
	descending_offsets.add_col({Matrix::Basis{0, Field{1}}});
	descending_offsets.add_col({});
	descending_offsets.col_ptr_[1] = 1;
	descending_offsets.col_ptr_[2] = 0;
	ok &= expect(throws_invalid_argument(descending_offsets),
		"descending CSC offsets are rejected");

	Matrix trailing_entries(2);
	trailing_entries.add_col({Matrix::Basis{0, Field{1}}});
	trailing_entries.col_ptr_.back() = 0;
	ok &= expect(throws_invalid_argument(trailing_entries),
		"entries outside every CSC column are rejected");
	return ok;
}

bool test_random_matrices_against_dense_elimination() {
	std::mt19937_64 rng(0x51A9EULL);
	std::uniform_int_distribution<int> dimension(0, 7);
	std::uniform_int_distribution<int> coefficient(-3, 3);

	for (int sample = 0; sample < 300; ++sample) {
		const std::size_t rows = static_cast<std::size_t>(dimension(rng));
		const std::size_t columns = static_cast<std::size_t>(dimension(rng));
		DenseMatrix dense(rows, std::vector<Field>(columns));
		for (std::size_t row = 0; row < rows; ++row) {
			for (std::size_t column = 0; column < columns; ++column) {
				dense[row][column] = Field{coefficient(rng)};
			}
		}

		const Matrix matrix = compressed_from_dense(dense, columns);
		const std::size_t expected_rank = dense_rank(dense, columns);
		const Rank actual(matrix);
		if (actual.rank() != expected_rank
			|| actual.nullity() != columns - expected_rank) {
			std::cerr << "FAIL: random matrix sample " << sample
				<< " (" << rows << " x " << columns << ") expected rank "
				<< expected_rank << " but got " << actual.rank() << '\n';
			return false;
		}
	}
	return true;
}

} // namespace

int main() {
	bool ok = true;
	ok &= test_empty_and_zero_matrices();
	ok &= test_rectangular_matrix();
	ok &= test_unsorted_duplicates_and_explicit_zeros();
	ok &= test_invalid_row_is_rejected();
	ok &= test_invalid_compressed_structure_is_rejected();
	ok &= test_random_matrices_against_dense_elimination();

	if (!ok) {
		return 1;
	}
	std::cout << "sparse rank tests passed\n";
	return 0;
}
