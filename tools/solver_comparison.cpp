#include <chrono>
#include <cstdlib>
#include <iostream>
#include <optional>
#include <random>
#include <string>
#include <tuple>
#include <vector>

#if defined(GC_PERF_WITH_LINBOX)
#if __has_include(<givaro/modular.h>) && __has_include(<linbox/matrix/sparse-matrix.h>) && __has_include(<linbox/solutions/solve.h>)
#include <givaro/modular.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/solve.h>
#include <linbox/vector/vector.h>
#else
#error "GC_PERF_WITH_LINBOX was set, but LinBox/Givaro headers were not found"
#endif
#endif

#include "types.hpp"
#include "VectorSpace/lil_matrix.hpp"
#include "VectorSpace/wiedemann_helper.hpp"

namespace {

struct MatrixEntry {
	std::size_t row;
	std::size_t col;
	fieldType value;
};

template <typename Fn>
double time_seconds(Fn&& fn) {
	const auto start = std::chrono::steady_clock::now();
	fn();
	const auto stop = std::chrono::steady_clock::now();
	return std::chrono::duration<double>(stop - start).count();
}

std::vector<MatrixEntry> random_entries(
	std::size_t image_dim,
	std::size_t domain_dim,
	std::size_t n_entries,
	std::uint64_t seed
) {
	std::mt19937_64 rng(seed);
	std::uniform_int_distribution<std::size_t> row_dist(0, image_dim - 1);
	std::uniform_int_distribution<std::size_t> col_dist(0, domain_dim - 1);
	std::uniform_int_distribution<std::uint64_t> coeff_dist(1, fieldType::modulus() - 1);

	std::vector<MatrixEntry> entries;
	entries.reserve(n_entries);
	for (std::size_t i = 0; i < n_entries; ++i) {
		entries.push_back({
			row_dist(rng),
			col_dist(rng),
			fieldType{coeff_dist(rng)}
		});
	}
	return entries;
}

lil_matrix<fieldType> make_matrix(
	std::size_t image_dim,
	std::size_t domain_dim,
	const std::vector<MatrixEntry>& entries
) {
	lil_matrix<fieldType> matrix(domain_dim);
	for (const auto& entry : entries) {
		matrix.add_element(entry.row, entry.col, entry.value);
	}
	if (image_dim > 0 && domain_dim > 0) {
		matrix.add_element(image_dim - 1, domain_dim - 1, fieldType{0});
	}
	return matrix;
}

std::unique_ptr<fieldType[]> random_domain_vector(std::size_t domain_dim, std::uint64_t seed) {
	std::mt19937_64 rng(seed);
	std::uniform_int_distribution<std::uint64_t> coeff_dist(0, fieldType::modulus() - 1);
	auto result = std::make_unique<fieldType[]>(domain_dim);
	for (std::size_t i = 0; i < domain_dim; ++i) {
		result[i] = fieldType{coeff_dist(rng)};
	}
	return result;
}

bool vectors_match(
	const std::unique_ptr<fieldType[]>& lhs,
	const std::unique_ptr<fieldType[]>& rhs,
	std::size_t size
) {
	for (std::size_t i = 0; i < size; ++i) {
		if (lhs[i] != rhs[i]) {
			return false;
		}
	}
	return true;
}

int run_solver_comparison(
	std::size_t image_dim,
	std::size_t domain_dim,
	std::size_t n_entries,
	std::uint64_t seed
) {
	const auto entries = random_entries(image_dim, domain_dim, n_entries, seed);
	auto matrix = make_matrix(image_dim, domain_dim, entries);
	const auto x_expected = random_domain_vector(domain_dim, seed + 1);
	const auto y = matrix.evaluate_from_dense(x_expected);

	std::cout << "solver comparison\n";
	std::cout << "image_dim = " << image_dim << '\n';
	std::cout << "domain_dim = " << domain_dim << '\n';
	std::cout << "requested entries = " << n_entries << '\n';
	std::cout << "stored entries = " << matrix.size() << '\n';

	std::optional<std::unique_ptr<fieldType[]>> own_solution;
	const double own_seconds = time_seconds([&]() {
		VectorSpace::wiedemann_solver<fieldType> solver(matrix);
		own_solution = solver.solve_MX_equals_y(y);
	});
	std::cout << "own wiedemann time = " << own_seconds << " s\n";
	std::cout << "own solution found = " << (own_solution.has_value() ? "yes" : "no") << '\n';
	if (own_solution.has_value()) {
		const auto y_check = matrix.evaluate_from_dense(*own_solution);
		std::cout << "own solution validates = "
		          << (vectors_match(y, y_check, image_dim) ? "yes" : "no") << '\n';
	}

#if defined(GC_PERF_WITH_LINBOX)
	using Field = Givaro::Modular<double>;
	using Matrix = LinBox::SparseMatrix<Field>;
	using Vector = LinBox::DenseVector<Field>;

	Field field(fieldType::modulus());
	std::cout << "linbox field = GF(" << fieldType::modulus() << ")\n";
	Matrix linbox_matrix(field, image_dim, domain_dim);
	for (const auto& entry : entries) {
		Field::Element value;
		field.init(value, entry.value.value());
		linbox_matrix.setEntry(entry.row, entry.col, value);
	}

	Vector linbox_rhs(field, image_dim);
	for (std::size_t i = 0; i < image_dim; ++i) {
		field.init(linbox_rhs[i], y[i].value());
	}
	Vector linbox_solution(field, domain_dim);

	const double linbox_seconds = time_seconds([&]() {
		LinBox::solve(
			linbox_solution,
			linbox_matrix,
			linbox_rhs,
			LinBox::RingCategories::ModularTag(),
			LinBox::Method::Wiedemann()
		);
	});
	std::cout << "linbox wiedemann time = " << linbox_seconds << " s\n";
#else
	std::cout << "linbox wiedemann time = disabled; rebuild this tool with GC_PERF_WITH_LINBOX\n";
#endif

	return EXIT_SUCCESS;
}

} // namespace

int main(int argc, char** argv) {
	const std::size_t image_dim = argc >= 2 ? std::stoull(argv[1]) : 5000;
	const std::size_t domain_dim = argc >= 3 ? std::stoull(argv[2]) : 5000;
	const std::size_t n_entries = argc >= 4 ? std::stoull(argv[3]) : 50000;
	const std::uint64_t seed = argc >= 5 ? std::stoull(argv[4]) : 17;

	if (image_dim == 0 || domain_dim == 0) {
		std::cerr << "usage: " << argv[0] << " [image_dim] [domain_dim] [n_entries] [seed]\n";
		return EXIT_FAILURE;
	}

	return run_solver_comparison(image_dim, domain_dim, n_entries, seed);
}
