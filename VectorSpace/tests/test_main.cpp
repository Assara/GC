#include "../wiedemann_helper.hpp"
#include "../block_wiedemann.hpp"
#include "../types.hpp"
#include <iostream>
#include <random>
#include <cstdint>
#include <chrono>
#include <cstdlib>
#include <string>

using namespace std;
using namespace VectorSpace;


lil_matrix<fieldType> create_diagonal_id_matrix(size_t dim, fieldType lambda) {
	lil_matrix<fieldType> M;


	for (size_t i = 0; i < dim ; i++) {

		LinComb<size_t, fieldType> tmp_col;
		tmp_col.append_in_basis_order(i, lambda);
		M.add_col(std::move(tmp_col)); 
	}

	return M;
}


lil_matrix<fieldType> create_triangular(size_t dim, fieldType lambda) {
	lil_matrix<fieldType> M;

	for (size_t i = 0; i < dim ; i++) {

		LinComb<size_t, fieldType> tmp_col;

		for (size_t j =0; j <= i; j++) {	
			tmp_col.append_in_basis_order(j, lambda);
		}

		M.add_col(std::move(tmp_col));
	}

	return M;	
}


bool test_wiedemann_solver_id() {
	int size = 17;

	std::unique_ptr<fieldType[]> vec = std::make_unique<fieldType[]>(size);

	for (int i = 0; i < size; i++) {
		vec[i] = fieldType(i);
	}

	lil_matrix<fieldType> id = create_diagonal_id_matrix(size, fieldType{1});


	wiedemann_solver<fieldType> solver(id);

	auto X = solver.solve_MX_equals_y(id.make_dense_image_vec_random());


	if (X.has_value()) {
		std::cout << "Wiedemann solver works for id matrix" << std::endl;
	}
	return X.has_value();
}


bool test_wiedemann_solver_triangular() {
	int size = 3;

	std::unique_ptr<fieldType[]> vec = std::make_unique<fieldType[]>(size);

	for (int i = 0; i < size; i++) {
		vec[i] = fieldType(i);
	}

	lil_matrix<fieldType> id = create_triangular(size, fieldType{2});


	wiedemann_solver<fieldType> solver(id);

	auto X = solver.solve_MX_equals_y(id.make_dense_image_vec_random());


	if (X.has_value()) {
		std::cout << "Wiedemann solver works for triangular matrix" << std::endl;
	}
	return X.has_value();
}

// k must have method modulus()
template<typename k>
lil_matrix<k> generate_random_matrix(size_t image_dim, size_t domain_dim, size_t n_entries, std::uint64_t seed = 17) {
	std::mt19937_64 rng(seed);

	lil_matrix<k> M;

	if (image_dim == 0 || domain_dim == 0) {
		return M;
	}

	std::uniform_int_distribution<std::uint64_t>
		coeff_dist(0, k::modulus() - 1);

	std::uniform_int_distribution<std::size_t>
		image_dist(0, image_dim - 1),
		domain_dist(0, domain_dim - 1);

	// Pin the requested dimensions so timings do not drift with random sparsity.
	for (std::size_t j = 0; j < domain_dim; ++j) {
		k value;
		do {
			value = k{coeff_dist(rng)};
		} while (value == k{0});

		M.add_element(image_dist(rng), j, value);
	}
	M.add_element(image_dim - 1, 0, k{1});

	for (std::size_t t = domain_dim; t < n_entries; ++t) {
		std::size_t i = image_dist(rng);
		std::size_t j = domain_dist(rng);

		k value;
		do {
			value = k{coeff_dist(rng)};
		} while (value == k{0});  // optional but recommended

		M.add_element(i, j, value);
	}

	M.sort_cols();

	return M;
}

template<typename k>
bool dense_vectors_equal(const std::unique_ptr<k[]>& lhs,
		const std::unique_ptr<k[]>& rhs,
		std::size_t size) {
	for (std::size_t i = 0; i < size; ++i) {
		if (lhs[i] != rhs[i]) {
			return false;
		}
	}
	return true;
}



bool performance_test_wiedemann_solver(size_t image_dim, size_t domain_dim, size_t n_entries, std::uint64_t seed = 17) {

	lil_matrix<fieldType> M = generate_random_matrix<fieldType>(image_dim, domain_dim, n_entries, seed);

	std::unique_ptr<fieldType[]> pre_x = M.make_dense_domain_vec_random(seed + 1);
	std::unique_ptr<fieldType[]> y = M.evaluate_from_dense(pre_x);



	wiedemann_solver<fieldType> solver(M);

	std::cout << "Wiedemann scalar performance test: image_dim="
		<< M.image_dim()
		<< " domain_dim=" << M.domain_dim()
		<< " nnz=" << M.size()
		<< " seed=" << seed
		<< std::endl;

	using clock = std::chrono::steady_clock;

	const auto t0 = clock::now();
	std::optional<std::unique_ptr<fieldType[]>> X = solver.solve_MX_equals_y(y);
	const auto t1 = clock::now();


	const std::chrono::duration<double> elapsed = t1 - t0;

	std::cout << "Wiedemann solve time: "
		<< elapsed.count() << " s\n";


	if (X.has_value()) {	
		std::unique_ptr<fieldType[]> y_found = M.evaluate_from_dense(*X);
		if (dense_vectors_equal(y, y_found, M.image_dim())) {
			std::cout << "solved and verified" << std::endl;
			return true;
		}
		std::cout << "solution failed verification" << std::endl;
		return false;
	}

	cout << "solution not found " << std::endl;


	return false;
}


bool performance_test_block_wiedemann_solver(size_t image_dim, size_t domain_dim, size_t n_entries, std::uint64_t seed = 17) {

	lil_matrix<fieldType> M = generate_random_matrix<fieldType>(image_dim, domain_dim, n_entries, seed);

	std::unique_ptr<fieldType[]> pre_x = M.make_dense_domain_vec_random(seed + 1);
	std::unique_ptr<fieldType[]> y = M.evaluate_from_dense(pre_x);

	VectorSpace::block_wiedemann_solver<fieldType> solver(M, 16);


	using clock = std::chrono::steady_clock;

	const auto t0 = clock::now();
	std::vector<std::unique_ptr<fieldType[]>> X = solver.solve_MX_equals_y(y);
	const auto t1 = clock::now();


	const std::chrono::duration<double> elapsed = t1 - t0;

	std::cout << "Wiedemann solve time: "
		<< elapsed.count() << " s\n";

	cout << "found " << X.size() << "solutions" << std::endl;

	return true;
}



std::size_t parse_size_arg(char* arg) {
	return std::strtoull(arg, nullptr, 10);
}

int main(int argc, char** argv) {
	if (!test_wiedemann_solver_id()) return 1;
	if (!test_wiedemann_solver_triangular()) return 1;

	if (argc > 1 && std::string(argv[1]) == "--block") {
		const std::size_t image_dim = argc > 2 ? parse_size_arg(argv[2]) : 1000;
		const std::size_t domain_dim = argc > 3 ? parse_size_arg(argv[3]) : 1000;
		const std::size_t n_entries = argc > 4 ? parse_size_arg(argv[4]) : 8000;
		const std::uint64_t seed = argc > 5 ? parse_size_arg(argv[5]) : 17;
		return performance_test_block_wiedemann_solver(image_dim, domain_dim, n_entries, seed) ? 0 : 1;
	}

	const std::size_t image_dim = argc > 1 ? parse_size_arg(argv[1]) : 18557;
	const std::size_t domain_dim = argc > 2 ? parse_size_arg(argv[2]) : 9846;
	const std::size_t n_entries = argc > 3 ? parse_size_arg(argv[3]) : 65017;
	const std::uint64_t seed = argc > 4 ? parse_size_arg(argv[4]) : 17;

	return performance_test_wiedemann_solver(image_dim, domain_dim, n_entries, seed) ? 0 : 1;
}
