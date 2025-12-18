
#include "../wiedemann_helper.hpp"
#include "../types.hpp"
#include <iostream>

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
		
		lil_matrix<fieldType> id = create_diagonal_id_matrix(size, 1);
		
		
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
		
		lil_matrix<fieldType> id = create_triangular(size, 2);
		
		
		wiedemann_solver<fieldType> solver(id);
		
		auto X = solver.solve_MX_equals_y(id.make_dense_image_vec_random());
		
		
		if (X.has_value()) {
			std::cout << "Wiedemann solver works for triangular matrix" << std::endl;
		}
		return X.has_value();
}


int main() {
		test_wiedemann_solver_id();
		
		test_wiedemann_solver_triangular();
}

