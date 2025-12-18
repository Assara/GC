
#pragma once

#include "LinComb.hpp"
#include "wiedemann_helper.hpp"


#include <vector>
#include <unordered_map>
#include <optional>
#include <memory>
#include <algorithm>
#include <iostream>

namespace VectorSpace {

template<typename A, typename B, typename k>
class wiedemann_primitive_finder {


private:
    std::vector<B> domain_space_enumeration;               // row index -> B basis element  
    lil_matrix<k> map_representative;
                           // map B -> A
    std::unordered_map<A, std::size_t> image_space_enumeration; // A basis -> column index
 

    using DenseDomainVec = std::unique_ptr<k[]>; //Assumed to have the size of the map_representative matrix
    using DenseImageVec = std::unique_ptr<k[]>; 
    
  
	std::size_t image_space_to_number(const A& a) {
		auto [it, inserted] =
			image_space_enumeration.try_emplace(a, image_space_enumeration.size());

		return it->second;
	}
	
	LinComb<size_t, k> map_to_enumeration_basis(const LinComb<A,k>& input) {
			LinComb<size_t, k> result;
			
			for (const auto& be : input ) {
				result.append_in_basis_order(BasisElement(image_space_to_number(be.getValue()), be.getCoefficient()));
			}
			
			
			result.sort_without_deduplicate();
			return result;
	}
	
	
	DenseImageVec map_to_enumeration_basis_dense(const LinComb<A,k>& input) {
			DenseImageVec result = map_representative.make_dense_image_vec_zero();
			
			for (const auto& be : input ) {
				result[image_space_enumeration[be.getValue()]] += be.getCoefficient();	
				result[image_space_enumeration[be.getValue()]] += be.getCoefficient();	
			}
			return result;
	}
	
	
	LinComb<B, k> domain_enumeration_inverse_from_dense(const DenseDomainVec& input) {
			LinComb<B, k> result;
			
			
			for (size_t i = 0; i< map_representative.domain_dim(); i++) {
				result.append_in_basis_order(BasisElement(domain_space_enumeration[i], input[i]));
			}
				
			//result.sort_without_deduplicate(); //this does not have to be sorted for current applications
			return result;
	}
	

public:

    // -----------------------------------------------------
    // Constructor: build sparse matrix for δ : B → LinComb<A,k>
    // -----------------------------------------------------
    wiedemann_primitive_finder(std::unordered_map<B, LinComb<A,k>>& delta) {
        domain_space_enumeration.reserve(delta.size());
        size_t n_matrix_entries = 0;
   
        for (const auto& entry : delta) {
	
            domain_space_enumeration.push_back(entry.first);   // each B becomes a row
			
			auto col = map_to_enumeration_basis(entry.second);
			
			n_matrix_entries += col.size();
			map_representative.add_col(std::move(col));	
			
        }
  
		cout << "created sparse solver: domain_dim = " << domain_space_enumeration.size() << endl
				<< "image_space_dim = " <<  image_space_enumeration.size() << endl
				<< "number of matrix entries =" <<  n_matrix_entries << endl;
	
    }

	std::optional<LinComb<B,k>> find_primitive_or_empty(LinComb<A,k> y) {
		cout << "using sparse_primitive_finder for LinComb:" << endl;
		
		auto y_enumerated = map_to_enumeration_basis_dense(y);
		
		cout << "mapped to enumeration!!!" << endl;
		
		//we do not need this anymore, but we do need heap space.
		image_space_enumeration.clear();
		wiedemann_solver<k> solver(map_representative);

		std::optional<DenseImageVec> X_enumerated = solver.solve_MX_equals_y(y_enumerated);
		
		cout << "solved!!!" << endl;
		
		if (!X_enumerated.has_value()) {
				cout << "no solution!" << endl;
				
				return std::nullopt;
		}
	
		LinComb<B,k> X = domain_enumeration_inverse_from_dense(*X_enumerated);
		
		return std::optional(std::move(X));
	}

};

} // namespace VectorSpace
