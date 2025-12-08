#pragma once

#include "LinComb.hpp"
#include "lil_matrix.hpp"

#include <vector>
#include <unordered_map>
#include <optional>
#include <memory>
#include <algorithm>
#include <iostream>

namespace VectorSpace {

template<typename A, typename B, typename k>
class sparse_primitive_finder {



private:
    std::vector<B> domain_space_enumeration;               // row index -> B basis element
    lil_matrix<k> map_representative;                           // map B -> A
    std::unordered_map<A, std::size_t> image_space_enumeration; // A basis -> column index
    
  
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
	
	LinComb<B, k> domain_enumeration_inverse(const LinComb<size_t, k>& input) {
			LinComb<B, k> result;
			
			for (const auto& be : input) {
				if(be.getValue() >= domain_space_enumeration.size()) {
						cout << "BAD VALUE IN domain_enumeration_inverse: " << be.getValue()<< " max = " << domain_space_enumeration.size() << endl;
						
						//domain_space_enumeration[be.getValue()].print();
				} 
				
				result.append_in_basis_order(BasisElement(domain_space_enumeration[be.getValue()], be.getCoefficient()));
			}			
			//result.sort_without_deduplicate(); //this does not have to be sorted for current applications
			return result;
	}
	

public:

    // -----------------------------------------------------
    // Constructor: build sparse matrix for δ : B → LinComb<A,k>
    // -----------------------------------------------------
    sparse_primitive_finder(std::unordered_map<B, LinComb<A,k>>& delta) {
        domain_space_enumeration.reserve(delta.size());
        
        for (const auto& entry : delta) {
	
            domain_space_enumeration.push_back(entry.first);   // each B becomes a row
  
			map_representative.add_col(map_to_enumeration_basis(entry.second));
			
        }
  
  
		cout << "finished contructing sparse matrix and enumeration maps " << endl; 
    }


	std::optional<LinComb<B,k>> find_primitive_or_empty(LinComb<A,k> y) {
		cout << "using sparse_primitive_finder for LinComb:" << endl;
		
		
		auto y_enumerated = map_to_enumeration_basis(y);
		
		cout << "mapped to enumeration!!!" << endl;
		
		
		auto X_enumerated = map_representative.solve_MX_equals_y(y_enumerated);
		
		
		cout << "solved!!!" << endl;
		
		if (!X_enumerated.has_value()) {
				cout << "no solution!" << endl;
				
				return std::nullopt;
		}
		
		
		LinComb<B,k> X = domain_enumeration_inverse(*X_enumerated);
		
		return std::optional(std::move(X));
	}

};

} // namespace VectorSpace
