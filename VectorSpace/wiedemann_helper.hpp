#pragma once

#include <vector>
#include <optional>
#include <random>
#include <cstddef>
#include <utility>
#include <iostream>

#include "lil_matrix.hpp"
#include "berlekamp_massey_state.hpp"

namespace VectorSpace {

template<typename k>
class wiedemann_solver {
public:
    using ImageVec  = typename lil_matrix<k>::DenseImageVec;   // expected: std::unique_ptr<k[]>
    using DomainVec = typename lil_matrix<k>::DenseDomainVec;  // expected: std::unique_ptr<k[]>

private:
    const lil_matrix<k>& M;
    ImageVec random_vector; // length = M.image_dim()

private:
    ImageVec M_MT(const ImageVec& y) const {
        // Requires lil_matrix to accept const refs for dense vectors.
        DomainVec v = M.evaluate_transpose_dense(y);
        return M.evaluate_from_dense(v);
    }

    DomainVec MT_M(const DomainVec& x) const {
        ImageVec v = M.evaluate_from_dense(x);
        return M.evaluate_transpose_dense(v);
    }

    k get_signature(const ImageVec& y) const {
        k result = k{0};
        for (std::size_t i = 0; i < M.image_dim(); ++i) {
            result += y[i] * random_vector[i];
        }
        return result;
    }

    static void add_scaled_inplace(std::unique_ptr<k[]>& dst,
                                   const std::unique_ptr<k[]>& src,
                                   const std::size_t size,
                                   const k scalar) {
        for (std::size_t i = 0; i < size; ++i)
            dst[i] += src[i] * scalar;
    }
    
    
     static std::unique_ptr<k[]> scaled_copy(const std::unique_ptr<k[]>& vec,
                                   const std::size_t size,
                                   const k scalar) {
									   
		std::unique_ptr<k[]> result = std::make_unique<k[]>(size);						   
									
        for (std::size_t i = 0; i < size; ++i) {
            result[i] = scalar * vec[i];
        }
        
        return result;
    }
    
    std::optional<DomainVec>
recompute_solution(const ImageVec& y0, const std::vector<k>& connection_poly) const {

    std::cout << "recomputing solution connection_poly.size() = "
              << connection_poly.size() << std::endl;

    // BM returns: [1, c1, ..., cL]
    if (connection_poly.size() < 2 || connection_poly[0] == k{0})
        return std::nullopt;

		const std::size_t L = connection_poly.size() - 1;

		// Constant term of operator polynomial f(t) = t^L + c1 t^{L-1} + ... + cL
		const k cL = connection_poly.back();
		if (cL == k{0})
			return std::nullopt; // (your current policy)

		const k scale = -(cL.inv()); // -1/cL

		// Build v[j] = A^j y0 for j=0..L-1, A = M M^T
		std::vector<ImageVec> v;
		v.reserve(L);
		v.push_back(scaled_copy(y0, image_dim(), k{1})); // v0 = y0
		for (std::size_t j = 1; j < L; ++j) {
			v.push_back(M_MT(v.back()));                 // vj = A * v_{j-1}
		}

		// y_-1 = scale * ( v[L-1] + c1*v[L-2] + ... + c_{L-1}*v[0] )
		ImageVec y_minus_1 = scaled_copy(v[L - 1], image_dim(), scale);

		for (std::size_t j = 1; j <= L - 1; ++j) {
			// connection_poly[j] is c_j, multiply v[L-1-j]
			add_scaled_inplace(y_minus_1, v[L - 1 - j], image_dim(), connection_poly[j] * scale);
		}

		// Lift to domain: x = M^T y_-1
		DomainVec result = M.evaluate_transpose_dense(y_minus_1);

		bool verified = verify_result(y0, result);

		if (!verified) {
			ImageVec y_found = M.evaluate_from_dense(result);

			std::cout << "False solution. connection_poly =" << std::endl;

			for (auto c : connection_poly) {
				std::cout << c << ", ";
			}
			std::cout << std::endl;

			std::cout << "y0 =    y_found =   " << std::endl;

			for (size_t i = 0; i < M.image_dim(); i++) {
				std::cout << y0[i] << "    " << y_found[i] << std::endl;
			}

			return std::nullopt;
		}

		return std::make_optional(std::move(result));
	}


      
    
    bool verify_result(const ImageVec& y0, const DomainVec& X) const {
			ImageVec y_found = M.evaluate_from_dense(X);
			
			k lambda = 0;
			
			size_t i = 0;
			
			for (; i< M.image_dim(); i++) {
					if (y0[i] != 0 && y_found[i] != 0) {
						lambda = y0[i]/y_found[i];
						break;
					} else if (y0[i] != 0 || y_found[i] != 0) {
						std::cout << "false solution type 1!" << std::endl;
						return false;
					}
			}
			
			for ( ;i< M.image_dim(); i++) {
					
					if (y0[i] != lambda* y_found[i]) {
						std::cout << "false solution type 2!" << std::endl;
					
						return false;
					}
			}
			
			std::cout << "Solution validated!" << std::endl;
			return true;
		
	} 


public:
    explicit wiedemann_solver(const lil_matrix<k>& matrix)
        : M(matrix) {
        // Allocate and fill random_vector
        random_vector = M.make_dense_image_vec_random();
    }
    

    std::optional<DomainVec> solve_MX_equals_y(const ImageVec& y0) {
        // Build scalar sequence s_i = <r, y_i>, y_{i+1} = (M M^T) y_i
        std::vector<k> signatures;
        signatures.reserve(64);

        signatures.emplace_back(get_signature(y0));

        ImageVec yi = M_MT(y0);
        signatures.emplace_back(get_signature(yi));

        berlekamp_massey_state<k> bm_state(signatures); // adjust template as needed
        bm_state.process_all_new();

        std::size_t slack = 8;
        std::size_t more_needed = bm_state.more_needed(slack);

        while (more_needed > 0) {
            for (std::size_t i = 0; i < more_needed; ++i) {
                yi = M_MT(yi);
                signatures.emplace_back(get_signature(yi));
            }
            bm_state.process_all_new();
            more_needed = bm_state.more_needed(slack);
        }
		
        return recompute_solution(y0, bm_state.connection_poly());
    }

    std::size_t domain_dim() const { return M.domain_dim(); }
    std::size_t image_dim()  const { return M.image_dim(); }
};

} // namespace VectorSpace
