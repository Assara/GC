#pragma once

#include <vector>
#include <optional>
#include <random>
#include <cstddef>
#include <utility>
#include <iostream>
#include <chrono>

#include "lil_matrix.hpp"
#include "berlekamp_massey_state.hpp"


namespace VectorSpace {

template<typename k>
class block_wiedemann_solver {
public:
    using ImageVec  = typename lil_matrix<k>::DenseImageVec;   // expected: std::unique_ptr<k[]>
    using DomainVec = typename lil_matrix<k>::DenseDomainVec;  // expected: std::unique_ptr<k[]>
    
    
    struct timing_t {
    timer_accum eval_M;      // M * x
    timer_accum eval_MT;     // M^T * y
    timer_accum bm;          // Berlekamp–Massey 
};

	mutable timing_t timing;

	void print_timing() const {
		timing.eval_M.print("M * x");
		timing.eval_MT.print("M^T * y");
		timing.bm.print("berlekamp massey computations:");
	}

private:
    const lil_matrix<k>& M;
    const std::vector<ImageVec> random_vectors; // length = M.image_dim()
    
    
    const compressed_sparse_matrix<k> compressed_M;
    const compressed_sparse_matrix<k> compressed_transpose_M;
    
public:
    explicit block_wiedemann_solver(const lil_matrix<k>& matrix, size_t block_size)
        : M(matrix),
        random_vectors(M.make_dense_image_vectors_random(block_size)),
        compressed_M(M.to_compressed_sparse_matrix()),
        compressed_transpose_M(M.transpose().to_compressed_sparse_matrix()) {
    }

    

private:
	ImageVec M_MT(const ImageVec& y) const {
		
		DomainVec v;
		{
			timer_accum::guard g(timing.eval_MT);
			v = compressed_M.evaluate_transpose(y);
		}
		{
			timer_accum::guard g(timing.eval_M);
			return compressed_transpose_M.evaluate_transpose(v);
		}
	}
	
	DomainVec MT_M(const DomainVec& y) const {
		ImageVec v;
		{
			timer_accum::guard g(timing.eval_MT);
			v = compressed_transpose_M.evaluate_transpose(y);
		}
		{
			timer_accum::guard g(timing.eval_M);
			return compressed_M.evaluate_transpose(v);
		}
	}



    k get_signature(const ImageVec& y, size_t i) const {
        k result = k{0};
        const auto& random_vector = random_vectors[i];
        
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
    
	
		
	std::vector<DomainVec>
	recompute_solutions(const ImageVec& y0,
						const std::vector<std::vector<k>>& polys) const
	{
		std::cout << "recomputing solutions (domain-side): " << polys.size() << "\n";
		if (polys.empty()) return {};

		const std::size_t n = polys.size();

		// x0 := M^T y0  (computed once)
		DomainVec x0;
		{
			timer_accum::guard g(timing.eval_MT);
			x0 = compressed_M.evaluate_transpose(y0); // Image -> Domain
		}

		// degrees and scales
		std::vector<std::size_t> L(n, 0);
		std::vector<k> scale(n, k{0});
		std::size_t Lmax = 0;

		for (std::size_t a = 0; a < n; ++a) {
			const auto& C = polys[a];
			if (C.size() < 2) {
				std::cout << "poly[" << a << "] too small (size=" << C.size() << ")\n";
				continue;
			}

			const std::size_t La = C.size() - 1;
			const k cL = C.back();
			if (cL == k{0}) {
				std::cout << "poly[" << a << "] singular (cL=0)\n";
				continue;
			}

			L[a] = La;
			Lmax = std::max(Lmax, La);

			// same scaling convention as your original routine
			scale[a] = -(cL.inv()); // -1/cL
		}

		if (Lmax == 0) return {};

		// One accumulator per poly in domain space
		std::vector<DomainVec> acc;
		acc.reserve(n);
		for (std::size_t a = 0; a < n; ++a) {
			acc.push_back(scaled_copy(x0, domain_dim(), k{0})); // zero
		}

		// Stream z = (M^T M)^t x0
		DomainVec z = scaled_copy(x0, domain_dim(), k{1}); // t = 0

		for (std::size_t t = 0; t < Lmax; ++t) {
			for (std::size_t a = 0; a < n; ++a) {
				const std::size_t La = L[a];
				if (La == 0 || t >= La) continue;

				const auto& C = polys[a];

				// Coefficient mapping consistent with your Horner semantics:
				// acc = B^{La-1}x0 + C[1]B^{La-2}x0 + ... + C[La-1]x0
				const k coeff = (t == La - 1) ? k{1} : C[La - 1 - t];

				if (coeff != k{0}) {
					add_scaled_inplace(acc[a], z, domain_dim(), coeff);
				}
			}

			if (t + 1 < Lmax) {
				z = MT_M(z); // z <- (M^T M) z
			}
		}

		// Final candidates: x_a = scale[a] * acc[a]
		std::vector<DomainVec> out;
		out.reserve(n);

		for (std::size_t a = 0; a < n; ++a) {
			if (L[a] == 0) continue;

			DomainVec x = scaled_copy(acc[a], domain_dim(), scale[a]);

			if (!verify_solution(y0, x)) {
				std::cout << "False solution for poly[" << a << "]\n";
				continue;
			}

			out.push_back(std::move(x));
		}

		return out;
	}


    bool verify_solution(const ImageVec& y0, const DomainVec& X) const {
			ImageVec y_found = M.evaluate_from_dense(X);
		
			size_t i = 0;
			
			
			
			for ( ;i< M.image_dim(); i++) {
					
					if (y0[i] != y_found[i]) {
						std::cout << "false solution type 2!" << std::endl;
					
						return false;
					}
			}
			
			std::cout << "Solution validated!" << std::endl;
			
			
			return true;	
	} 


	size_t block_size() const {
			return random_vectors.size();
	}
	
	
	void add_signatures(ImageVec& y_i, std::vector<std::vector<k>>& signature_vectors) const {
			for (size_t i = 0; i<block_size(); i++) {
				signature_vectors[i].emplace_back(get_signature(y_i, i));
			}
	}
			
			
		

public:
    std::vector<DomainVec> solve_MX_equals_y(const ImageVec& y0) {
        // Build scalar sequence s_i = <r, y_i>, y_{i+1} = (M M^T) y_i
        if (image_dim() == 0) {
				std::vector<DomainVec> out;
				out.emplace_back(M.make_dense_domain_vec_zero());
				
				return out;
		}
		
	
		
		std::cout << "using solve_MX_equals_y. ";
		std::cout << "image dim = " << image_dim();
		std::cout << " domain dim = " << domain_dim() << std::endl;
		
		bool is_zero = true;
		
		for (std::size_t i = 0; i< image_dim(); i++) {
				if (y0[i] != 0) {
					is_zero = false;
					break;
				} 
		}
		
		if (is_zero) {
				std::vector<DomainVec> out;
				out.emplace_back(M.make_dense_domain_vec_zero());
				return out;
		};
		
		std::cout << std::endl;
        std::vector<std::vector<k>> signatures_vector(block_size());
        
        
        std::size_t slack = 8;
        ImageVec yi = M_MT(y0);
        
        std::vector<berlekamp_massey_state<k>> bm_states;
        
        for (size_t i = 0; i<block_size(); i++) {	
			signatures_vector[i].reserve(2* M.image_dim() + slack);
			signatures_vector[i].emplace_back(get_signature(y0, i));
			
			signatures_vector[i].emplace_back(get_signature(yi, i));
			
			
			bm_states.emplace_back(berlekamp_massey_state<k>(signatures_vector[i]));
			bm_states.back().process_all_new();
		}
		
  
        std::size_t more_needed = bm_states[0].more_needed(slack);

        while (more_needed > 0) {
            for (std::size_t i = 0; i < more_needed; ++i) {
                yi = M_MT(yi);
                add_signatures(yi, signatures_vector);    // add mask for completed signature vectors?
            }
			{
				timer_accum::guard g(timing.bm);
				more_needed = 0;
				for (auto& bm_state : bm_states) {
					bm_state.process_all_new();
					more_needed = std::max(more_needed, bm_state.more_needed(slack));
				}
			}
        }
        
        std::vector<std::vector<k>> connection_polynomials;
        connection_polynomials.reserve(block_size());
        for (auto& bm_state: bm_states) {
				connection_polynomials.emplace_back(bm_state.connection_poly());
		} 
		
        return recompute_solutions(y0, connection_polynomials);
    }

    std::size_t domain_dim() const { return M.domain_dim(); }
    std::size_t image_dim()  const { return M.image_dim(); }
};

} // namespace VectorSpace
