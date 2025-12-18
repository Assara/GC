#pragma once

#include <vector>
#include <cstddef>

template<typename k>
class berlekamp_massey_state {
public:
    using scalar_type = k;

    // S_ref must outlive this object
    explicit berlekamp_massey_state(std::vector<k>& S_ref)
        : S(S_ref),
          C{ k{1} },  // C(x) = 1
          B{ k{1} },  // B(x) = 1
          L(0),
          m(1),
          b(k{1}),
          n_processed(0)
    {}

    // Process all currently available elements in S
    void process_all_new() {
        process_up_to(S.size());
    }

    // Process S[n_processed .. N-1]; N must be <= S.size()
    void process_up_to(std::size_t N) {
        if (N > S.size()) N = S.size();
        if (n_processed >= N) return;

        for (std::size_t n = n_processed; n < N; ++n) {
            // discrepancy d at position n:
            // d = S[n] + sum_{i=1}^L C[i] * S[n-i]
            k d = S[n];
            for (std::size_t i = 1; i <= L; ++i) {
                d = d + C[i] * S[n - i];
            }

            if (d == k{}) {
                // Recurrence is consistent at n
                ++m;
                continue;
            }

            // Save current C
            std::vector<k> T = C;

            // Factor for update
            k factor = d / b;

            // Ensure C is large enough to accommodate B shifted by m
            if (C.size() < B.size() + m) {
                C.resize(B.size() + m, k{});
            }

            // C[j + m] -= factor * B[j]
            for (std::size_t j = 0; j < B.size(); ++j) {
                C[j + m] = C[j + m] - factor * B[j];
            }

            // Decide whether to increase recurrence length
            if (2 * L <= n) {
                // Increase recurrence length
                std::size_t L_new = n + 1 - L;
                B = std::move(T);
                b = d;
                L = L_new;
                m = 1;
            } else {
                ++m;
            }
        }

        n_processed = N;
    }
    
	std::size_t more_needed(std::size_t slack = 0) const noexcept {
		// No recurrence yet: force at least one more sample
		if (L == 0) {
			return 1;
		}

		std::size_t target = 2 * L + slack;

		if (n_processed >= target) {
			return 0;
		}

		return target - n_processed;
	}


    // Return current connection polynomial trimmed to degree L:
    // C[0..L], C[0] = 1
    std::vector<k> connection_poly() const {
        std::vector<k> res = C;
        if (res.size() > L + 1) {
            res.resize(L + 1);
        }
        return res;
    }

    // Order of the recurrence (degree)
    std::size_t order() const noexcept {
        return L;
    }

    // Number of sequence elements incorporated so far
    std::size_t processed_count() const noexcept {
        return n_processed;
    }

private:
    std::vector<k>& S;   // reference to the sequence

    std::vector<k> C;    // current connection polynomial
    std::vector<k> B;    // backup polynomial from last "big" update
    std::size_t   L;     // current recurrence length
    std::size_t   m;     // steps since last update of B
    k             b;     // last nonzero discrepancy (for scaling)
    std::size_t   n_processed; // how many S[n] have been processed
};
