
#include <array>
#include <vector>
#include <iostream>
#include <numeric>
#include <cstring>
#include "types.hpp"

using namespace std;

namespace combutils {
    inline bool nextSubset(vector<Int>& S, Int A_max_index) {
        Int i = 0;
        for (; i < S.size(); ++i) {
            if (S[i] < A_max_index - i) {
                S[i]++;
                break;
            }
        }
        if (i == S.size()) return false;
        while(i > 0) {
            --i;
            S[i] = S[i+1] + 1;
    
        }
        return true;
    }

    inline vector<Int> firstSubset(Int startIndex, Int size) {
        vector<Int> subset(size);
        for (Int i = 0; i < size; ++i) {
            subset[i] = startIndex + (size - 1 - i);
        }
        return subset;
    }


    inline bigInt n_splits(bigInt n_adjacent) {
        if (n_adjacent < 4) {
            return 0;
        }

        return (bigInt(1) << (n_adjacent - 1)) - n_adjacent - 1;
    }


    template<size_t N>
    inline signedInt compareHalfEdges(const std::array<Int, N>& a, const std::array<Int, N>& b) {
        return std::memcmp(a.data(), b.data(), N * sizeof(Int));

    }


    template <std::size_t N>
    [[nodiscard]] inline constexpr bool
    lessThan(const std::array<Int, N>& a, const std::array<Int, N>& b) noexcept {
        for (std::size_t i = 0; i < N; ++i) {
            if (a[i] < b[i]) return true;   // first differing element decides
            if (a[i] > b[i]) return false;
        }
        return false; // equal
    }

    template<typename T>
    vector<T> intersection(unordered_set<T> const& a,
                                unordered_set<T> const& b){
        vector<T> result;

        // iterate over the smaller set for efficiency
        auto const& small = (a.size() < b.size()) ? a : b;
        auto const& large = (a.size() < b.size()) ? b : a;

        for (auto const& x : small)
            if (large.contains(x))
                result.push_back(x);

        return result;
    }
    
}
