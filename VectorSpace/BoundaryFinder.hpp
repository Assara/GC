#pragma once
#include "LinComb.hpp"
#include "FiniteMatrix.hpp"   // VectorSpace::FiniteMatrix<k> (single-block, column-major)
#include <algorithm>
#include <vector>
#include <unordered_map>

using namespace std;

namespace VectorSpace {

template<typename A, typename B, typename k>
class BoundaryFinder {

    using FM  = FiniteMatrix<k>;
    using Row = typename FM::Row;  // length = cols()
    using Col = typename FM::Col; 

    private:
    vector<B> coBoundaryRepMap;

    FM deltaRep;
    unordered_map<A, size_t> boundaryRepMap;

    public:

    BoundaryFinder(unordered_map<B, VectorSpace::LinComb<A,k>> delta) {
        coBoundaryRepMap.reserve(delta.size());
        size_t a = 0;

        //construct the in out maps
        for (const auto& entry : delta) {
                coBoundaryRepMap.push_back(entry.first);
                  
                for (const BasisElement<A,k>& elem : entry.second) {
                        if (auto [it, inserted] = boundaryRepMap.try_emplace(elem.getValue(), a); inserted) {
                            ++a;
                        }
                }
        }

        // create the internal Map 
        deltaRep = FM(coBoundaryRepMap.size(), boundaryRepMap.size());
        for (size_t i = 0; i < coBoundaryRepMap.size(); ++i) {
            const LinComb<A,k>& linComb = delta[coBoundaryRepMap[i]];
            for (const BasisElement<A,k>& elem : linComb) {
                    deltaRep.set(i, boundaryRepMap.at(elem.getValue()), elem.getCoefficient());
            }
        }
    }

    std::optional<LinComb<B,k>> find_coboundary_or_empty(LinComb<A,k> linComb) {
            Row rep = std::make_unique<k[]>(boundaryRepMap.size());

            for (const BasisElement<A,k>& elem : linComb) {
                    rep[boundaryRepMap[elem.getValue()]] = elem.getCoefficient();
            }

            optional<Row> coBoundaryRep = deltaRep.solve_exact_then_mul_T(rep);

            if (!coBoundaryRep.has_value()) {
                 return std::nullopt;
            }

            k* coeffs = coBoundaryRep->get(); 

            cout  <<"coeffs = " ;
            for (size_t i = 0 ; i < coBoundaryRepMap.size(); i++) {
                cout <<  coeffs[i] << ", ";

            }
            cout << endl;

            LinComb<B, k> coBoundary;

            for (size_t i = 0; i < coBoundaryRepMap.size(); ++i) {
                //not without standardizing or deduplicating
                coBoundary.append_in_basis_order(coBoundaryRepMap[i], -coeffs[i]);
            }

            //sort coBoundary but no need to deduplicate
            coBoundary.sort_without_deduplicate();

            return std::move(coBoundary);
    }



};

};