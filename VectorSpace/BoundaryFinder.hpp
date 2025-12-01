#pragma once
#include "LinComb.hpp"
#include "FiniteMatrix.hpp"   // VectorSpace::FiniteMatrix<k> (single-block, column-major)
#include <algorithm>
#include <vector>
#include <unordered_map>

#include "../GC.hpp"

using namespace std;

namespace VectorSpace {

template<typename A, typename B, typename k>
class BoundaryFinder {

    using FM  = FiniteMatrix<k>;
    using Row = typename FM::Row;  // length = deltaRep.cols()
    using Col = typename FM::Col;   // length = deltaRep.rows()

    private:
    vector<B> coBoundaryRepMap;

    FM deltaRep;
    unordered_map<A, size_t> boundaryRepMap;

    public:

    BoundaryFinder(unordered_map<B, VectorSpace::LinComb<A,k>>& delta) {
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


            for (size_t i = 0; i < deltaRep.rows(); ++i) {
                    const LinComb<A,k>& linComb = delta[coBoundaryRepMap[i]];

                    for (const BasisElement<A,k>& elem : linComb) {
                            deltaRep.set(i, boundaryRepMap.at(elem.getValue()), elem.getCoefficient());
                    }
            }
    }


        BoundaryFinder(unordered_map<B, VectorSpace::LinComb<A,k>>& delta, unordered_set<A> filter) {
                coBoundaryRepMap.reserve(delta.size());
                size_t a = 0;

                cout << "filter size =" << filter.size() << endl;

                cout << "delta size =" << delta.size()<< endl;


                for (auto& graph : filter) {
                        if (!boundaryRepMap.contains(graph)) {
                                cout << "the mapped part of Delta is insufficient for the filter! " << endl;
                        }
                }

                //construct the in out maps
                for (const auto& entry : delta) {
                        coBoundaryRepMap.push_back(entry.first);             
        
                        for (const BasisElement<A,k>& elem : entry.second) {
                                if (!filter.contains(elem.getValue())) continue;
                                
                                if (auto [it, inserted] = boundaryRepMap.try_emplace(elem.getValue(), a); inserted) {
                                        ++a;
                                } 
                        }

            }

            // create the internal Map 
            deltaRep = FM(coBoundaryRepMap.size(), boundaryRepMap.size());


            for (size_t i = 0; i < deltaRep.rows(); ++i) {
                    const LinComb<A,k>& linComb = delta[coBoundaryRepMap[i]];

                    for (const BasisElement<A,k>& elem : linComb) {
                            const auto& key = elem.getValue();      
                            if (auto it = boundaryRepMap.find(key); it != boundaryRepMap.end()) {
                                    deltaRep.set(i, it->second, elem.getCoefficient());
                            }
                    }
            }

            cout << "Created filtered boundary finder with dimensions. rows: "<< deltaRep.rows() << " cols: " << deltaRep.cols() <<endl;
    }

    std::optional<LinComb<B,k>> find_primitive_or_empty(const LinComb<A,k>& linComb) {
            Row rep = std::make_unique<k[]>(deltaRep.cols());

            for (const BasisElement<A,k>& elem : linComb) {
                     const auto& key = elem.getValue(); 
                     if (auto it = boundaryRepMap.find(key); it != boundaryRepMap.end()) {
                            rep[it->second] = elem.getCoefficient();
                     }
            }

            optional<Col> coBoundaryRep = deltaRep.solve_exact_MT(rep);

            if (!coBoundaryRep.has_value()) {
                cout << "COULD NOT FIND PRIMITIVE!" << endl;
                 return std::nullopt;
            }

            Col coeffs = std::move(coBoundaryRep.value()); 
     
            LinComb<B, k> primitive;

            for (size_t i = 0; i < deltaRep.rows(); ++i) {
                //not without standardizing or deduplicating
                if (coeffs[i] == 0) continue;
                primitive.append_in_basis_order(coBoundaryRepMap[i], -coeffs[i]);
            }
            //not necessarly sorted, but that is fine
            return primitive;
    }
};

};