
#pragma once

#include "types.hpp"
#include <utility>
#include <vector>
#include <unordered_map>
#include <unordered_set>

#include "VectorSpace/LinComb.hpp"

using namespace::std;


template<typename A, typename B, typename k> 
class MetaGraph {
        using beA = BasisElement<A, k>;
        using lcB = VectorSpace::LinComb<B,k>;
        public:

        struct Edge {
            A x, y;
            unordered_set<B> common_boundary_components; // perhaps unnessesary 
        };

        vector<Edge> edges;
        
        vector<A> hair;

        auto begin()       { return edges.begin(); }
        auto end()         { return edges.end(); }

        MetaGraph(unordered_map<A, lcB> delta_of_remainder) {


            for (auto it1 = delta_of_remainder.begin(); it1 != delta_of_remainder.end(); ++it1) {
				
				if (it1 -> second.empty()) {
					hair.push_back(it1 -> first);
					continue;
				}
				
                auto it2 = it1;
                it2++;
                    for (; it2 != delta_of_remainder.end(); ++it2) {
                        unordered_set<B> common_boundary_components = intersection(it1 -> second, it2 -> second);
                        if (common_boundary_components.empty()) continue;

                        this-> edges.push_back(Edge{it1 -> first, it2 -> first, common_boundary_components});

                    }
            }
        }

        // two pointer style. Linear complexity. Can be changed to interval halving approach if we ever want large X,Y here 
        static unordered_set<B> intersection(lcB const& X, lcB const& Y) {
            unordered_set<B> result;

            auto itX  = X.begin();
            auto endX = X.end();
            auto itY  = Y.begin();
            auto endY = Y.end();

            while (itX != endX && itY != endY) {
                auto const& ex = *itX; // BasisElement<B,k>
                auto const& ey = *itY; // BasisElement<B,k>

                int cmp = ex.compare(ey);

                if (cmp < 0) {
                    // ex < ey
                    ++itX;
                } else if (cmp > 0) {
                    // ex > ey
                    ++itY;
                } else {
                    // ex == ey
                    result.insert(ex.getValue());   // take the underlying B
                    ++itX;
                    ++itY;
                }
            }

            return result;
    }


	unordered_set<A> component_containing(A a) {
		bool added_new = true;
		
		unordered_set<A> result;
		result.insert(a);
		
		while (added_new) {
			added_new = false;
			for (const Edge& edge: edges) {
				if (result.contains(edge.x)) {
						added_new |= result.insert(edge.y).second;
				} else if(result.contains(edge.y)) {
						added_new |= result.insert(edge.x).second;
				} 
			}
		}
		
		return result;
		
	}
};
