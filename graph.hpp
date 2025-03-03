

#pragma once
#include <array>
#include<vector>
#include <iostream>
#include "CombinatorialUtils.cpp"
#include "Types.hpp"
#include <memory>

using namespace std;

template <
    Int N_VERTICES,
    Int N_EDGES,
    Int N_OUT_HAIR,
    Int N_IN_HAIR
>
class Graph
{
public:
    using SplitGraph = Graph<
        (N_VERTICES + 1),
        (N_EDGES + 1),
        N_OUT_HAIR,
        N_IN_HAIR 
    >;
    static constexpr Int SIZE = N_IN_HAIR + N_OUT_HAIR + 2 * N_EDGES;


    Graph() = default;
//private:
    // Two arrays, each with fixed size. The first N_EDGES entries will be our edges.
    // outArray[i] is the "from" vertex for edge i
    // inArray[i]  is the "to"   vertex for edge i
    array<Int, SIZE> half_edges{};

    explicit Graph(const array<Int, SIZE>& arr)
        : half_edges(arr)
    {}

    
    vector<Int> adjacent(Int v) const
    {
        vector<Int> result;
        result.reserve(SIZE);

        for (Int i = 0; i < SIZE; ++i) {
            if (half_edges[i] == v) {
                result.push_back(i);
            }
        }
        return result;
    }


    vector<unique_ptr<SplitGraph>> split_vertex_differential() const
    {
        vector<vector<Int>> adjRepresentation;
        adjRepresentation.reserve(N_VERTICES);

        bigInt resultSize = 0;

        for (Int v = 0; v < N_VERTICES; v++) {
            adjRepresentation.push_back(adjacent(v));
            resultSize += combutils::n_splits(adjRepresentation.back().size());
        }
        cout << "Should get " << resultSize << " splits"<< endl;

        vector<unique_ptr<SplitGraph>>  result;
        result.reserve(resultSize);

        for (Int v = 0; v < N_VERTICES; v++) {
            split_vertex(v, adjRepresentation[v], result);
        }

        return result;
    }


    void split_vertex(Int split_vertex, vector<Int> &adjacent, vector<unique_ptr<SplitGraph>> &result) const
    {
        if (adjacent.size() < 4) {
            cout << "BAJS!";
            return;
        }

        Int max_index = adjacent.size()-1;

        for (Int i = 2; i < max_index; i++) {
            vector<Int> S = combutils::firstSubset(1, i);

            do {
                result.push_back(splitGraph(split_vertex, adjacent[0], S));
                cout << endl;
                for (auto s : S) {
                    cout << s << ", " ;

                }
                cout << endl;
            } while (combutils::nextSubset(S, max_index));

        }
       
    }

    unique_ptr<SplitGraph> splitGraph(Int split_vertex, Int firstIndex, vector<Int> &S) const {
        auto sg = make_unique<SplitGraph>();

        // Copy half_edges from firstIndex onward.
        for (Int i = firstIndex; i < SIZE; ++i) {
            sg->half_edges[i] = half_edges[i];
        }

        // For each index in S, set that position in half_edges to N_VERTICES.
        for (auto s : S) {
            sg->half_edges[s] = N_VERTICES;
        }

        // Set the last two entries in half_edges:
        // - Second-to-last: split_vertex
        // - Last: N_VERTICES
        sg->half_edges[SplitGraph::SIZE - 2] = split_vertex;
        sg->half_edges[SplitGraph::SIZE - 1] = N_VERTICES;

        return sg;
    }


};