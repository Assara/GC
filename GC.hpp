#pragma once

#include "GraphStandardizer.hpp"
#include "VectorSpace/LinComb.hpp"
#include "VectorSpace/FiniteSubSpace.hpp"
#include "VectorSpace/BasisElement.hpp"

#include "VectorSpace/BoundaryFinder.hpp"

#include "MetaGraph.hpp"

template <
    Int N_VERTICES, Int N_EDGES,
    Int N_OUT_HAIR, Int N_IN_HAIR,
    signedInt c, signedInt d
>
class GC {

public:
        using ThisGC = GC<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d>;
        using GraphType = Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d>;
        using SplitGC = GC<N_VERTICES + 1, N_EDGES + 1, N_OUT_HAIR, N_IN_HAIR, c, d>;
        using ContGC = GC<N_VERTICES - 1, N_EDGES - 1, N_OUT_HAIR, N_IN_HAIR, c, d>;

        using ExtraEdgeGC = GC<N_VERTICES, N_EDGES + 1, N_OUT_HAIR, N_IN_HAIR, c, d>;

        using SplitGraphType = typename GraphType::SplitGraph;
        using ContGraphType = typename GraphType::ContGraph;
        using ExtraEdgeGraphType = typename GraphType::ExtraEdgeGraph;

        using L      = VectorSpace::LinComb<GraphType, fieldType>;
        using ContL  = VectorSpace::LinComb<ContGraphType, fieldType>;
        using SplitL = VectorSpace::LinComb<SplitGraphType, fieldType>;


private:
        VectorSpace::LinComb<GraphType, fieldType> vec;

public:
        // Default constructor
        GC() = default;

        // Construct from a single BasisElement
        explicit GC(const BasisElement<GraphType, fieldType>& b) : vec(b) {}

        explicit GC(const BasisElement<GraphType, fieldType>& b, AssumeBasisOrderTag) : vec(b, AssumeBasisOrderTag{}) {}

        // Construct from a GraphType (copying it into a BasisElement)
        explicit GC(const GraphType& G) : vec(G) {}

        explicit GC(const GraphType& G, AssumeBasisOrderTag) : vec(G, AssumeBasisOrderTag{}) {}


        // Construct from a ptr to GraphType (moving it into a BasisElement)
        explicit GC(std::unique_ptr<GraphType>&& g, fieldType coeff = 1.0f)
        : vec(BasisElement<GraphType, fieldType>(std::move(g), coeff)) {}

        // construct from a LinComb
        explicit GC(const VectorSpace::LinComb<GraphType, fieldType>& v) : vec(v)  {}

        // Construct from a vector of unique_ptrs to GraphType and a coefficient
        explicit GC(std::vector<std::unique_ptr<GraphType>>&& graphs, fieldType coeff)
                        : vec(std::move(graphs), coeff) {} 

        explicit GC(std::vector<BasisElement<GraphType, fieldType>>&& elems)
                        : vec(std::move(elems)) {} 

        explicit GC(std::vector<BasisElement<GraphType, fieldType>>&& elems,
                        AssumeBasisOrderTag)                    
                : vec(std::move(elems), AssumeBasisOrderTag{})  
        {}

        // Access the underlying vector
        const VectorSpace::LinComb<GraphType, fieldType>& data() const {
                return vec;
        }

        void standatdize_all() {
                vec.standardize_all();
        }

        // += operator for summing GCs
        ThisGC& operator+=(const ThisGC& other) {
                vec += other.vec;
                return *this;
        }

        // Standardize and sort elements
        void standardize_all() { vec.standardize_all(); }
        void sort_elements()   { vec.sort_elements(); }
        void standardize_and_sort() { vec.standardize_and_sort(); }

        BasisElement<GraphType, fieldType>& back() {
                return vec.back();
        }

        // Compute delta for the whole GC
        SplitGC delta() const {
                return delta_recursive(0, vec.raw_elements().size());
        }

        // Static method to compute delta of a single basis element
        static SplitGC delta(const BasisElement<GraphType, fieldType>& G) {
                return SplitGC(G.getValue().split_vertex_differential(), G.getCoefficient());
        }

        ExtraEdgeGC add_edge_differential() {
                return add_edge_differential_recursive(0, vec.raw_elements().size());
        }

        static ExtraEdgeGC add_edge_differential(const BasisElement<GraphType, fieldType>& G) {
                return ExtraEdgeGC(G.getValue().add_edge_differential(), G.getCoefficient());
        }


        vector<BasisElement<typename GraphType::ContGraph, fieldType>> d_contraction_without_sort() const {
                using ContGraph = typename GraphType::ContGraph;
                vector<BasisElement<typename GraphType::ContGraph, fieldType>> result;
                result.reserve(vec.raw_elements().size() * GraphType::N_EDGES_);

                for (const auto& elem : vec.raw_elements()) {
                        const auto& be = elem;
                        for (Int i = 0; i < GraphType::N_EDGES_; ++i) {
                                // Contract the i-th edge
                                BasisElement<ContGraph, fieldType> contracted = GraphType::contract_edge(be, i);
                                
                                if (contracted.getCoefficient() != 0) {
                                        result.push_back(std::move(contracted));;
                                }
                        }
                }
                return result;
        }
        
        vector<BasisElement<typename GraphType::ContGraph, fieldType>> d_even_contraction_without_sort() const {
                using ContGraph = typename GraphType::ContGraph;
                vector<BasisElement<typename GraphType::ContGraph, fieldType>> result;
                result.reserve(vec.raw_elements().size() * GraphType::N_EDGES_);

                for (const auto& elem : vec.raw_elements()) {
                        const BasisElement<GraphType, fieldType>& be = elem;

                        auto valence_array = be.getValue().valence_array();
                        for (Int i = 0; i < GraphType::N_EDGES_; ++i) {
                                auto edge = be.getValue().getEdge(i);
                                
                                if(valence_array[edge.first] %2 == 1 && valence_array[edge.second] %2 == 1) {
                                        continue;
                                }

                                // Contract the i-th edge
                                BasisElement<ContGraph, fieldType> contracted = GraphType::contract_edge(be, i);
                                
                                if (contracted.getCoefficient() != 0) {
                                        result.push_back(std::move(contracted));;
                                }
                        }
                }
                return result;
        }


        ContGC d_contraction() {
                std::vector<BasisElement<typename GraphType::ContGraph, fieldType>> elems = d_contraction_without_sort();
                ContGC dThis(std::move(elems));
                return dThis;
        }

        ContGC add_contractions_to_set(unordered_set<ContGraphType>& seenGraphs) {
                std::vector<BasisElement<ContGraphType, fieldType>> elems = d_contraction_without_sort();
                ContGC dThis(std::move(elems), AssumeBasisOrderTag{});
                
                dThis.standardize_all();

                for (auto b : dThis.data()) {
                        seenGraphs.insert(b.getValue());

                }
                return dThis;
        }

        ContGC d_contraction_with_recording_seen_gaphs(unordered_set<ContGraphType>& seenGraphs) {
                ContGC dThis = add_contractions_to_set(seenGraphs);
                dThis.sort_elements();
                return dThis;
        }

 


        ContGC d_even_contraction_with_recording_seen_gaphs(unordered_set<ContGraphType>& seenGraphs) {
                std::vector<BasisElement<ContGraphType, fieldType>> elems = d_even_contraction_without_sort();
                ContGC dThis(std::move(elems), AssumeBasisOrderTag{});
                
                dThis.standardize_all();

                for (auto b : dThis.data()) {
                        seenGraphs.insert(b.getValue());

                }



                dThis.sort_elements();
                return dThis;
        }

        
        ContGC d_even_contraction() {
                std::vector<BasisElement<typename GraphType::ContGraph, fieldType>> elems = d_even_contraction_without_sort();
                ContGC dThis(std::move(elems));
                return dThis;
        }

        void add(BasisElement<GraphType, fieldType>&& elem) {
                vec.add(VectorSpace::LinComb<GraphType, fieldType>(std::move(elem)));
        }


        std::optional<ContGC> try_find_split_primitive() {
                unordered_set<ContGraphType> seen_graphs;
                add_contractions_to_set(seen_graphs);
                unordered_map<ContGraphType, L> coboundary_map;

                for (auto& gamma : seen_graphs) {
                        coboundary_map.emplace(gamma, ContGC(gamma, AssumeBasisOrderTag{}).delta().data());
                }

                VectorSpace::BoundaryFinder solver(coboundary_map);
                       
                std::optional<ContL> primitive_optional = solver.find_primitive_or_empty(this -> data());

                return primitive_optional.transform([](ContL lin_comb) { 
                        return ContGC(lin_comb);
                });
        }

        bigInt size() {
                return vec.size();
        }

        Int frontValence() {
                return vec.front().getValue().valence_array()[0];
        }

        void expand_map(std::unordered_map<SplitGraphType, L>& boundary_map,
                        std::unordered_set<GraphType>&            cousins,
                        std::unordered_set<GraphType>&            has_split){
                std::unordered_set<SplitGraphType> to_contract;
                
                for (const GraphType& graph : cousins) {
                        if (has_split.contains(graph)) continue;

                        graph.add_even_splits_to_set(to_contract);
                        has_split.insert(graph);

                        for (const SplitGraphType& sg : to_contract) {
                                if (boundary_map.contains(sg)) continue;

                                // keep a local “seen” set in the contracted world
                                std::unordered_set<GraphType> seen_contractions;

                                // use the correct tag type for SplitGC
                                auto contracted = SplitGC(sg, AssumeBasisOrderTag{})
                                                        .d_even_contraction_with_recording_seen_gaphs(seen_contractions)
                                                        .data();


                                // construct the value (ThisGC) in-place
                                boundary_map.try_emplace(sg, std::move(contracted));
                        }
                        to_contract.clear();
                }
        }

        static void expand_map2(unordered_map<SplitGraphType, L>& boundary_map, L remainder) {

                unordered_map<GraphType, ContL> delta_remainder;
                delta_remainder.reserve(remainder.size());

                unordered_map<GraphType, unordered_set<SplitGraphType>> split_cache;
                split_cache.reserve(remainder.size());
                
                for (const auto& be: remainder) {
                        delta_remainder.emplace(be.getValue(), 
                                ThisGC(be.getValue(), AssumeBasisOrderTag{}).d_even_contraction().data());

                        be.getValue().add_even_splits_to_set(split_cache[be.getValue()]);
                }

                //add single graph cycle splits 
                for (auto entry : delta_remainder) {
                        if(entry.second.size() == 0) {
                                add_splits_to_boundary_map(boundary_map, BasisElement<GraphType, fieldType>(entry.first));
                        }
                }
     
                MetaGraph<GraphType, ContGraphType, fieldType> metaGraph(delta_remainder);

                for (const auto& metaEdge : metaGraph) {
                        auto common_splits = combutils::intersection(split_cache[metaEdge.x], split_cache[metaEdge.y]);

                        for (auto split : common_splits) {
                                add_graph_to_boundary_map(boundary_map, split);
                        }
                        if (!common_splits.empty()) continue;

                        //create upper delta
                        
                        unordered_map<SplitGraphType, L> upper_delta;
                        upper_delta.clear();
                        for (auto split : split_cache[metaEdge.x]) {
                                upper_delta.emplace(split, 
                                        SplitGC(split, AssumeBasisOrderTag{}).d_even_contraction().data());

                        }
                        for (auto split : split_cache[metaEdge.y]) {
                                upper_delta.emplace(split, 
                                        SplitGC(split, AssumeBasisOrderTag{}).d_even_contraction().data());

                        }

                        MetaGraph<SplitGraphType, GraphType, fieldType> upperMetaGraph(upper_delta);

                        for (auto edge : upperMetaGraph) {
                                add_graph_to_boundary_map(boundary_map, edge.x);
                                add_graph_to_boundary_map(boundary_map, edge.y);
                        }

                }
        }


        static bool add_splits_if_cycle(unordered_map<SplitGraphType, L>& boundary_map, const BasisElement<GraphType, fieldType>& be) {
                ContGC d_even_be = ThisGC(be, AssumeBasisOrderTag{}).d_even_contraction();
                if (d_even_be.size() != 0) {
                        return false;
                }

                add_splits_to_boundary_map(boundary_map, be);
                return true;
        }

        static void add_splits_to_boundary_map(unordered_map<SplitGraphType, L>& boundary_map, const BasisElement<GraphType, fieldType>& be) {
                std::unordered_set<SplitGraphType> to_contract;
                be.getValue().add_even_splits_to_set(to_contract);

                for (auto G : to_contract) {
                        add_graph_to_boundary_map(boundary_map, G);
                }
        }

        static void add_graph_to_boundary_map(unordered_map<SplitGraphType, L>& boundary_map, const SplitGraphType& G) {
                if (boundary_map.contains(G)) {
                        return;
                }
                
                boundary_map.emplace(G, SplitGC(G, AssumeBasisOrderTag{}).d_even_contraction().data());
        }


        ThisGC reduce2() {
                unordered_map<SplitGraphType, L> boundary_map;
                unordered_set<GraphType> seen_graphs;

        
                signedInt grade = vec.front().getValue().custom_filter();

   
                size_t max_depth = 10;

                ThisGC boundary;
                for (size_t i = 0; i < max_depth ; i++) {
                        cout << "depth = " <<  i << endl;
          

                        auto it = std::find_if(
                                vec.begin(), vec.end(),
                                [grade](auto const& elem) {
                                        return elem.getValue().custom_filter() < grade;
                                });


                        L top_grade_comb(vec.begin(), it);

                                
                        cout << "top_grade_comb.front() grade : " << top_grade_comb.front().getValue().custom_filter() << endl;

                        cout << "top_grade_comb.back() grade : " << top_grade_comb.back().getValue().custom_filter() << endl;


                        ContGC shouldBe0 = ThisGC(top_grade_comb).d_even_contraction();

                        

                        if(shouldBe0.size() > 0) {

                                cout << endl << endl;
                                cout << "WARNING the top grad combination does not satisfy requierments, something is wrong!" << endl;
                                shouldBe0.print();

                                cout << "see graph above" << endl;
                                


                        } else {
                                cout << "Top grade comb is good!" << endl;
                        }

                        cout << "top_grade_comb.size() = " << top_grade_comb.size() << endl;
                        top_grade_comb.print();


                        expand_map2(boundary_map, top_grade_comb);

                        bool exists = true;

                        for (const auto& be : top_grade_comb) {
                                exists = false;
                                for (const auto& entry : boundary_map) {
                                        exists = std::binary_search(entry.second.begin(), entry.second.end(), be);
                                        if (exists) break;
                                }
                                if (!exists) {
                                         cout << "insufficient map!!!!!! " << endl;
                                         cout << "missing graph: " << endl;
                                         be.getValue().print();
                                }
                        }

                        
                        cout << "boundary_map.size() = " << boundary_map.size() <<  endl;
                        VectorSpace::BoundaryFinder solver(boundary_map);
                  
                        auto primitive = solver.find_primitive_or_empty(top_grade_comb);
                       
                        if (primitive.has_value()) {
                                primitive -> print();
                
                                SplitGC primitive_as_GC = SplitGC(std::move(*primitive));
                            
                                
                                primitive_as_GC.print();
                    
                                boundary += primitive_as_GC.d_contraction();
                                boundary.print();
                        
                                return *this += boundary;
                        } 


                        while (true) {
                                //Find a boundary that covers the seen graphs
                                for (auto& be : top_grade_comb) {
                                        seen_graphs.insert(be.getValue());
                                }
                                VectorSpace::BoundaryFinder filtered_solver(boundary_map, seen_graphs);
                                primitive = filtered_solver.find_primitive_or_empty(top_grade_comb);
                                
                                
                                if (!primitive.has_value()) {
                                      break;
                                }
                                
                                SplitGC primitive_as_GC = SplitGC(std::move(*primitive));
                                boundary = primitive_as_GC.d_contraction();

                                if(boundary.d_contraction().size() > 0) {

                                        cout << endl << endl;
                                        cout << "WARNING boundary is not a cycle!" << endl;
                                        boundary.d_contraction().print();

                                        cout << "see graph above" << endl;
                                        

                                } else {
                                        cout << "Boundary comb is good!" << endl;
                                }


                                *this += boundary; 
                                for (const auto& be : this -> data()) {
                                        if (seen_graphs.contains(be.getValue())) {
                                                be.getValue().print();
                                                cout << "coefficient = " << be.getCoefficient() << endl;
                                                cout << "MATH ERROR" << endl;
                                        } 
                                }

                                auto it = std::find_if(
                                        vec.begin(), vec.end(),
                                        [grade](auto const& elem) {
                                                return elem.getValue().custom_filter() < grade;
                                        });

                                top_grade_comb = L(vec.begin(), it);
                        }
                         
                }

                return *this;
        
                
        }



        ThisGC reduce() {
                if (vec.raw_elements().empty()) {
                        return ThisGC{};
                }
         
                ThisGC yhat = ThisGC{};
                unordered_set<SplitGraphType> split_set;
                unordered_set<SplitGraphType> has_added_boundary;
                unordered_set<GraphType> to_split;
                unordered_set<GraphType> has_split;
                
                vector<VectorSpace::LinComb<GraphType, fieldType>> boundaries;


                signedInt grade = vec.front().getValue().custom_filter();


                for (auto b : vec) {
                        if (b.getValue().custom_filter() != grade) {
                                break;
                        }
                
                        to_split.insert(b.getValue());
                }
                
                const Int maxDepth = 2;
                int depth = 0;

                while (++depth <= maxDepth) {
                        /*for (VectorSpace::LinComb<GraphType, fieldType> boundary : boundaries) {
                                for (BasisElement<GraphType, fieldType> be : boundary) {
                                        if (vec.front().getValue() < be.getValue() && !has_split.contains(be.getValue())) {
                                                to_split.insert(be.getValue());
                                        }
                                }
                        }*/
                        cout << "to_split.size() = " << to_split.size() << endl;
                        cout << "depth = " << depth << endl << endl;
                        for (GraphType graph : to_split) {
                                if (!has_split.contains(graph) && graph.custom_filter() == grade) {
                                        graph.add_even_splits_to_set(split_set);
                                        has_split.insert(graph);
                                }
                        }

                        to_split.clear();
                


                        boundaries.reserve(split_set.size());

                        for (auto& split : split_set) {
                                 if(!has_added_boundary.contains(split)) {
                                        has_added_boundary.insert(split);
                                        SplitGC singleSplit = ThisGC::SplitGC(std::move(split));
                                        boundaries.emplace_back(singleSplit.d_contraction_with_recording_seen_gaphs(to_split).data());
                                        
                                }
                                // split_set.erase(split);
                        }
                        split_set.clear();
    
                        VectorSpace::FiniteSubSpace<GraphType, fieldType> subSpace(boundaries);
                        yhat = ThisGC(subSpace.project_solve_then_map(vec));
                        cout << endl << "yhat.size() = " << yhat.vec.size() <<endl;
                        
                        *this += yhat;

                        if (vec.front().getValue().custom_filter() < grade) {
                                break;
                        }
                }
                return *this;
              
        }

        GC& scalar_multiply(fieldType scalar) {
                vec.scalar_multiply(scalar);
                return *this;
        }


private:

        // Recursive delta helper
        SplitGC delta_recursive(size_t start, size_t end) const {
                if (end - start == 0) {
                        return SplitGC(); // empty
                }
                if (end - start == 1) {
                        const auto& be = vec.raw_elements()[start];
                        return delta(be);
                }

                size_t mid = start + (end - start) / 2;
                SplitGC left = delta_recursive(start, mid);
                SplitGC right = delta_recursive(mid, end);
                left += right;
                return left;
        }


        ExtraEdgeGC add_edge_differential_recursive(size_t start, size_t end) const {
                if (end - start == 0) {
                        return ExtraEdgeGC(); // empty
                }
                if (end - start == 1) {
                        const auto& be = vec.raw_elements()[start];
                        return add_edge_differential(be);
                }

                size_t mid = start + (end - start) / 2;
                ExtraEdgeGC left = add_edge_differential_recursive(start, mid);
                ExtraEdgeGC right = add_edge_differential_recursive(mid, end);
                left += right;
                return left;
        }

public: //debug
        void print() const {
                const auto& elems = vec.raw_elements();
                for (const auto& elem : elems) {
                        elem.getValue().print();
                        std::cout << "Coefficient: " << elem.getCoefficient() << "\n\n";
                }
        }

        void printFront() const {
                const auto& elems = vec.raw_elements();
                if (elems.size() == 0) {
                        cout << "0";
                }

                elems.front().getValue().print();
                std::cout << "Coefficient: " << elems.front().getCoefficient() << "\n\n";
        }

};