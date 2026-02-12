#pragma once

#include "VectorSpace/LinComb.hpp"
#include "VectorSpace/FiniteSubSpace.hpp"
#include "VectorSpace/BasisElement.hpp"
#include "VectorSpace/wiedemann_primitive_finder.hpp"

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
        using GraphType = Graph<N_VERTICES, N_EDGES, N_OUT_HAIR, N_IN_HAIR, c, d, fieldType>;
        using SplitGC = GC<N_VERTICES + 1, N_EDGES + 1, N_OUT_HAIR, N_IN_HAIR, c, d>;
        using ContGC = GC<N_VERTICES - 1, N_EDGES - 1, N_OUT_HAIR, N_IN_HAIR, c, d>;

        using ExtraEdgeGC = GC<N_VERTICES, N_EDGES + 1, N_OUT_HAIR, N_IN_HAIR, c, d>;

        using SplitGraphType = typename GraphType::SplitGraph;
        using ContGraphType = typename GraphType::ContGraph;
        using ExtraEdgeGraphType = typename GraphType::ExtraEdgeGraph;

        using L      = VectorSpace::LinComb<GraphType, fieldType>;
        using ContL  = VectorSpace::LinComb<ContGraphType, fieldType>;
        using SplitL = VectorSpace::LinComb<SplitGraphType, fieldType>;
        using ExtraEdgeL = VectorSpace::LinComb<ExtraEdgeGraphType, fieldType>;

		using Base      = BasisElement<GraphType, fieldType>;

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

        // construct from a LinComb
        explicit GC(VectorSpace::LinComb<GraphType, fieldType>& v) : vec(std::move(v))  {}

        explicit GC(std::vector<Base>&& elems)
                        : vec(std::move(elems)) {} 

        explicit GC(std::vector<Base>&& elems,
                        AssumeBasisOrderTag)                    
                : vec(std::move(elems), AssumeBasisOrderTag{})  
        {}

        // Access the underlying vector
        const VectorSpace::LinComb<GraphType, fieldType>& data() const {
                return vec;
        }
        
        VectorSpace::LinComb<GraphType, fieldType>& data_mutable() const {
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
                SplitL lin_comb = G.getValue().split_vertex_differential(G.getCoefficient());

                return SplitGC(lin_comb);
        }

        ExtraEdgeGC add_edge_differential() const {
                return add_edge_differential_recursive(0, vec.raw_elements().size());
        }

        static ExtraEdgeGC add_edge_differential(const BasisElement<GraphType, fieldType>& G) {
                ExtraEdgeL lin_comb = G.getValue().add_edge_differential(G.getCoefficient());


                return ExtraEdgeGC(lin_comb);
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

        void add_contractions_to_set(unordered_set<ContGraphType>& seenGraphs) {
                std::vector<BasisElement<ContGraphType, fieldType>> elems = d_contraction_without_sort();
                ContGC dThis(std::move(elems), AssumeBasisOrderTag{});
                
                dThis.standardize_all();

                for (auto b : dThis.data()) {	
						if (b.getCoefficient != fieldType{}) seenGraphs.insert(b.getValue());
                }
                return dThis;
        }
        
        
        unordered_map<ContGraphType, bigInt> contractions_count_map() const {
                std::vector<BasisElement<ContGraphType, fieldType>> elems = d_contraction_without_sort();
                ContGC dThis(std::move(elems), AssumeBasisOrderTag{});
                
                dThis.standardize_all();
				unordered_map<ContGraphType, bigInt> count_map;
				count_map.reserve(dThis.size());
				
                for (auto b : dThis.data()) {
                        count_map[b.getValue()]++;

                }
                return count_map;
        }
        
        
        
       ThisGC filter_on_odd_pairs(signedInt n_pairs) {
			vector<Base> filtered;
			
			for (const auto &be : data()) {
					if (be.getValue().odd_pairs() > n_pairs) {
							filtered.push_back(be);
					}
			}
			
			return ThisGC(filtered);		   
		}
        
                
        unordered_map<ContGraphType, bigInt> contractions_count_map_with_grading() const {
                std::vector<BasisElement<ContGraphType, fieldType>> elems = d_contraction_without_sort();
                
              
                ContGC dThis(std::move(elems), AssumeBasisOrderTag{});
                
                dThis.standardize_all();
				unordered_map<ContGraphType, bigInt> count_map;
				count_map.reserve(dThis.size());
				
                for (auto b : dThis.data()) {
                        count_map[b.getValue()]++;

                }
                return count_map;
        }
        

        ContGC d_contraction_with_recording_seen_gaphs(unordered_set<ContGraphType>& seenGraphs) {
				ContGC dThis = d_contraction();
			
				for (const auto& be :dThis.data()) {
					seenGraphs.emplace(be.getValue());
					
				}
				
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
        
        
        
        std::optional<ContGC> try_find_split_primitive2() {
				unordered_map<GraphType, unordered_set<ContGraphType>> contraction_cache;
				unordered_map<GraphType, SplitL> split_diff;
				
				contraction_cache.reserve(size());
				
				
				for (const auto& be : data()) {
						contraction_cache.emplace(be.getValue(), be.getValue().contraction_set());
						split_diff.emplace(be.getValue(), be.getValue().split_vertex_differential(be.getCoefficient()));
				}
				
				cout << "created contraction cache and split diff" << endl;
				
				MetaGraph split_meta_graph(split_diff);
				
				cout << "Created meta graph" << endl;
				split_diff.clear();
				
				unordered_map<ContGraphType, L> upper_split_diff;
				for (const auto& meta_edge : split_meta_graph) {						
						for (const ContGraphType& g : contraction_cache[meta_edge.x]) {
									if (!contraction_cache[meta_edge.y].contains(g)) { 
										continue;
									}
									//we are in the intersection
									
									if (!upper_split_diff.contains(g)) {
										upper_split_diff.emplace(g, g.split_vertex_differential(fieldType{1}));
									}
									break;
					
							}
				} 
				
				for (const GraphType& h: split_meta_graph.hair) {
						for (const ContGraphType& g: contraction_cache[h]) {
								upper_split_diff.emplace(g, g.split_vertex_differential(fieldType{1}));
						}
				}
				cout << "created upper split differential. size = " << upper_split_diff.size();
				contraction_cache.clear();
				
				
				//todo: also clear split_meta_graph
				
				
				VectorSpace::wiedemann_primitive_finder solver(upper_split_diff);
                       
       
                cout << "created solver" << endl;
                std::optional<ContL> primitive_optional = solver.find_primitive_or_empty(this -> data());
                
				
				cout << "solved "<< endl;
                return primitive_optional.transform([](ContL lin_comb) { 
                        return ContGC(lin_comb);
                });
				 
		}
        
       
		static std::pair<L, L> split_L_by_n_odd_vertex_pairs(const L& x, signedInt n) {
			L same;
			L rest;

			for (const auto& be : x) {
				if (be.getValue().n_odd_pairs() == n) same.append_in_basis_order(be);
				else rest.append_in_basis_order(be);
			}
			return {same, rest};
		}
			
		
		vector<L> graded_by_odd_vertex_pairs() const {
				constexpr int max_pairs = static_cast<int>(N_VERTICES)/2 + 1;	
				
				vector<L> split(max_pairs);
				
				
				for (const auto& be : data()) {
						signedInt n = be.getValue().n_odd_pairs();
						
						
						split[n].append_in_basis_order(be);
					
				}
				
				return split;
		}
			
	std::optional<ContGC> try_find_split_primitive_graded() const
	{
		std::cout << "using graded solver\n";

		auto contraction_counts = contractions_count_map();
		std::cout << "created contraction counts. total size: " << contraction_counts.size() << "\n";

		std::size_t map_size = 0;
		for (const auto& gamma : contraction_counts) {
			if (gamma.second > 1) ++map_size;
		}

		std::unordered_map<ContGraphType, L> full_coboundary_map;
		full_coboundary_map.reserve(map_size);

		for (const auto& gamma : contraction_counts) {
			if (gamma.second < 2) continue;
			full_coboundary_map.emplace(
				gamma.first,
				gamma.first.split_vertex_differential(fieldType{1})
			);
		}
		std::cout << "created full coboundary map!\n";

		// RHS split by grade
		std::vector<L> ass_graded = graded_by_odd_vertex_pairs();

		constexpr std::size_t max_pairs = static_cast<std::size_t>(N_VERTICES) / 2 + 1;

		// δ0 columns (grade i -> i)
		std::vector<std::unordered_map<ContGraphType, L>> same_to_same(max_pairs);

		// δ1 columns (grade i -> i+1), used only for propagation
		std::vector<std::unordered_map<ContGraphType, L>> up_maps(max_pairs);

		// full δ columns from grade i-1, available as auxiliary variables in grade i solve
		// constraint_maps[0] intentionally empty
		std::vector<std::unordered_map<ContGraphType, L>> constraint_maps(max_pairs);

		auto split_n_and_nplus1 = [](const L& x, signedInt n) -> std::pair<L, L> {
			L diag;
			L up;
			for (const auto& be : x) {
				const signedInt m = be.getValue().n_odd_pairs();
				if (m == n) {
					diag.append_in_basis_order(be);
				} else if (m == n + 1) {
					up.append_in_basis_order(be);
				} else {
					std::cout << "WARNING: δ produced grade " << m
							  << " from grade " << n << "\n";
				}
			}
			return {diag, up};
		};

		// Build the three maps
		for (const auto& [g, d_g] : full_coboundary_map) {
			const signedInt n = g.n_odd_pairs();
			if (n < 0 || static_cast<std::size_t>(n) >= max_pairs) {
				std::cout << "warning: n_odd_pairs() out of range: " << n << "\n";
				continue;
			}

			auto [diag, up] = split_n_and_nplus1(d_g, n);

			same_to_same[static_cast<std::size_t>(n)].emplace(g, std::move(diag));

			if (static_cast<std::size_t>(n + 1) < max_pairs) {
				up_maps[static_cast<std::size_t>(n)].emplace(g, std::move(up));

				// Constraint columns for the NEXT grade solve (i = n+1): store FULL δ(g)
				// Copy is intended (do not move).
				constraint_maps[static_cast<std::size_t>(n + 1)].emplace(g, d_g);
			}
		}

		std::cout << "built same_to_same / up_maps / constraint_maps\n";

		std::vector<ContL> primitive_by_grade(max_pairs);

		for (std::size_t i = 0; i < max_pairs; ++i) {
			std::cout << "-------------------------------------------------------------------\n";
			std::cout << "Trying to find primitive for " << i << "\n";

			// Constrained solve at grade i:
			// variables from V_i via same_to_same[i] (δ0),
			// plus auxiliary vars from V_{i-1} via constraint_maps[i] (full δ columns).
			VectorSpace::wiedemann_primitive_finder solver(same_to_same[i], constraint_maps[i]);

			auto prim_opt = solver.find_primitive_or_empty(ass_graded[i]);
			if (!prim_opt) {
				std::cout << "failed to find primitive\n";
				
				
				// try with full matrix
				
				std::cout<< "trying with full matrix" << endl;
				
				VectorSpace::wiedemann_primitive_finder full_solver(full_coboundary_map);
				
				
				auto primitive = full_solver.find_primitive_or_empty(data());
				
				if (!primitive) {
					std::cout << "Also failed with full matrix. " << std::endl;	
					return std::nullopt;
				}
				std::cout << "Worked with full solver. Must be a bug in the decomposition" << std::endl;
				
				return ContGC(*primitive);
				
			}

			primitive_by_grade[i] = std::move(*prim_opt);

			// Propagate upward using ONLY the +1 block: y_{i+1} <- y_{i+1} - δ1(x_i)
			if (i + 1 >= max_pairs) continue;

			L rest;
			for (const auto& be : primitive_by_grade[i]) {
				auto it = up_maps[i].find(be.getValue());
				if (it == up_maps[i].end()) continue;
				rest = rest.add_scaled(it->second, -be.getCoefficient());
			}
			ass_graded[i + 1] += rest;
		}

		ContL final_primitive;
		for (const auto& p : primitive_by_grade) final_primitive += p;

		return ContGC(final_primitive);
	}

                
        std::optional<ContGC> try_find_split_primitive() const {        
                unordered_map<ContGraphType, bigInt> contraction_counts = contractions_count_map();
                cout << "created contraction counts. total size: " << contraction_counts.size() <<endl;
                unordered_map<ContGraphType, L> coboundary_map;
			
				size_t map_size = 0;
				for (const auto& gamma : contraction_counts) {
						if (gamma.second > 1) map_size++;
				} 
				
				
				cout << "map_size = " << map_size <<  endl;
				coboundary_map.reserve(map_size);

                for (const auto& gamma : contraction_counts) {
						if (gamma.second < 2) continue;
                        coboundary_map.emplace(gamma.first, ContGC(gamma.first, AssumeBasisOrderTag{}).delta().data());
                }
                
                cout << "created couboundary map!" << endl;

                VectorSpace::wiedemann_primitive_finder solver(coboundary_map);
                       
       
                cout << "created solver" << endl;
                std::optional<ContL> primitive_optional = solver.find_primitive_or_empty(this -> data());
                
                if (!primitive_optional.has_value()) {
					//try with old dense solver
		
					cout << "trying with dense solver" << endl;

					VectorSpace::BoundaryFinder solver(coboundary_map);
						   
		   
					cout << "created solver" << endl;
					primitive_optional = solver.find_primitive_or_empty(this -> data());
				
					if (primitive_optional.has_value()) {
							cout << "sparse solver is bugged" << endl; 
					}
				}
                
				cout << "solved "<< endl;
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
        
        unordered_map<GraphType, SplitL> map_split_differential() {
			unordered_map<GraphType, SplitL> delta;
			for (const auto& be : vec) {
					delta.emplace(be.getValue(), be.getValue().split_vertex_differential(fieldType{1}));
			}
			return delta;
		}
		
		ThisGC filtered(unordered_set<GraphType>& filter) {
				vector<Base> filtered;
				
				for (const auto& be : vec) {
						if (!filter.contains(be.getValue())) continue;
						filtered.push_back(be);	
				}
	
				return ThisGC(std::move(filtered), AssumeBasisOrderTag{});
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
        
        
        
        static void add_to_even_splits_map(const GraphType& graph, unordered_map<SplitGraphType, L>& boundary_map,
											unordered_set<GraphType>& next_to_add) {
				
				SplitL delta_graph = graph.split_vertex_differential_even(fieldType{1});
				
				
				for (const auto& be : delta_graph) {
						if (!boundary_map.contains(be.getValue())) {
								auto contracted = SplitGC(be.getValue(), AssumeBasisOrderTag{})
                                                        .d_even_contraction_with_recording_seen_gaphs(next_to_add)
                                                        .data();
								
								boundary_map.emplace(be.getValue(), std::move(contracted));
								
						} 		
				}
		}
		
		        
        static void add_to_splits_map(const GraphType& graph, unordered_map<SplitGraphType, L>& boundary_map,
											unordered_set<GraphType>& next_to_add) {
				
				SplitL delta_graph = graph.split_vertex_differential(fieldType{1});
				
				
				for (const auto& be : delta_graph) {
						if (!boundary_map.contains(be.getValue())) {
								auto contracted = SplitGC(be.getValue(), AssumeBasisOrderTag{})
                                                        .d_contraction_with_recording_seen_gaphs(next_to_add)
                                                        .data();
								
								boundary_map.emplace(be.getValue(), std::move(contracted));
								
						} 
						
				}
			
		}
        
        std::optional<SplitGC> try_find_even_cont_primitive() const {
				unordered_map<SplitGraphType, L> boundary_map;
				unordered_set<GraphType> recently_added[2];
				unordered_set<GraphType> already_added;
						
				size_t max_depth = 4;
			
				for (const auto& be : data()) {
						recently_added[0].insert(be.getValue());
				}
				
				for (size_t i= 0; i< max_depth; ++i) {
						std::cout << "------------- depth: " << i << " ------------------"<<  std::endl;
					
						size_t i_mod2 = i%2;
						size_t i_plus_mod2 = (i+1)%2;
						
						
						for (const auto& graph : recently_added[i_mod2]) {
								if (already_added.contains(graph)) continue;
							
								add_to_even_splits_map(graph, boundary_map, recently_added[i_plus_mod2]);
								
						}
						
						already_added.insert(std::make_move_iterator(recently_added[i_mod2].begin()),
											 std::make_move_iterator(recently_added[i_mod2].end()));
											 							 
						recently_added[i_mod2].clear();
						
						VectorSpace::wiedemann_primitive_finder solver(boundary_map);
						
						
						auto primitive = solver.find_primitive_or_empty(data());
						
						
						if (primitive.has_value()) {
								return std::optional(SplitGC(*primitive));
						}
					
				}
				
				std::cout << "reach max depth " << max_depth << " did not find solution" << std::endl;
				return std::nullopt;
		}
		
		
		std::optional<SplitGC> try_find_cont_primitive() const {
			vector<SplitGraphType> empty;
			
			return try_find_cont_primitive(empty);
		
		}
		
		
		
		
		
	
		std::optional<SplitGC> try_find_cont_primitive(vector<SplitGraphType>& graps_to_excluce) const {
		
				unordered_map<SplitGraphType, L> boundary_map;
				
				for (const auto& g : graps_to_excluce) {
						boundary_map[g] = L{};
				}
	
				
				unordered_set<GraphType> recently_added[2];
				unordered_set<GraphType> already_added;
			
				size_t max_depth = 4;
			
				for (const auto& be : data()) {
						recently_added[0].insert(be.getValue());
				}
				
				for (size_t i= 0; i< max_depth; ++i) {
						std::cout << "------------- depth: " << i << " ------------------"<<  std::endl;
					
						size_t i_mod2 = i%2;
						size_t i_plus_mod2 = (i+1)%2;
						
						
						for (const auto& graph : recently_added[i_mod2]) {
								if (already_added.contains(graph)) continue;
							
								add_to_splits_map(graph, boundary_map, recently_added[i_plus_mod2]);
								
						}
						
						already_added.insert(std::make_move_iterator(recently_added[i_mod2].begin()),
											 std::make_move_iterator(recently_added[i_mod2].end()));
											 							 
						recently_added[i_mod2].clear();
						
						VectorSpace::wiedemann_primitive_finder solver(boundary_map);
						auto primitive = solver.find_primitive_or_empty(data());
						
						if (primitive.has_value()) {
								return std::optional(SplitGC(*primitive));
						}
					
				}
				
				std::cout << "reach max depth " << max_depth << " did not find solution" << std::endl;
				
				return std::nullopt;
				
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
        void print(std::ostream& out = std::cout) const {
                out<< "GC print of size : " << vec.size() << endl;
                out<< "Coefficient type : " << TypeName<fieldType>::name() << endl;
                const auto& elems = vec.raw_elements();
                for (const auto& elem : elems) {
                        elem.getValue().print(out);
                        out << "Coefficient: " << elem.getCoefficient() << "\n\n";
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
