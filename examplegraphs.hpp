#pragma once

#include "graph.hpp"
#include "GC.hpp"


using namespace std;

namespace graphaliases {

	// Alias for a graph with odd degree zero condition.
	template <Int N>
		using OddGraphdegZero = Graph<N, 2 * N - 2, 0, 0, 0, 1, fieldType>;

	template <Int N>
		using OddGCdegZero = GC<N, 2 * N - 2, 0, 0, 0, 1>;


	template <Int N>
		using OddLoopGraphType = Graph<N, N, 0, 0, 0, 1, fieldType>;


}

using namespace graphaliases;

template<Int N> 
 OddGraphdegZero<N + 1> wheel_graph() {
	array<Int, 4*N> arr{};

	OddGraphdegZero<N+1>  W(arr);

	//add spokes
	for (Int i = 0; i<N; ++i) {
		W.setEdge(i, 0, i+1);
	}
	//add rim
	for (Int i = 0; i <N-1; ++i) {
		W.setEdge(N+i, i+1, i+2); 

	}
	W.setEdge(2*N-1, 1, N); 

	return W;
}



template<Int N> 
OddGraphdegZero<N + 1> truss_graph() {
	array<Int, 4*N > arr{};

	OddGraphdegZero<N+1>  T(arr);
	// N + 1 vertices,
	// 2 *N edges
	Int edge_it = 0;
	for (Int vertex_it = 0; vertex_it < N-1; vertex_it++) {
		T.setEdge(edge_it++, vertex_it, vertex_it+1);
		T.setEdge(edge_it++, vertex_it, vertex_it+2);
	}

	T.setEdge(edge_it++, N-1, N);
	T.setEdge(edge_it, 0, N);

	return T;
}

template<Int N> 
OddLoopGraphType<N> loop_graph() {
	array<Int, 2*N> arr{};
	OddLoopGraphType<N>  loop(arr);

	for (Int i = 0; i < N-1; ++i) {
		loop.setEdge(i, i, i+1); 
	}

	loop.setEdge(N-1, N-1, 0); 
	return loop;
}

template<Int N> 
OddGraphdegZero<N + 1> V_graph_failed(vector<bool>& left_right_sequence) {
	Int k = N - 2 * left_right_sequence.size(); 

	OddGraphdegZero<N+1> V;
	Int e = 0;

	for (Int i = 1; i < k; i++) {
		V.setEdge(e++, 0, i);
		if (i != k-1) {
			V.setEdge(e++, i, i+1);
		}
	} 

	V.setEdge(e++, 0, k);

	Int root = k;
	Int left_root = 1;
	Int right_root = k-1;

	for (auto d : left_right_sequence) {
		V.setEdge(e++, root, root + 1);
		V.setEdge(e++, root, root + 2);

		if (d) {
			V.setEdge(e++, root, right_root);
			V.setEdge(e++, root+1 , right_root);
			
			right_root = root+ 1;
			
		} else {
			V.setEdge(e++, root, left_root);
			V.setEdge(e++, root+1, left_root);
			
			left_root = root+ 1;


		}
		root = root + 2; 
	}
	V.setEdge(e++, root , left_root);
	V.setEdge(e++, root , right_root);

	return V;
}

template<Int N> 
OddGraphdegZero<N + 1> V_graph(vector<bool>& left_right_sequence) {
	OddGraphdegZero<N+1> V; 
	Int e = 0;

	V.setEdge(e++, 0,1);
	V.setEdge(e++, 0,2);
	V.setEdge(e++, 0,3);
	V.setEdge(e++, 1,2);
	V.setEdge(e++, 2,3);


	Int left = 1; 
	Int center = 2;
	Int right = 3;

	Int next_vertex = 4;

	for (auto d : left_right_sequence) {
		V.setEdge(e++, center , next_vertex);
		center = next_vertex;
		next_vertex++;
		if (d) {
			V.setEdge(e++, center, right);
			V.setEdge(e++, right, next_vertex);
			V.setEdge(e++, center, next_vertex);
			right = next_vertex;
			next_vertex++;
		} else {
			V.setEdge(e++, center, left);
			V.setEdge(e++, left, next_vertex);
			V.setEdge(e++, center, next_vertex);
			left = next_vertex;
			next_vertex++;
		}
	}
	
	Int arc_pos = left;
	Int arc_end = right;

	while (next_vertex< N + 1) {
		V.setEdge(e++, arc_pos, next_vertex);
		V.setEdge(e++, center, next_vertex);

		arc_pos = next_vertex;
		next_vertex++;
	}

	V.setEdge(e++, arc_pos, right);
	return V;
}

template<typename GraphType>
void append_V_sequence_block(
	GraphType& graph,
	Int& edge_it,
	Int& next_vertex,
	Int& left,
	Int& center,
	Int& right,
	const vector<bool>& left_right_sequence
) {
	for (const bool d : left_right_sequence) {
		graph.setEdge(edge_it++, center, next_vertex);
		center = next_vertex;
		next_vertex++;
		if (d) {
			graph.setEdge(edge_it++, center, right);
			graph.setEdge(edge_it++, right, next_vertex);
			graph.setEdge(edge_it++, center, next_vertex);
			right = next_vertex;
			next_vertex++;
		} else {
			graph.setEdge(edge_it++, center, left);
			graph.setEdge(edge_it++, left, next_vertex);
			graph.setEdge(edge_it++, center, next_vertex);
			left = next_vertex;
			next_vertex++;
		}
	}
}

template<typename GraphType>
void extend_V_seed_block(
	GraphType& graph,
	Int& edge_it,
	Int& next_vertex,
	Int& left,
	Int& right
) {
	const Int new_0 = next_vertex++;
	const Int new_L = next_vertex++;
	const Int new_R = next_vertex++;

	graph.setEdge(edge_it++, new_0, new_L);
	graph.setEdge(edge_it++, new_0, new_R);
	graph.setEdge(edge_it++, 0, new_0);
	graph.setEdge(edge_it++, left, new_L);
	graph.setEdge(edge_it++, right, new_R);

	left = new_L;
	right = new_R;
}

inline vector<array<Int, 4>> all_V_seed_attachment_patterns(
	Int zero,
	Int left,
	Int center,
	Int right
) {
	const array<Int, 4> trivalent_vertices{zero, left, center, right};
	vector<array<Int, 4>> patterns;
	patterns.reserve(24);

	for (Int omitted_index = 0; omitted_index < 4; ++omitted_index) {
		array<Int, 3> attached_vertices{};
		Int attached_it = 0;
		for (Int vertex_index = 0; vertex_index < 4; ++vertex_index) {
			if (vertex_index == omitted_index) {
				continue;
			}
			attached_vertices[attached_it++] = trivalent_vertices[vertex_index];
		}

		for (Int i = 0; i < 3; ++i) {
			for (Int j = 0; j < 3; ++j) {
				if (j == i) {
					continue;
				}
				for (Int k = 0; k < 3; ++k) {
					if (k == i || k == j) {
						continue;
					}
					patterns.push_back({
						attached_vertices[i],
						attached_vertices[j],
						attached_vertices[k],
						trivalent_vertices[omitted_index]
					});
				}
			}
		}
	}

	return patterns;
}

template<typename GraphType>
void extend_V_seed_block_general(
	GraphType& graph,
	Int& edge_it,
	Int& next_vertex,
	Int attach_to_new_0,
	Int attach_to_new_L,
	Int attach_to_new_R,
	Int new_center,
	Int& left,
	Int& center,
	Int& right
) {
	const Int new_0 = next_vertex++;
	const Int new_L = next_vertex++;
	const Int new_R = next_vertex++;

	graph.setEdge(edge_it++, new_0, new_L);
	graph.setEdge(edge_it++, new_0, new_R);
	graph.setEdge(edge_it++, attach_to_new_0, new_0);
	graph.setEdge(edge_it++, attach_to_new_L, new_L);
	graph.setEdge(edge_it++, attach_to_new_R, new_R);

	left = new_L;
	center = new_center;
	right = new_R;
}

template<typename GraphType>
void extend_V_lcr_block(
	GraphType& graph,
	Int& edge_it,
	Int& next_vertex,
	Int& left,
	Int& center,
	Int& right
) {
	const Int new_left = next_vertex++;
	const Int new_center = next_vertex++;
	const Int new_right = next_vertex++;

	graph.setEdge(edge_it++, left, new_left);
	graph.setEdge(edge_it++, center, new_center);
	graph.setEdge(edge_it++, right, new_right);
	graph.setEdge(edge_it++, new_left, new_center);
	graph.setEdge(edge_it++, new_center, new_right);

	left = new_left;
	center = new_center;
	right = new_right;
}

template<Int N>
OddGraphdegZero<N + 1> iterated_V_graph(const vector<vector<bool>>& sequence_blocks) {
	const std::size_t number_of_blocks = sequence_blocks.size();
	const std::size_t total_sequence_length = [&]() {
		std::size_t total = 0;
		for (const auto& block : sequence_blocks) {
			total += block.size();
		}
		return total;
	}();

	const std::size_t expected_vertices = 4 + 2 * total_sequence_length
		+ (number_of_blocks == 0 ? 0 : 3 * (number_of_blocks - 1));
	assert(expected_vertices == static_cast<std::size_t>(N + 1));

	OddGraphdegZero<N + 1> V;
	Int e = 0;

	V.setEdge(e++, 0, 1);
	V.setEdge(e++, 0, 2);
	V.setEdge(e++, 0, 3);
	V.setEdge(e++, 1, 2);
	V.setEdge(e++, 2, 3);

	Int left = 1;
	Int center = 2;
	Int right = 3;
	Int next_vertex = 4;

	for (std::size_t block_index = 0; block_index < number_of_blocks; ++block_index) {
		append_V_sequence_block(V, e, next_vertex, left, center, right, sequence_blocks[block_index]);
		V.setEdge(e++, left, right);

		if (block_index + 1 < number_of_blocks) {
			extend_V_seed_block(V, e, next_vertex, left, right);
		}
	}

	assert(next_vertex == static_cast<Int>(N + 1));
	assert(e == static_cast<Int>(2 * N));
	return V;
}

template<Int N>
OddGraphdegZero<N + 1> iterated_V_graph_lcr_extension(const vector<vector<bool>>& sequence_blocks) {
	const std::size_t number_of_blocks = sequence_blocks.size();
	const std::size_t total_sequence_length = [&]() {
		std::size_t total = 0;
		for (const auto& block : sequence_blocks) {
			total += block.size();
		}
		return total;
	}();

	const std::size_t expected_vertices = 4 + 2 * total_sequence_length
		+ (number_of_blocks == 0 ? 0 : 3 * (number_of_blocks - 1));
	assert(expected_vertices == static_cast<std::size_t>(N + 1));

	OddGraphdegZero<N + 1> V;
	Int e = 0;

	V.setEdge(e++, 0, 1);
	V.setEdge(e++, 0, 2);
	V.setEdge(e++, 0, 3);
	V.setEdge(e++, 1, 2);
	V.setEdge(e++, 2, 3);

	Int left = 1;
	Int center = 2;
	Int right = 3;
	Int next_vertex = 4;

	for (std::size_t block_index = 0; block_index < number_of_blocks; ++block_index) {
		append_V_sequence_block(V, e, next_vertex, left, center, right, sequence_blocks[block_index]);
		V.setEdge(e++, left, right);

		if (block_index + 1 < number_of_blocks) {
			extend_V_lcr_block(V, e, next_vertex, left, center, right);
		}
	}

	assert(next_vertex == static_cast<Int>(N + 1));
	assert(e == static_cast<Int>(2 * N));
	return V;
}


template<Int N> 
typename OddGraphdegZero<N + 1>::ContGraph V_con_graph(vector<bool>& left_right_sequence, Int separation_point) {
	typename OddGraphdegZero<N + 1>::ContGraph  V; 
	Int e = 0;

	V.setEdge(e++, 0,1);
	V.setEdge(e++, 0,2);
	V.setEdge(e++, 0,3);
	V.setEdge(e++, 1,2);
	V.setEdge(e++, 2,3);

	Int left = 1; 
	Int center = 2;
	Int right = 3;

	Int edge_to_contract;

	Int next_vertex = 4;

	for (Int i =0 ; i< left_right_sequence.size(); i++) {
		if (i != separation_point) {
			V.setEdge(e++, center , next_vertex);
			center = next_vertex;
			next_vertex++;
		}
		if (left_right_sequence[i]) {
			V.setEdge(e++, center, right);
			V.setEdge(e++, right, next_vertex);
			V.setEdge(e++, center, next_vertex);
			right = next_vertex;
			next_vertex++;
		} else {
			V.setEdge(e++, center, left);
			V.setEdge(e++, left, next_vertex);
			V.setEdge(e++, center, next_vertex);
			left = next_vertex;
			next_vertex++;
		}
	}
	
	Int arc_pos = left;
	Int arc_end = right;

	while (next_vertex< N + 1) {
		V.setEdge(e++, arc_pos, next_vertex);
		V.setEdge(e++, center, next_vertex);

		arc_pos = next_vertex;
		next_vertex++;
	}

	V.setEdge(e++, arc_pos, right);
	return V;
}



template<Int N> 
typename OddGraphdegZero<N + 1>::SplitGraph U_graph(vector<bool>& left_right_sequence) {
	typename OddGraphdegZero<N + 1>::SplitGraph  U; 
	Int e = 0;

	U.setEdge(e++, 0,1);
	U.setEdge(e++, 0,2);
	U.setEdge(e++, 0,3);
	U.setEdge(e++, 1,2);
	U.setEdge(e++, 2,3);

	Int left = 1; 
	Int center = 2;
	Int right = 3;

	Int next_vertex = 4;

	for (auto d : left_right_sequence) {
		U.setEdge(e++, center , next_vertex);
		center = next_vertex;
		next_vertex++;
		if (d) {
			U.setEdge(e++, center, right);
			U.setEdge(e++, right, next_vertex);
			U.setEdge(e++, center, next_vertex);
			right = next_vertex;
			next_vertex++;
		} else {
			U.setEdge(e++, center, left);
			U.setEdge(e++, left, next_vertex);
			U.setEdge(e++, center, next_vertex);
			left = next_vertex;
			next_vertex++;
		}
	}

	U.setEdge(e++, center, next_vertex);
	center = next_vertex;
	next_vertex++;
	
	Int arc_pos = left;
	Int arc_end = right;

	while (next_vertex< N + 2) {
		U.setEdge(e++, arc_pos, next_vertex);
		U.setEdge(e++, center, next_vertex);

		arc_pos = next_vertex;
		next_vertex++;
	}

	U.setEdge(e++, arc_pos, right);
	return U;
}

inline OddGraphdegZero<10> python_G9_1_graph() {
	OddGraphdegZero<10> g;
	g.setEdge(0, 4, 5);
	g.setEdge(1, 3, 4);
	g.setEdge(2, 2, 3);
	g.setEdge(3, 1, 2);
	g.setEdge(4, 0, 1);
	g.setEdge(5, 0, 5);
	g.setEdge(6, 5, 7);
	g.setEdge(7, 4, 7);
	g.setEdge(8, 3, 7);
	g.setEdge(9, 3, 8);
	g.setEdge(10, 2, 8);
	g.setEdge(11, 2, 9);
	g.setEdge(12, 1, 9);
	g.setEdge(13, 1, 6);
	g.setEdge(14, 0, 6);
	g.setEdge(15, 7, 8);
	g.setEdge(16, 8, 9);
	g.setEdge(17, 6, 9);
	return g;
}

inline OddGraphdegZero<10> python_G9_2_graph() {
	OddGraphdegZero<10> g;
	g.setEdge(0, 0, 5);
	g.setEdge(1, 4, 5);
	g.setEdge(2, 3, 4);
	g.setEdge(3, 2, 3);
	g.setEdge(4, 1, 2);
	g.setEdge(5, 0, 1);
	g.setEdge(6, 0, 7);
	g.setEdge(7, 5, 7);
	g.setEdge(8, 4, 7);
	g.setEdge(9, 4, 8);
	g.setEdge(10, 3, 8);
	g.setEdge(11, 3, 9);
	g.setEdge(12, 2, 9);
	g.setEdge(13, 1, 6);
	g.setEdge(14, 0, 6);
	g.setEdge(15, 7, 8);
	g.setEdge(16, 8, 9);
	g.setEdge(17, 6, 9);
	return g;
}

inline OddGraphdegZero<10> python_G9_3_graph() {
	OddGraphdegZero<10> g;
	g.setEdge(0, 0, 5);
	g.setEdge(1, 4, 5);
	g.setEdge(2, 3, 4);
	g.setEdge(3, 2, 3);
	g.setEdge(4, 1, 2);
	g.setEdge(5, 0, 1);
	g.setEdge(6, 0, 7);
	g.setEdge(7, 5, 7);
	g.setEdge(8, 4, 7);
	g.setEdge(9, 4, 8);
	g.setEdge(10, 3, 8);
	g.setEdge(11, 3, 6);
	g.setEdge(12, 2, 6);
	g.setEdge(13, 1, 9);
	g.setEdge(14, 0, 9);
	g.setEdge(15, 7, 8);
	g.setEdge(16, 8, 9);
	g.setEdge(17, 6, 9);
	return g;
}

inline OddGraphdegZero<10> python_G9_4_graph() {
	OddGraphdegZero<10> g;
	g.setEdge(0, 0, 5);
	g.setEdge(1, 4, 5);
	g.setEdge(2, 3, 4);
	g.setEdge(3, 2, 3);
	g.setEdge(4, 1, 2);
	g.setEdge(5, 0, 1);
	g.setEdge(6, 0, 7);
	g.setEdge(7, 5, 7);
	g.setEdge(8, 4, 7);
	g.setEdge(9, 4, 8);
	g.setEdge(10, 3, 8);
	g.setEdge(11, 2, 6);
	g.setEdge(12, 1, 6);
	g.setEdge(13, 1, 9);
	g.setEdge(14, 0, 9);
	g.setEdge(15, 7, 8);
	g.setEdge(16, 8, 9);
	g.setEdge(17, 6, 9);
	return g;
}

inline OddGCdegZero<10> python_G9_class() {
	using G = OddGraphdegZero<10>;
	using B = BasisElement<G, fieldType>;
	std::vector<B> elems;
	elems.reserve(4);
	elems.emplace_back(python_G9_1_graph(), fieldType{-2});
	elems.emplace_back(python_G9_2_graph(), fieldType{-2});
	elems.emplace_back(python_G9_3_graph(), fieldType{2});
	elems.emplace_back(python_G9_4_graph(), fieldType{2});
	return OddGCdegZero<10>(std::move(elems));
}
