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
OddGraphdegZero<N+1> wheel_graph() {
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
