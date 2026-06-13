#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "examplegraphs.hpp"
#include "OCGraph_test.hpp"
#include "GraphStandardizer.hpp"

template <typename GraphType>
GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_> make_graph_isomorphism_from_vertex_perm(
	const GraphType& graph,
	const std::array<Int, GraphType::N_VERTICES_>& vertex_perm
) {
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	IsoType iso;
	iso.vertex_permutation_data() = vertex_perm;

	for (Int source_edge = 0; source_edge < GraphType::N_EDGES_; ++source_edge) {
		auto [u, v] = graph.getEdge(source_edge);
		u = vertex_perm[u];
		v = vertex_perm[v];

		const Int target_edge = graph.find_edge_index(u, v);
		iso.edge_permutation_data()[source_edge] = target_edge;

		const auto [tu, tv] = graph.getEdge(target_edge);
		iso.edge_flip_data()[source_edge] = (tu != u || tv != v);
	}

	iso.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();

	return iso;
}

template <Int N>
bool check_wheel_ocgraph_delta_e_isomorphism_equivariance(const char* label) {
	using GraphType = OddGraphdegZero<N + 1>;
	using DirectionType = GraphDirections<GraphType>;
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	const GraphType source = wheel_graph<N>();

	DirectionType directions;
	directions.fill(true);

	std::vector<IsoType> isomorphisms;

	std::array<Int, GraphType::N_VERTICES_> rotation{};
	rotation[0] = 0;
	for (Int v = 1; v < GraphType::N_VERTICES_; ++v) {
		rotation[v] = (v == N) ? 1 : static_cast<Int>(v + 1);
	}
	isomorphisms.push_back(make_graph_isomorphism_from_vertex_perm(source, rotation));

	std::array<Int, GraphType::N_VERTICES_> reflection{};
	reflection[0] = 0;
	reflection[1] = 1;
	for (Int v = 2; v < GraphType::N_VERTICES_; ++v) {
		reflection[v] = static_cast<Int>(GraphType::N_VERTICES_ + 1 - v);
	}
	isomorphisms.push_back(make_graph_isomorphism_from_vertex_perm(source, reflection));

	return check_ocgraph_delta_e_isomorphism_equivariance_on_data(label, source, directions, isomorphisms);
}

template <Int N, signedInt C, signedInt D>
Graph<N + 1, 2 * N, 0, 0, C, D, fieldType> generic_wheel_graph() {
	using GraphType = Graph<N + 1, 2 * N, 0, 0, C, D, fieldType>;

	std::array<Int, 4 * N> arr{};
	GraphType wheel(arr);

	for (Int i = 0; i < N; ++i) {
		wheel.setEdge(i, 0, i + 1);
	}
	for (Int i = 0; i < N - 1; ++i) {
		wheel.setEdge(N + i, i + 1, i + 2);
	}
	wheel.setEdge(2 * N - 1, 1, N);

	return wheel;
}

template <Int N, signedInt C, signedInt D>
bool check_wheel_ocgc_delta_v_squares_to_zero(const char* label) {
	using GraphType = Graph<N + 1, 2 * N, 0, 0, C, D, fieldType>;
	using DirectionType = GraphDirections<GraphType>;
	using OCGraphType = OCGraph<GraphType>;
	using OCGCType = OCGC<GraphType>;

	const GraphType source = generic_wheel_graph<N, C, D>();

	DirectionType d0;
	d0.fill(true);

	DirectionType d1;
	d1.fill(true);
	d1[0] = false;
	d1[GraphType::SIZE / 2] = false;

	typename OCGraphType::DirectionComb directions;
	directions.append_in_basis_order(d0, fieldType{1});
	directions.append_in_basis_order(d1, fieldType{2});

	OCGCType ocgc(OCGraphType(source, directions));
	auto delta1 = ocgc.delta_v();
	delta1.standardize_and_sort();
	auto delta2 = delta1.delta_v();
	delta2.standardize_and_sort();

	const bool ok = delta2.size() == 0;
	std::cout << label << ": wheel ocgc delta_v^2 -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <Int N, signedInt C, signedInt D>
bool check_wheel_ocgc_delta_anticommute(const char* label) {
	using GraphType = Graph<N + 1, 2 * N, 0, 0, C, D, fieldType>;
	using DirectionType = GraphDirections<GraphType>;
	using OCGraphType = OCGraph<GraphType>;
	using OCGCType = OCGC<GraphType>;

	const GraphType source = generic_wheel_graph<N, C, D>();

	DirectionType d0;
	d0.fill(true);

	DirectionType d1;
	d1.fill(true);
	d1[0] = false;
	d1[GraphType::SIZE / 2] = false;

	typename OCGraphType::DirectionComb directions;
	directions.append_in_basis_order(d0, fieldType{1});
	directions.append_in_basis_order(d1, fieldType{2});

	OCGCType ocgc(OCGraphType(source, directions));

	auto dve = ocgc.delta_v();
	dve.standardize_and_sort();
	dve = dve.delta_e();
	dve.standardize_and_sort();

	auto dev = ocgc.delta_e();
	dev.standardize_and_sort();
	auto dev_v = dev.delta_v();
	dev_v.standardize_and_sort();

	dve += dev_v;
	dve.standardize_and_sort();

	const bool ok = dve.size() == 0;
	std::cout << label << ": wheel ocgc delta_v delta_e + delta_e delta_v -> "
		  << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <Int N, signedInt C = 0, signedInt D = 1>
bool check_full_unsorted_splits_double(const char* label) {
	using GraphType = Graph<N + 1, 2 * N, 0, 0, C, D, fieldType>;

	const GraphType source = generic_wheel_graph<N, C, D>();

	auto ordinary = source.unsorted_splits(fieldType{1});
	auto full = source.unsorted_full_splits(fieldType{1});

	const bigInt ordinary_raw_size = ordinary.size();
	const bigInt full_raw_size = full.size();
	const bool raw_size_ok = full_raw_size == 2 * ordinary_raw_size;

	ordinary.standardize_and_sort();
	full.standardize_and_sort();

	const auto doubled = ordinary * fieldType{2};
	const bool standardized_ok = full == doubled;
	const bool ok = raw_size_ok && standardized_ok;

	std::cout << label
		<< ": raw=" << ordinary_raw_size
		<< " full_raw=" << full_raw_size
		<< " -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

template <typename GraphType>
std::vector<GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>> sample_graph_isomorphisms(
	const GraphType&
) {
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	std::vector<IsoType> isomorphisms;

	IsoType a;
	for (Int v = 0; v < GraphType::N_VERTICES_; ++v) {
		a.vertex_permutation_data()[v] = static_cast<Int>((v + 1) % GraphType::N_VERTICES_);
	}
	for (Int e = 0; e < GraphType::N_EDGES_; ++e) {
		a.edge_permutation_data()[e] = static_cast<Int>((e + 1) % GraphType::N_EDGES_);
		a.edge_flip_data()[e] = (e % 2) == 0;
	}
	a.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();
	isomorphisms.push_back(a);

	IsoType b;
	for (Int v = 0; v < GraphType::N_VERTICES_; ++v) {
		b.vertex_permutation_data()[v] = static_cast<Int>((GraphType::N_VERTICES_ - 1 - v));
	}
	for (Int e = 0; e < GraphType::N_EDGES_; ++e) {
		b.edge_permutation_data()[e] = static_cast<Int>((GraphType::N_EDGES_ - 1 - e));
		b.edge_flip_data()[e] = (e % 3) == 1;
	}
	b.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();
	isomorphisms.push_back(b);

	IsoType c;
	for (Int v = 0; v < GraphType::N_VERTICES_; ++v) {
		c.vertex_permutation_data()[v] = static_cast<Int>((v + 2) % GraphType::N_VERTICES_);
	}
	for (Int e = 0; e < GraphType::N_EDGES_; ++e) {
		c.edge_permutation_data()[e] = static_cast<Int>((e + 3) % GraphType::N_EDGES_);
		c.edge_flip_data()[e] = (e % 2) == 1;
	}
	c.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();
	isomorphisms.push_back(c);

	return isomorphisms;
}

template <Int N>
GraphIsomorphism<OddGraphdegZero<N + 1>::N_VERTICES_, OddGraphdegZero<N + 1>::N_EDGES_>
wheel_swap_spoke_and_rim_labels_isomorphism() {
	using GraphType = OddGraphdegZero<N + 1>;
	using IsoType = GraphIsomorphism<GraphType::N_VERTICES_, GraphType::N_EDGES_>;

	IsoType iso;
	for (Int edge = 0; edge < N; ++edge) {
		iso.edge_permutation_data()[edge] = static_cast<Int>(N + edge);
		iso.edge_permutation_data()[N + edge] = edge;
	}
	iso.template compute_signs<GraphType::N_OUT_HAIR_, GraphType::N_IN_HAIR_, GraphType::C_, GraphType::D_, fieldType>();
	return iso;
}

template <Int N>
std::vector<GraphIsomorphism<OddGraphdegZero<N + 1>::N_VERTICES_, OddGraphdegZero<N + 1>::N_EDGES_>>
wheel_OCG_F0_test_isomorphisms() {
	using GraphType = OddGraphdegZero<N + 1>;

	auto isomorphisms = sample_graph_isomorphisms<GraphType>(wheel_graph<N>());
	isomorphisms.push_back(wheel_swap_spoke_and_rim_labels_isomorphism<N>());
	return isomorphisms;
}

template <Int N>
bool check_wheel_OCG_F0_exact_equivariance(const char* label) {
	using GraphType = OddGraphdegZero<N + 1>;
	const GraphType source = wheel_graph<N>();
	return check_OCG_F0_exact_equivariance_on_graph(label, source, wheel_OCG_F0_test_isomorphisms<N>());
}

template <typename GraphType>
VectorSpace::LinComb<GraphType, fieldType> standardize4_and_sort(
	std::vector<BasisElement<GraphType, fieldType>> elems
) {
	GraphStandardizer<
		GraphType::N_VERTICES_,
		GraphType::N_EDGES_,
		GraphType::N_OUT_HAIR_,
		GraphType::N_IN_HAIR_,
		GraphType::C_,
		GraphType::D_,
		fieldType
	> standardizer;

	for (auto& be : elems) {
		be = standardizer.standardize4(be);
	}

	VectorSpace::LinComb<GraphType, fieldType> result(std::move(elems), AssumeBasisOrderTag{});
	result.sort_elements();
	return result;
}

template <typename GraphType>
VectorSpace::LinComb<typename GraphType::SplitGraph, fieldType> split_standardize4(
	const GraphType& graph,
	fieldType coeff = fieldType{1}
) {
	using SplitGraph = typename GraphType::SplitGraph;
	const auto raw = graph.unsorted_splits(coeff).raw_elements();
	std::vector<BasisElement<SplitGraph, fieldType>> elems(raw.begin(), raw.end());
	return standardize4_and_sort<SplitGraph>(std::move(elems));
}

template <typename SplitGraph>
VectorSpace::LinComb<typename SplitGraph::ContGraph, fieldType> contract_standardize4(
	const SplitGraph& graph,
	fieldType coeff = fieldType{1}
) {
	using ContGraph = typename SplitGraph::ContGraph;
	std::vector<BasisElement<ContGraph, fieldType>> elems;
	elems.reserve(SplitGraph::N_EDGES_);
	for (Int edge = 0; edge < SplitGraph::N_EDGES_; ++edge) {
		auto be = graph.contract_edge(edge, coeff);
		if (be.getCoefficient() != fieldType{}) {
			elems.push_back(std::move(be));
		}
	}
	return standardize4_and_sort<ContGraph>(std::move(elems));
}

template <typename GraphType>
bigInt automorphism_group_size4_gc_test(const GraphType& input_graph) {
	using Standardizer = GraphStandardizer<
		GraphType::N_VERTICES_,
		GraphType::N_EDGES_,
		GraphType::N_OUT_HAIR_,
		GraphType::N_IN_HAIR_,
		GraphType::C_,
		GraphType::D_,
		fieldType
	>;
	using IsomorphismType = typename Standardizer::IsomorphismType;

	Standardizer standardizer;
	auto [attempts, valid_attempts] = standardizer.create_final_attempts4(input_graph);

	std::vector<IsomorphismType> minimizers;
	typename GraphType::Basis best_basis;
	bool have_best = false;

	for (const std::size_t attempt_index : valid_attempts) {
		typename GraphType::ThisGraph graph;
		IsomorphismType iso;
		const signedInt sign = graph.assignPermutedDirectedSortedEdgesWithIsomorphism(
			input_graph,
			attempts[attempt_index].create_vertex_permutation(),
			iso
		);
		if (sign == 0) {
			continue;
		}

		typename GraphType::Basis attempt_basis(std::move(graph), fieldType{sign});
		if (!have_best) {
			best_basis = attempt_basis;
			minimizers.clear();
			minimizers.push_back(iso);
			have_best = true;
			continue;
		}

		const signedInt comparison = best_basis.compare(attempt_basis);
		if (comparison < 0) {
			best_basis = attempt_basis;
			minimizers.clear();
			minimizers.push_back(iso);
		} else if (comparison == 0) {
			minimizers.push_back(iso);
		}
	}

	return minimizers.size();
}

template <typename GCType>
GCType split_contract_step_gc_test(const GCType& gamma) {
	using GraphType = typename GCType::GraphType;
	using SplitGraph = typename GCType::SplitGraphType;
	using ContGraph = typename SplitGraph::ContGraph;

	std::vector<BasisElement<SplitGraph, fieldType>> split_terms;
	for (const auto& be : gamma.data()) {
		const auto raw = be.getValue().unsorted_splits(be.getCoefficient()).raw_elements();
		split_terms.insert(split_terms.end(), raw.begin(), raw.end());
	}
	auto split = standardize4_and_sort<SplitGraph>(std::move(split_terms));

	std::vector<BasisElement<ContGraph, fieldType>> contract_terms;
	for (const auto& be : split) {
		for (Int edge = 0; edge < SplitGraph::N_EDGES_; ++edge) {
			auto cont = be.getValue().contract_edge(edge, be.getCoefficient());
			if (cont.getCoefficient() != fieldType{}) {
				contract_terms.push_back(std::move(cont));
			}
		}
	}
	auto contracted = standardize4_and_sort<ContGraph>(std::move(contract_terms));
	return GCType(contracted);
}

template <Int N>
bool check_scaled_transpose_matches_contraction_rounds4(const char* label) {
	using GCType = OddGCdegZero<N + 1>;
	using GraphType = typename GCType::GraphType;
	using SplitGraph = typename GCType::SplitGraphType;

	typename GraphType::Basis wheel_basis(wheel_graph<N>(), fieldType{1});
	{
		GraphStandardizer<
			GraphType::N_VERTICES_,
			GraphType::N_EDGES_,
			GraphType::N_OUT_HAIR_,
			GraphType::N_IN_HAIR_,
			GraphType::C_,
			GraphType::D_,
			fieldType
		> standardizer;
		wheel_basis = standardizer.standardize4(wheel_basis);
	}
	GCType current(wheel_basis.getValue(), AssumeBasisOrderTag{});

	std::unordered_set<GraphType> seen_y_set;
	std::unordered_set<SplitGraph> seen_x_set;
	seen_y_set.insert(wheel_basis.getValue());

	const int rounds = 4;
	for (int i = 0; i <= rounds; ++i) {
		std::vector<BasisElement<SplitGraph, fieldType>> split_terms;
		for (const auto& be : current.data()) {
			const auto raw = be.getValue().unsorted_splits(be.getCoefficient()).raw_elements();
			split_terms.insert(split_terms.end(), raw.begin(), raw.end());
		}
		auto split = standardize4_and_sort<SplitGraph>(std::move(split_terms));
		for (const auto& be : split) {
			seen_x_set.insert(be.getValue());
		}
		if (i == rounds) {
			break;
		}
		current = split_contract_step_gc_test(current);
		for (const auto& be : current.data()) {
			seen_y_set.insert(be.getValue());
		}
	}

	std::vector<GraphType> y_basis(seen_y_set.begin(), seen_y_set.end());
	std::vector<SplitGraph> x_basis(seen_x_set.begin(), seen_x_set.end());
	std::sort(y_basis.begin(), y_basis.end());
	std::sort(x_basis.begin(), x_basis.end());

	std::unordered_map<GraphType, std::size_t> y_index;
	std::unordered_map<SplitGraph, std::size_t> x_index;
	for (std::size_t i = 0; i < y_basis.size(); ++i) {
		y_index.emplace(y_basis[i], i);
	}
	for (std::size_t i = 0; i < x_basis.size(); ++i) {
		x_index.emplace(x_basis[i], i);
	}

	std::vector<std::vector<std::pair<std::size_t, fieldType>>> split_columns;
	split_columns.reserve(y_basis.size());
	for (const auto& y_graph : y_basis) {
		std::vector<std::pair<std::size_t, fieldType>> col;
		auto split = split_standardize4<GraphType>(y_graph);
		col.reserve(split.size());
		for (const auto& be : split) {
			col.emplace_back(x_index.at(be.getValue()), be.getCoefficient());
		}
		split_columns.push_back(std::move(col));
	}

	std::vector<std::vector<fieldType>> contraction_matrix(
		y_basis.size(),
		std::vector<fieldType>(x_basis.size(), fieldType{})
	);
	for (std::size_t i = 0; i < x_basis.size(); ++i) {
		auto contracted = contract_standardize4<SplitGraph>(x_basis[i]);
		for (const auto& be : contracted) {
			auto it = y_index.find(be.getValue());
			if (it != y_index.end()) {
				contraction_matrix[it->second][i] = be.getCoefficient();
			}
		}
	}

	std::vector<fieldType> y_aut(y_basis.size(), fieldType{});
	std::vector<fieldType> x_aut(x_basis.size(), fieldType{});
	for (std::size_t j = 0; j < y_basis.size(); ++j) {
		y_aut[j] = fieldType{automorphism_group_size4_gc_test(y_basis[j])};
	}
	for (std::size_t i = 0; i < x_basis.size(); ++i) {
		x_aut[i] = fieldType{automorphism_group_size4_gc_test(x_basis[i])};
	}

	std::size_t mismatches = 0;
	for (std::size_t i = 0; i < x_basis.size(); ++i) {
		for (std::size_t j = 0; j < y_basis.size(); ++j) {
			fieldType scaled_entry{};
			for (const auto& [row, coeff] : split_columns[j]) {
				if (row == i) {
					scaled_entry = coeff * x_aut[i] / y_aut[j];
					break;
				}
			}
			if (scaled_entry != contraction_matrix[j][i]) {
				++mismatches;
			}
		}
	}

	const bool ok = mismatches == 0;
	std::cout << label << ": W7 rounds=4 scaled transpose vs d_contraction -> "
	          << (ok ? "ok" : "failed")
	          << " (mismatches=" << mismatches << ")\n";
	return ok;
}

bool check_w9_OCG_F0_unique_direction_data(const char* label) {
	using GraphType = OddGraphdegZero<10>;
	const GraphType source = wheel_graph<9>();
	return check_OCG_F0_unique_direction_on_graph(label, source, wheel_OCG_F0_test_isomorphisms<9>());
}

template <Int N>
bool check_wheel_OCG_F0_delta_homotopy_equivalence(const char* label) {
	using GraphType = OddGraphdegZero<N + 1>;
	const GraphType source = wheel_graph<N>();
	return check_OCG_F0_delta_homotopy_equivalence_on_graph(label, source);
}

std::vector<std::vector<bool>> left_right_sequences_all(int length) {
	std::vector<std::vector<bool>> result;
	if (length < 0) {
		return result;
	}

	const int count = 1 << length;
	result.reserve(count);
	for (int mask = 0; mask < count; ++mask) {
		std::vector<bool> sequence;
		sequence.reserve(length);
		for (int bit = 0; bit < length; ++bit) {
			sequence.push_back(((mask >> bit) & 1) != 0);
		}
		result.push_back(std::move(sequence));
	}
	return result;
}

template <Int N>
bool check_v_graph_sort_matches_quick(const char* label) {
	using GraphType = OddGraphdegZero<N + 1>;

	const int maximal_v_sequence_length = (N - 3) / 2;
	const auto sequences = left_right_sequences_all(maximal_v_sequence_length);

	std::size_t checked = 0;
	for (const auto& sequence : sequences) {
		const GraphType source = V_graph<N>(const_cast<std::vector<bool>&>(sequence));

		GraphType quick_sorted = source;
		GraphType radix_sorted = source;

		const signedInt quick_sign = quick_sorted.directEdges() * quick_sorted.sortEdgesQuick();
		const signedInt radix_sign = radix_sorted.directAndSortEdges();

		if (quick_sorted != radix_sorted || quick_sign != radix_sign) {
			std::cout << label << ": mismatch on sequence index " << checked << '\n';
			std::cout << "quick sign = " << quick_sign << ", radix sign = " << radix_sign << '\n';
			return false;
		}
		++checked;
	}

	std::cout << label << ": checked " << checked << " graphs -> ok\n";
	return true;
}

template <Int N>
bool check_v_graph_permuted_sort_matches_old_path(const char* label) {
	using GraphType = OddGraphdegZero<N + 1>;
	using PermType = Permutation<GraphType::N_VERTICES_>;

	const int maximal_v_sequence_length = (N - 3) / 2;
	const auto sequences = left_right_sequences_all(maximal_v_sequence_length);

	std::array<PermType, 3> perms{};
	for (Int v = 0; v < GraphType::N_VERTICES_; ++v) {
		perms[0].p[v] = v;
	}
	perms[1].p[0] = 0;
	for (Int v = 1; v < GraphType::N_VERTICES_; ++v) {
		perms[1].p[v] = (v == N) ? 1 : static_cast<Int>(v + 1);
	}
	perms[2].p[0] = 0;
	perms[2].p[1] = 1;
	for (Int v = 2; v < GraphType::N_VERTICES_; ++v) {
		perms[2].p[v] = static_cast<Int>(GraphType::N_VERTICES_ + 1 - v);
	}

	std::size_t checked = 0;
	for (const auto& sequence : sequences) {
		const GraphType source = V_graph<N>(const_cast<std::vector<bool>&>(sequence));

		for (const auto& perm : perms) {
			GraphType old_path = source;
			GraphType fused_path;

			const signedInt old_sign =
				old_path.permuteVertices(perm) *
				old_path.directEdges() *
				old_path.sortEdgesQuick();
			const signedInt fused_sign =
				fused_path.assignPermutedDirectedSortedEdges(source, perm);

			if (old_path != fused_path || old_sign != fused_sign) {
				std::cout << label << ": mismatch on graph index " << checked << '\n';
				std::cout << "old sign = " << old_sign << ", fused sign = " << fused_sign << '\n';
				return false;
			}
			++checked;
		}
	}

	std::cout << label << ": checked " << checked << " permuted graphs -> ok\n";
	return true;
}

template <typename GCType>
bool check_odd_even_contraction_split(const GCType& input, const char* label) {
	GCType gc = input;
	auto d_total = gc.d_contraction();
	auto d_even = gc.d_even_contraction();
	auto d_odd = gc.d_odd_contraction();
	auto sum = d_even;
	sum += d_odd;
	sum.standardize_all();
	d_total.standardize_all();

	const bool ok = (sum.data() == d_total.data());
	std::cout << label
		  << ": total=" << d_total.size()
		  << " even=" << d_even.size()
		  << " odd=" << d_odd.size()
		  << " -> " << (ok ? "ok" : "failed") << '\n';
	return ok;
}

static OddGraphdegZero<10> w9_term_1() {
	OddGraphdegZero<10> g;
	int e = 0;
	for (auto [u, v] : std::initializer_list<std::pair<int, int>>{
		{0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,5}, {2,4}, {2,6},
		{3,5}, {3,7}, {4,6}, {4,8}, {5,7}, {6,9}, {7,8}, {7,9}, {8,9}
	}) {
		g.setEdge(e++, u, v);
	}
	return g;
}

static OddGraphdegZero<10> w9_term_2() {
	OddGraphdegZero<10> g;
	int e = 0;
	for (auto [u, v] : std::initializer_list<std::pair<int, int>>{
		{0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,5}, {2,4}, {2,6},
		{3,5}, {4,7}, {5,6}, {5,8}, {6,8}, {6,9}, {7,8}, {7,9}, {8,9}
	}) {
		g.setEdge(e++, u, v);
	}
	return g;
}

static OddGraphdegZero<10> w9_term_3() {
	OddGraphdegZero<10> g;
	int e = 0;
	for (auto [u, v] : std::initializer_list<std::pair<int, int>>{
		{0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {1,3}, {1,5}, {2,4}, {2,6},
		{3,7}, {3,8}, {4,6}, {5,6}, {5,7}, {5,9}, {6,9}, {7,8}, {8,9}
	}) {
		g.setEdge(e++, u, v);
	}
	return g;
}

int main() {
	bool ok = true;
	OddGCdegZero<10> combined(w9_term_1());
	combined += OddGCdegZero<10>(w9_term_2());
	combined += OddGCdegZero<10>(w9_term_3());

	ok &= check_odd_even_contraction_split(OddGCdegZero<10>(w9_term_1()), "w9_term_1");
	ok &= check_odd_even_contraction_split(OddGCdegZero<10>(w9_term_2()), "w9_term_2");
	ok &= check_odd_even_contraction_split(OddGCdegZero<10>(w9_term_3()), "w9_term_3");
	ok &= check_odd_even_contraction_split(combined, "w9_sum");
	ok &= check_graph_isomorphism_action<OddLoopGraphType<3>>("triangle_iso");
	ok &= check_graph_isomorphism_composition<OddLoopGraphType<3>>("triangle_iso_comp");
	ok &= check_graph_isomorphism_inverse<OddLoopGraphType<3>>("triangle_iso_inv");
	ok &= check_graph_isomorphism_graph_sign<OddLoopGraphType<3>>("triangle_iso_graph_sign");
	ok &= check_generated_isomorphism_standardization_matches_standardize("triangle_iso_std_compare", loop_graph<3>());
	ok &= check_graph_directions_shape<OddLoopGraphType<3>>("triangle_directions");
	ok &= check_graph_directions_isomorphism_action<OddLoopGraphType<3>>("triangle_directions_iso");
	ok &= check_graph_directions_order<OddLoopGraphType<3>>("triangle_directions_order");
	ok &= check_graph_directions_standardize<OddLoopGraphType<3>>("triangle_directions_std");
	ok &= check_ocgraph_standardize_and_sort<OddLoopGraphType<3>>("triangle_ocgraph_std");
	ok &= check_ocgraph_filters_non_covering_directions<OddLoopGraphType<3>>("triangle_ocgraph_filter");
	ok &= check_ocgraph_OCG_F0_same_degree("triangle_ocgraph_OCG_F0_same");
	ok &= check_wheel_ocgraph_delta_e_isomorphism_equivariance<5>("wheel5_ocgraph_delta_e_iso");
	ok &= check_wheel_ocgraph_delta_e_isomorphism_equivariance<7>("wheel7_ocgraph_delta_e_iso");
	ok &= check_full_unsorted_splits_double<7>("wheel7_full_unsorted_splits");
	ok &= check_wheel_ocgc_delta_v_squares_to_zero<7, 0, 0>("wheel7_ocgc_delta_v_squared_c00");
	ok &= check_wheel_ocgc_delta_v_squares_to_zero<7, 0, 1>("wheel7_ocgc_delta_v_squared_c01");
	ok &= check_wheel_ocgc_delta_v_squares_to_zero<7, 1, 0>("wheel7_ocgc_delta_v_squared_c10");
	ok &= check_wheel_ocgc_delta_v_squares_to_zero<7, 1, 1>("wheel7_ocgc_delta_v_squared_c11");
	std::cout << "wheel7_ocgc_delta_anticommute_c00: skipped\n";
	ok &= check_wheel_ocgc_delta_anticommute<7, 0, 1>("wheel7_ocgc_delta_anticommute_c01");
	ok &= check_wheel_ocgc_delta_anticommute<7, 1, 0>("wheel7_ocgc_delta_anticommute_c10");
	ok &= check_wheel_ocgc_delta_anticommute<7, 1, 1>("wheel7_ocgc_delta_anticommute_c11");
	ok &= check_wheel_OCG_F0_exact_equivariance<5>("wheel5_OCG_F0_exact");
	ok &= check_wheel_OCG_F0_exact_equivariance<7>("wheel7_OCG_F0_exact");
	ok &= check_wheel_OCG_F0_exact_equivariance<9>("wheel9_OCG_F0_exact");
	ok &= check_w9_OCG_F0_unique_direction_data("wheel9_OCG_F0_direction_data");
	ok &= check_wheel_OCG_F0_delta_homotopy_equivalence<5>("wheel5_OCG_F0_delta_homotopy");
	ok &= OCGraph_test<OddLoopGraphType<3>>::object();
	ok &= check_minimizing_isomorphisms("triangle_minimizers", loop_graph<3>(), 6);
	ok &= check_minimizing_isomorphisms_match_standardization("wheel_11_minimizers", wheel_graph<11>(), 22);
	ok &= check_generated_isomorphism_standardization_matches_standardize("wheel_11_iso_std_compare", wheel_graph<11>());
	ok &= check_scaled_transpose_matches_contraction_rounds4<7>("wheel7_scaled_transpose_matches_contraction_rounds4");
	ok &= check_v_graph_sort_matches_quick<25>("V25_sort_matches_quick");
	ok &= check_v_graph_permuted_sort_matches_old_path<25>("V25_permuted_sort_matches_old_path");

	if (!ok) {
		return EXIT_FAILURE;
	}

	std::cout << "all tests passed\n";
	return EXIT_SUCCESS;
}
