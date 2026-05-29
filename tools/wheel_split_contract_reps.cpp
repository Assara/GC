#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <optional>
#include <string>
#include <unordered_set>
#include <vector>

#include "examplegraphs.hpp"
#include "VectorSpace/row_based_lil_matrix.hpp"

namespace {

enum class ExperimentMode {
	GeneratedSupport,
	GeneratedConstrained,
	Krylov,
	KrylovDual,
	KrylovConstrained,
	KrylovConstrainedDual
};

template <typename Fn>
double time_seconds(Fn&& fn) {
	const auto start = std::chrono::steady_clock::now();
	fn();
	const auto stop = std::chrono::steady_clock::now();
	return std::chrono::duration<double>(stop - start).count();
}

template <Int N>
std::string cache_path(int rounds) {
	return (std::filesystem::path("output")
		/ "split_contract"
		/ "cache"
		/ ("gc_wheel_same_degree_candidates_W" + std::to_string(N)
			+ "_rounds" + std::to_string(rounds)
			+ "_int" + std::to_string(sizeof(Int))
			+ "_v1.bin")).string();
}

template <typename GraphType>
bool load_graph_cache(const std::string& path, std::vector<GraphType>& graphs) {
	std::ifstream in(path, std::ios::binary);
	if (!in) {
		return false;
	}

	std::uint64_t count = 0;
	in.read(reinterpret_cast<char*>(&count), sizeof(count));
	if (!in) {
		return false;
	}

	graphs.resize(static_cast<std::size_t>(count));
	for (auto& graph : graphs) {
		in.read(
			reinterpret_cast<char*>(graph.half_edges.data()),
			static_cast<std::streamsize>(graph.half_edges.size() * sizeof(Int))
		);
		if (!in) {
			graphs.clear();
			return false;
		}
	}

	return true;
}

template <typename GraphType>
void save_graph_cache(const std::string& path, const std::vector<GraphType>& graphs) {
	const std::filesystem::path cache_file(path);
	if (cache_file.has_parent_path()) {
		std::filesystem::create_directories(cache_file.parent_path());
	}

	std::ofstream out(path, std::ios::binary);
	if (!out) {
		std::cerr << "warning: failed to write cache " << path << '\n';
		return;
	}

	const std::uint64_t count = graphs.size();
	out.write(reinterpret_cast<const char*>(&count), sizeof(count));
	for (const auto& graph : graphs) {
		out.write(
			reinterpret_cast<const char*>(graph.half_edges.data()),
			static_cast<std::streamsize>(graph.half_edges.size() * sizeof(Int))
		);
	}
}

template <typename GraphType>
bool write_representatives(
	const std::filesystem::path& output_path,
	const std::vector<GraphType>& graphs
) {
	if (output_path.has_parent_path()) {
		std::filesystem::create_directories(output_path.parent_path());
	}

	std::ofstream out(output_path, std::ios::trunc);
	if (!out) {
		return false;
	}

	out << "graph_size: (" << GraphType::N_VERTICES_ << "," << GraphType::N_EDGES_ << ")\n";
	out << "field_type: " << TypeName<fieldType>::name() << "\n";
	out << "number_of_graphs: " << graphs.size() << "\n";
	for (const auto& graph : graphs) {
		for (Int edge = 0; edge < GraphType::N_EDGES_; ++edge) {
			const auto [u, v] = graph.getEdge(edge);
			out << "(" << +u << "," << +v << ")";
			if (edge + 1 < GraphType::N_EDGES_) {
				out << ", ";
			}
		}
		out << "\n";
	}

	return true;
}

std::filesystem::path round_snapshot_path(
	const std::filesystem::path* output_path,
	Int wheel_n,
	int round
) {
	if (output_path != nullptr) {
		const std::filesystem::path parent = output_path->has_parent_path()
			? output_path->parent_path()
			: std::filesystem::path{};
		const std::string stem = output_path->stem().string();
		return parent / (stem + "_round" + std::to_string(round) + ".txt");
	}

	return std::filesystem::path("output")
		/ "split_contract"
		/ ("gc_wheel_same_degree_candidates_W" + std::to_string(wheel_n)
			+ "_round" + std::to_string(round) + ".txt");
}

template <typename GCType>
bool write_class_file(
	const std::filesystem::path& output_path,
	const GCType& gamma
) {
	if (output_path.has_parent_path()) {
		std::filesystem::create_directories(output_path.parent_path());
	}

	std::ofstream out(output_path, std::ios::trunc);
	if (!out) {
		return false;
	}

	using GraphType = typename GCType::GraphType;
	out << "graph_size: (" << GraphType::N_VERTICES_ << "," << GraphType::N_EDGES_ << ")\n";
	out << "field_type: " << TypeName<fieldType>::name() << "\n";
	out << "number_of_graphs: " << gamma.data().size() << "\n";
	for (const auto& be : gamma.data()) {
		out << be.getCoefficient() << "; ";
		for (Int edge = 0; edge < GraphType::N_EDGES_; ++edge) {
			const auto [u, v] = be.getValue().getEdge(edge);
			out << "(" << +u << "," << +v << ")";
			if (edge + 1 < GraphType::N_EDGES_) {
				out << ", ";
			}
		}
		out << "\n";
	}

	return true;
}

template <typename GCType>
std::filesystem::path related_output_path(
	const std::filesystem::path* output_path,
	const std::string& suffix
) {
	if (output_path != nullptr) {
		const std::filesystem::path parent = output_path->has_parent_path()
			? output_path->parent_path()
			: std::filesystem::path{};
		return parent / (output_path->stem().string() + suffix + output_path->extension().string());
	}

	using GraphType = typename GCType::GraphType;
	return std::filesystem::path("output")
		/ "split_contract"
		/ "debug"
		/ ("gc_debug_"
			+ std::to_string(GraphType::N_VERTICES_)
			+ "_"
			+ std::to_string(GraphType::N_EDGES_)
			+ suffix
			+ ".txt");
}

template <typename GraphType>
bool canonicalize_nonzero(typename GraphType::Basis& basis) {
	if (basis.getCoefficient() == fieldType{}) {
		return false;
	}
	GraphType::std(basis);
	return basis.getCoefficient() != fieldType{};
}

template <Int N>
std::vector<OddGraphdegZero<N + 1>> generate_same_degree_candidates(
	int rounds,
	const std::filesystem::path* snapshot_output_path
) {
	using GraphType = OddGraphdegZero<N + 1>;

	std::unordered_set<GraphType> seen;
	std::vector<GraphType> all_graphs;
	std::vector<GraphType> frontier;

	typename GraphType::Basis wheel(wheel_graph<N>(), fieldType{1});
	GraphType::std(wheel);
	seen.insert(wheel.getValue());
	all_graphs.push_back(wheel.getValue());
	frontier.push_back(wheel.getValue());

	std::cout << "round 0: frontier = " << frontier.size()
	          << ", cumulative = " << all_graphs.size() << '\n';
	const auto round0_path = round_snapshot_path(snapshot_output_path, N, 0);
	if (write_representatives(round0_path, all_graphs)) {
		std::cout << "round 0 written = " << round0_path.string() << '\n';
	}

	for (int round = 1; round <= rounds; ++round) {
		std::unordered_set<GraphType> next_seen;
		std::size_t split_terms = 0;
		std::size_t raw_splits = 0;
		std::size_t raw_contracts = 0;
		std::size_t nonzero_contracts = 0;
		std::size_t new_candidates_before_round_dedupe = 0;

		const double seconds = time_seconds([&]() {
			for (const auto& graph : frontier) {
				const auto splits = graph.unsorted_splits(fieldType{1});
				split_terms += static_cast<std::size_t>(splits.raw_elements().size());
				for (const auto& split_be : splits.raw_elements()) {
					if (split_be.getCoefficient() == fieldType{}) {
						continue;
					}
					++raw_splits;

					for (Int edge = 0; edge < GraphType::SplitGraph::N_EDGES_; ++edge) {
						auto contracted = split_be.getValue().contract_edge(edge, split_be.getCoefficient());
						++raw_contracts;
						if (!canonicalize_nonzero<typename GraphType::SplitGraph::ContGraph>(contracted)) {
							continue;
						}
						++nonzero_contracts;

						const auto& candidate = contracted.getValue();
						if (!seen.contains(candidate)) {
							++new_candidates_before_round_dedupe;
							next_seen.insert(candidate);
						}
					}
				}
			}
		});

		std::cout << "round " << round
		          << " split step: input frontier = " << frontier.size()
		          << ", split terms = " << split_terms
		          << ", nonzero splits = " << raw_splits << '\n';
		std::cout << "round " << round
		          << " contraction step: raw contracts = " << raw_contracts
		          << ", nonzero contracts = " << nonzero_contracts
		          << ", new candidates before round dedupe = "
		          << new_candidates_before_round_dedupe
		          << ", new candidates after round dedupe = " << next_seen.size()
		          << '\n';

		frontier.assign(next_seen.begin(), next_seen.end());
		std::sort(frontier.begin(), frontier.end());
		for (const auto& graph : frontier) {
			seen.insert(graph);
			all_graphs.push_back(graph);
		}
		std::sort(all_graphs.begin(), all_graphs.end());

		std::cout << "round " << round
		          << ": raw_splits = " << raw_splits
		          << ", raw_contracts = " << raw_contracts
		          << ", nonzero_contracts = " << nonzero_contracts
		          << ", new frontier = " << frontier.size()
		          << ", cumulative = " << all_graphs.size()
		          << ", time = " << seconds << " s\n";

		const auto round_path = round_snapshot_path(snapshot_output_path, N, round);
		if (write_representatives(round_path, all_graphs)) {
			std::cout << "round " << round << " written = "
			          << round_path.string() << '\n';
		}

		if (frontier.empty()) {
			break;
		}
	}

	return all_graphs;
}

template <Int N>
std::optional<OddGCdegZero<N + 1>> solve_same_degree_representative(
	const std::vector<OddGraphdegZero<N + 1>>& support,
	const std::filesystem::path* representative_output_path
) {
	using GCType = OddGCdegZero<N + 1>;
	using GraphType = typename GCType::GraphType;
	using SplitGraph = typename GCType::SplitGraphType;
	using SplitL = typename GCType::SplitL;

	typename GraphType::Basis wheel_basis(wheel_graph<N>(), fieldType{1});
	GraphType::std(wheel_basis);
	const GraphType wheel = wheel_basis.getValue();

	SplitL target = GCType(wheel, AssumeBasisOrderTag{}).delta().data();

	std::unordered_map<GraphType, SplitL> split_map;
	split_map.reserve(support.size());
	std::size_t skipped_wheel = 0;

	const double map_seconds = time_seconds([&]() {
		for (const GraphType& graph : support) {
			if (graph == wheel) {
				++skipped_wheel;
				continue;
			}
			auto column = GCType(graph, AssumeBasisOrderTag{}).delta().data();
			if (!column.empty()) {
				split_map.emplace(graph, std::move(column));
			}
		}
	});

	std::cout << "split map domain = " << split_map.size() << '\n';
	std::cout << "skipped wheel columns = " << skipped_wheel << '\n';
	std::cout << "target delta(W) terms = " << target.size() << '\n';
	std::cout << "map build time = " << map_seconds << " s\n";

	VectorSpace::wiedemann_primitive_finder<SplitGraph, GraphType, fieldType> solver(split_map);

	std::optional<typename GCType::L> correction;
	const double solve_seconds = time_seconds([&]() {
		correction = solver.find_primitive_or_empty(target);
	});
	std::cout << "solve time = " << solve_seconds << " s\n";

	if (!correction.has_value()) {
		return std::nullopt;
	}

	auto correction_copy = *correction;
	GCType correction_gc(correction_copy);
	auto delta_correction = correction_gc.delta();
	SplitL direct_residual = target.add_scaled(delta_correction.data(), fieldType{-1});
	typename GCType::SplitGC direct_residual_gc(direct_residual);

	std::vector<typename GCType::Base> terms;
	terms.reserve(static_cast<std::size_t>(correction->size()) + 1);
	terms.emplace_back(wheel, fieldType{1});
	for (const auto& be : *correction) {
		terms.emplace_back(be.getValue(), -be.getCoefficient());
	}

	GCType representative(std::move(terms));
	auto delta_representative = representative.delta();
	std::cout << "correction terms = " << correction->size() << '\n';
	std::cout << "delta(correction) terms = " << delta_correction.data().size() << '\n';
	std::cout << "delta(correction) == delta(W): "
	          << ((delta_correction.data() == target) ? "yes" : "no") << '\n';
	std::cout << "direct delta(W)-delta(correction) terms = "
	          << direct_residual.size() << '\n';
	std::cout << "representative terms = " << representative.data().size() << '\n';
	std::cout << "delta(representative) terms = " << delta_representative.data().size() << '\n';

	if (!delta_representative.data().empty()) {
		const auto target_path =
			related_output_path<typename GCType::SplitGC>(representative_output_path, "_delta_wheel");
		const auto correction_path =
			related_output_path<typename GCType::SplitGC>(representative_output_path, "_delta_correction");
		const auto direct_residual_path =
			related_output_path<typename GCType::SplitGC>(representative_output_path, "_direct_delta_residual");
		const auto representative_residual_path =
			related_output_path<typename GCType::SplitGC>(representative_output_path, "_delta_residual");

		typename GCType::SplitGC target_gc(target);
		write_class_file(target_path, target_gc);
		write_class_file(correction_path, delta_correction);
		write_class_file(direct_residual_path, direct_residual_gc);
		write_class_file(representative_residual_path, delta_representative);

		std::cout << "saved delta(W) = " << target_path.string() << '\n';
		std::cout << "saved delta(correction) = " << correction_path.string() << '\n';
		std::cout << "saved direct residual = " << direct_residual_path.string() << '\n';
		std::cout << "saved representative residual = "
		          << representative_residual_path.string() << '\n';
		return std::nullopt;
	}

	return representative;
}

template <Int N>
int run_case(
	int rounds,
	const std::filesystem::path* representatives_output_path,
	const std::filesystem::path* representative_output_path
) {
	using GraphType = OddGraphdegZero<N + 1>;

	std::vector<GraphType> graphs;
	const std::string path = cache_path<N>(rounds);
	if (load_graph_cache(path, graphs)) {
		std::cout << "loaded cache = " << path << '\n';
	} else {
		graphs = generate_same_degree_candidates<N>(rounds, representatives_output_path);
		save_graph_cache(path, graphs);
		std::cout << "saved cache = " << path << '\n';
	}

	std::cout << "wheel = W" << +N << '\n';
	std::cout << "split-contract rounds = " << rounds << '\n';
	std::cout << "same-degree candidates = " << graphs.size() << '\n';

	if (representatives_output_path != nullptr) {
		if (!write_representatives(*representatives_output_path, graphs)) {
			std::cerr << "failed to write " << representatives_output_path->string() << '\n';
			return EXIT_FAILURE;
		}
		std::cout << "saved representatives = " << representatives_output_path->string() << '\n';
	}

	auto representative = solve_same_degree_representative<N>(graphs, representative_output_path);
	if (!representative.has_value()) {
		std::cout << "representative found = no\n";
		return EXIT_FAILURE;
	}

	std::cout << "representative found = yes\n";
	if (representative_output_path != nullptr) {
		if (!write_class_file(*representative_output_path, *representative)) {
			std::cerr << "failed to write " << representative_output_path->string() << '\n';
			return EXIT_FAILURE;
		}
		std::cout << "saved representative = "
		          << representative_output_path->string() << '\n';
	}

	return EXIT_SUCCESS;
}

template <typename GCType>
GCType split_contract_krylov_step(const GCType& gamma) {
	GCType next = gamma.delta().d_contraction();
	next.standardize_all();
	next.sort_elements();
	return next;
}

template <Int N>
std::optional<OddGCdegZero<N + 1>> solve_krylov_representative(
	int rounds,
	const std::filesystem::path* representative_output_path
) {
	using GCType = OddGCdegZero<N + 1>;
	using GraphType = typename GCType::GraphType;
	using SplitGraph = typename GCType::SplitGraphType;
	using SplitL = typename GCType::SplitL;

	typename GraphType::Basis wheel_basis(wheel_graph<N>(), fieldType{1});
	GraphType::std(wheel_basis);
	const GraphType wheel = wheel_basis.getValue();

	std::vector<GCType> krylov_vectors;
	krylov_vectors.reserve(static_cast<std::size_t>(rounds) + 1);
	krylov_vectors.emplace_back(wheel, AssumeBasisOrderTag{});

	std::cout << "krylov step 0 terms = "
	          << krylov_vectors.back().data().size() << '\n';
	for (int step = 1; step <= rounds; ++step) {
		const double seconds = time_seconds([&]() {
			krylov_vectors.push_back(split_contract_krylov_step(krylov_vectors.back()));
		});
		std::cout << "krylov step " << step
		          << " terms = " << krylov_vectors.back().data().size()
		          << ", time = " << seconds << " s\n";
	}

	SplitL target = GCType(wheel, AssumeBasisOrderTag{}).delta().data();
	std::unordered_map<std::size_t, SplitL> split_map;
	split_map.reserve(static_cast<std::size_t>(rounds));

	const double map_seconds = time_seconds([&]() {
		for (std::size_t index = 1; index < krylov_vectors.size(); ++index) {
			SplitL column = krylov_vectors[index].delta().data();
			std::cout << "delta(krylov step " << index
			          << ") terms = " << column.size() << '\n';
			if (!column.empty()) {
				split_map.emplace(index, std::move(column));
			}
		}
	});

	std::cout << "krylov map domain = " << split_map.size() << '\n';
	std::cout << "target delta(W) terms = " << target.size() << '\n';
	std::cout << "map build time = " << map_seconds << " s\n";

	VectorSpace::wiedemann_primitive_finder<SplitGraph, std::size_t, fieldType> solver(split_map);

	std::optional<VectorSpace::LinComb<std::size_t, fieldType>> correction_coefficients;
	const double solve_seconds = time_seconds([&]() {
		correction_coefficients = solver.find_primitive_or_empty(target);
	});
	std::cout << "solve time = " << solve_seconds << " s\n";

	if (!correction_coefficients.has_value()) {
		return std::nullopt;
	}

	GCType correction;
	for (const auto& be : *correction_coefficients) {
		const std::size_t index = be.getValue();
		if (index >= krylov_vectors.size()) {
			std::cerr << "internal error: Krylov coefficient index out of range\n";
			return std::nullopt;
		}
		GCType scaled = krylov_vectors[index];
		scaled.scalar_multiply(be.getCoefficient());
		correction += scaled;
	}
	correction.standardize_all();
	correction.sort_elements();

	GCType representative(wheel, AssumeBasisOrderTag{});
	GCType negative_correction = correction;
	negative_correction.scalar_multiply(fieldType{-1});
	representative += negative_correction;
	representative.standardize_all();
	representative.sort_elements();

	const auto delta_correction = correction.delta();
	const auto delta_representative = representative.delta();
	std::cout << "correction Krylov coefficients = "
	          << correction_coefficients->size() << '\n';
	std::cout << "correction graph terms = " << correction.data().size() << '\n';
	std::cout << "delta(correction) terms = "
	          << delta_correction.data().size() << '\n';
	std::cout << "delta(correction) == delta(W): "
	          << ((delta_correction.data() == target) ? "yes" : "no") << '\n';
	std::cout << "representative terms = "
	          << representative.data().size() << '\n';
	std::cout << "delta(representative) terms = "
	          << delta_representative.data().size() << '\n';

	if (!delta_representative.data().empty()) {
		const auto representative_residual_path =
			related_output_path<typename GCType::SplitGC>(representative_output_path, "_krylov_delta_residual");
		write_class_file(representative_residual_path, delta_representative);
		std::cout << "saved representative residual = "
		          << representative_residual_path.string() << '\n';
		return std::nullopt;
	}

	return representative;
}

template <Int N>
int run_krylov_case(
	int rounds,
	const std::filesystem::path* representative_output_path
) {
	std::cout << "wheel = W" << +N << '\n';
	std::cout << "krylov rounds = " << rounds << '\n';

	auto representative = solve_krylov_representative<N>(rounds, representative_output_path);
	if (!representative.has_value()) {
		std::cout << "krylov representative found = no\n";
		return EXIT_FAILURE;
	}

	std::cout << "krylov representative found = yes\n";
	if (representative_output_path != nullptr) {
		if (!write_class_file(*representative_output_path, *representative)) {
			std::cerr << "failed to write " << representative_output_path->string() << '\n';
			return EXIT_FAILURE;
		}
		std::cout << "saved representative = "
		          << representative_output_path->string() << '\n';
	}

	return EXIT_SUCCESS;
}

template <Int N>
int run_krylov_dual_case(
	int rounds,
	const std::filesystem::path* witness_output_path
) {
	using GCType = OddGCdegZero<N + 1>;
	using GraphType = typename GCType::GraphType;
	using SplitGraph = typename GCType::SplitGraphType;
	using SplitL = typename GCType::SplitL;

	typename GraphType::Basis wheel_basis(wheel_graph<N>(), fieldType{1});
	GraphType::std(wheel_basis);
	const GraphType wheel = wheel_basis.getValue();

	std::vector<GCType> krylov_vectors;
	krylov_vectors.reserve(static_cast<std::size_t>(rounds) + 1);
	krylov_vectors.emplace_back(wheel, AssumeBasisOrderTag{});

	std::cout << "wheel = W" << +N << '\n';
	std::cout << "krylov dual depth = " << rounds << '\n';
	std::cout << "krylov step 0 terms = "
	          << krylov_vectors.back().data().size() << '\n';
	for (int step = 1; step <= rounds; ++step) {
		const double seconds = time_seconds([&]() {
			krylov_vectors.push_back(split_contract_krylov_step(krylov_vectors.back()));
		});
		std::cout << "krylov step " << step
		          << " terms = " << krylov_vectors.back().data().size()
		          << ", time = " << seconds << " s\n";
	}

	SplitL target = GCType(wheel, AssumeBasisOrderTag{}).delta().data();
	std::vector<SplitL> columns;
	columns.reserve(static_cast<std::size_t>(rounds));
	for (std::size_t index = 1; index < krylov_vectors.size(); ++index) {
		columns.push_back(krylov_vectors[index].delta().data());
		std::cout << "delta(krylov step " << index
		          << ") terms = " << columns.back().size() << '\n';
	}

	std::unordered_map<SplitGraph, std::size_t> image_index;
	std::vector<SplitGraph> image_basis;
	auto get_image_index = [&](const SplitGraph& graph) {
		auto [it, inserted] = image_index.try_emplace(graph, image_basis.size());
		if (inserted) {
			image_basis.push_back(graph);
		}
		return it->second;
	};

	row_based_lil_matrix<fieldType> equations(columns.size() + 1);
	for (std::size_t column_index = 0; column_index < columns.size(); ++column_index) {
		for (const auto& be : columns[column_index]) {
			equations.add_element(
				column_index,
				get_image_index(be.getValue()),
				be.getCoefficient()
			);
		}
	}

	const std::size_t target_row = columns.size();
	for (const auto& be : target) {
		equations.add_element(
			target_row,
			get_image_index(be.getValue()),
			be.getCoefficient()
		);
	}
	equations.sort_rows();

	VectorSpace::LinComb<std::size_t, fieldType> rhs;
	rhs.append_in_basis_order(target_row, fieldType{1});

	std::optional<VectorSpace::LinComb<std::size_t, fieldType>> witness_indices;
	const double solve_seconds = time_seconds([&]() {
		witness_indices = equations.solve_MX_equals_y(rhs);
	});

	std::cout << "dual equations = " << (columns.size() + 1) << '\n';
	std::cout << "dual image basis = " << image_basis.size() << '\n';
	std::cout << "dual matrix entries = " << equations.size() << '\n';
	std::cout << "dual solve time = " << solve_seconds << " s\n";

	if (!witness_indices.has_value()) {
		std::cout << "dual witness found = no\n";
		return EXIT_FAILURE;
	}

	SplitL witness;
	for (const auto& be : *witness_indices) {
		witness.append_in_basis_order(
			image_basis[be.getValue()],
			be.getCoefficient()
		);
	}
	witness.sort_without_deduplicate();

	auto evaluate = [](const SplitL& functional, const SplitL& vector) {
		fieldType result{};
		std::size_t i = 0;
		std::size_t j = 0;
		while (i < functional.raw_elements().size() && j < vector.raw_elements().size()) {
			const auto& lhs = functional.raw_elements()[i];
			const auto& rhs_term = vector.raw_elements()[j];
			if (lhs.getValue() < rhs_term.getValue()) {
				++i;
			} else if (rhs_term.getValue() < lhs.getValue()) {
				++j;
			} else {
				result += lhs.getCoefficient() * rhs_term.getCoefficient();
				++i;
				++j;
			}
		}
		return result;
	};

	std::size_t nonzero_column_pairings = 0;
	for (const auto& column : columns) {
		if (evaluate(witness, column) != fieldType{}) {
			++nonzero_column_pairings;
		}
	}
	const fieldType target_pairing = evaluate(witness, target);

	std::cout << "dual witness terms = " << witness.size() << '\n';
	std::cout << "nonzero pairings with Krylov columns = "
	          << nonzero_column_pairings << '\n';
	std::cout << "pairing with delta(W) = " << target_pairing << '\n';
	std::cout << "dual witness found = "
	          << ((nonzero_column_pairings == 0 && target_pairing != fieldType{})
	              ? "yes"
	              : "no") << '\n';

	if (nonzero_column_pairings != 0 || target_pairing == fieldType{}) {
		return EXIT_FAILURE;
	}

	if (witness_output_path != nullptr) {
		typename GCType::SplitGC witness_gc(witness);
		if (!write_class_file(*witness_output_path, witness_gc)) {
			std::cerr << "failed to write " << witness_output_path->string() << '\n';
			return EXIT_FAILURE;
		}
		std::cout << "saved dual witness = "
		          << witness_output_path->string() << '\n';
	}

	return EXIT_SUCCESS;
}

template <typename GCType>
class CombinedSplitContractIndex {
public:
	using GraphType = typename GCType::GraphType;
	using SplitGraph = typename GCType::SplitGraphType;
	using ContGraph = typename GCType::ContGraphType;

	struct BasisRef {
		bool is_split;
		std::size_t local_index;
	};

	std::size_t split_index(const SplitGraph& graph) {
		auto [it, inserted] = split_indices_.try_emplace(graph, refs_.size());
		if (inserted) {
			refs_.push_back(BasisRef{true, split_basis_.size()});
			split_basis_.push_back(graph);
		}
		return it->second;
	}

	std::size_t contraction_index(const ContGraph& graph) {
		auto [it, inserted] = contraction_indices_.try_emplace(graph, refs_.size());
		if (inserted) {
			refs_.push_back(BasisRef{false, contraction_basis_.size()});
			contraction_basis_.push_back(graph);
		}
		return it->second;
	}

	std::size_t size() const {
		return refs_.size();
	}

	const std::vector<BasisRef>& refs() const {
		return refs_;
	}

	const std::vector<SplitGraph>& split_basis() const {
		return split_basis_;
	}

	const std::vector<ContGraph>& contraction_basis() const {
		return contraction_basis_;
	}

private:
	std::unordered_map<SplitGraph, std::size_t> split_indices_;
	std::unordered_map<ContGraph, std::size_t> contraction_indices_;
	std::vector<BasisRef> refs_;
	std::vector<SplitGraph> split_basis_;
	std::vector<ContGraph> contraction_basis_;
};

template <typename GCType>
typename GCType::ContGC contraction_of(GCType gamma) {
	auto result = gamma.d_contraction();
	result.standardize_all();
	result.sort_elements();
	return result;
}

template <Int N>
std::optional<OddGCdegZero<N + 1>> solve_constrained_krylov_representative(
	int rounds,
	const std::filesystem::path* representative_output_path
) {
	using GCType = OddGCdegZero<N + 1>;
	using GraphType = typename GCType::GraphType;

	typename GraphType::Basis wheel_basis(wheel_graph<N>(), fieldType{1});
	GraphType::std(wheel_basis);
	const GraphType wheel = wheel_basis.getValue();

	std::vector<GCType> krylov_vectors;
	krylov_vectors.reserve(static_cast<std::size_t>(rounds) + 1);
	krylov_vectors.emplace_back(wheel, AssumeBasisOrderTag{});

	std::cout << "wheel = W" << +N << '\n';
	std::cout << "constrained krylov depth = " << rounds << '\n';
	std::cout << "krylov step 0 terms = "
	          << krylov_vectors.back().data().size() << '\n';
	for (int step = 1; step <= rounds; ++step) {
		const double seconds = time_seconds([&]() {
			krylov_vectors.push_back(split_contract_krylov_step(krylov_vectors.back()));
		});
		std::cout << "krylov step " << step
		          << " terms = " << krylov_vectors.back().data().size()
		          << ", time = " << seconds << " s\n";
	}

	CombinedSplitContractIndex<GCType> image_index;
	row_based_lil_matrix<fieldType> matrix;
	for (std::size_t column_index = 1; column_index < krylov_vectors.size(); ++column_index) {
		const auto split_column = krylov_vectors[column_index].delta().data();
		const auto contraction_column = contraction_of(krylov_vectors[column_index]).data();
		std::cout << "column " << column_index
		          << ": d_split terms = " << split_column.size()
		          << ", d_contract terms = " << contraction_column.size() << '\n';

		for (const auto& be : split_column) {
			matrix.add_element(
				image_index.split_index(be.getValue()),
				column_index - 1,
				be.getCoefficient()
			);
		}
		for (const auto& be : contraction_column) {
			matrix.add_element(
				image_index.contraction_index(be.getValue()),
				column_index - 1,
				be.getCoefficient()
			);
		}
	}
	matrix.sort_rows();

	VectorSpace::LinComb<std::size_t, fieldType> target;
	const auto split_target = GCType(wheel, AssumeBasisOrderTag{}).delta().data();
	for (const auto& be : split_target) {
		target.append_in_basis_order(
			image_index.split_index(be.getValue()),
			be.getCoefficient()
		);
	}
	target.sort_without_deduplicate();

	std::optional<VectorSpace::LinComb<std::size_t, fieldType>> correction_coefficients;
	const double solve_seconds = time_seconds([&]() {
		correction_coefficients = matrix.solve_MX_equals_y(target);
	});

	std::cout << "combined image basis = " << image_index.size() << '\n';
	std::cout << "split image basis = " << image_index.split_basis().size() << '\n';
	std::cout << "contraction image basis = "
	          << image_index.contraction_basis().size() << '\n';
	std::cout << "matrix entries = " << matrix.size() << '\n';
	std::cout << "target d_split(W) terms = " << split_target.size() << '\n';
	std::cout << "target d_contract terms = 0\n";
	std::cout << "solve time = " << solve_seconds << " s\n";

	if (!correction_coefficients.has_value()) {
		return std::nullopt;
	}

	GCType correction;
	for (const auto& be : *correction_coefficients) {
		const std::size_t index = be.getValue() + 1;
		GCType scaled = krylov_vectors[index];
		scaled.scalar_multiply(be.getCoefficient());
		correction += scaled;
	}
	correction.standardize_all();
	correction.sort_elements();

	GCType representative(wheel, AssumeBasisOrderTag{});
	GCType negative_correction = correction;
	negative_correction.scalar_multiply(fieldType{-1});
	representative += negative_correction;
	representative.standardize_all();
	representative.sort_elements();

	const auto correction_split = correction.delta();
	const auto correction_contract = contraction_of(correction);
	const auto representative_split = representative.delta();
	std::cout << "correction Krylov coefficients = "
	          << correction_coefficients->size() << '\n';
	std::cout << "correction graph terms = " << correction.data().size() << '\n';
	std::cout << "d_split(correction) terms = "
	          << correction_split.data().size() << '\n';
	std::cout << "d_contract(correction) terms = "
	          << correction_contract.data().size() << '\n';
	std::cout << "d_split(correction) == d_split(W): "
	          << ((correction_split.data() == split_target) ? "yes" : "no") << '\n';
	std::cout << "d_split(representative) terms = "
	          << representative_split.data().size() << '\n';

	if (!correction_contract.data().empty() || !representative_split.data().empty()) {
		if (representative_output_path != nullptr) {
			write_class_file(
				related_output_path<typename GCType::ContGC>(representative_output_path, "_correction_contract"),
				correction_contract
			);
			write_class_file(
				related_output_path<typename GCType::SplitGC>(representative_output_path, "_representative_split"),
				representative_split
			);
		}
		return std::nullopt;
	}

	return representative;
}

template <Int N>
int run_constrained_krylov_case(
	int rounds,
	const std::filesystem::path* representative_output_path
) {
	auto representative = solve_constrained_krylov_representative<N>(
		rounds,
		representative_output_path
	);
	if (!representative.has_value()) {
		std::cout << "constrained krylov representative found = no\n";
		return EXIT_FAILURE;
	}

	std::cout << "constrained krylov representative found = yes\n";
	if (representative_output_path != nullptr) {
		if (!write_class_file(*representative_output_path, *representative)) {
			std::cerr << "failed to write " << representative_output_path->string() << '\n';
			return EXIT_FAILURE;
		}
		std::cout << "saved representative = "
		          << representative_output_path->string() << '\n';
	}
	return EXIT_SUCCESS;
}

template <Int N>
std::optional<OddGCdegZero<N + 1>> solve_generated_constrained_representative(
	const std::vector<OddGraphdegZero<N + 1>>& support,
	const std::filesystem::path* representative_output_path
) {
	using GCType = OddGCdegZero<N + 1>;
	using GraphType = typename GCType::GraphType;

	typename GraphType::Basis wheel_basis(wheel_graph<N>(), fieldType{1});
	GraphType::std(wheel_basis);
	const GraphType wheel = wheel_basis.getValue();

	std::vector<GraphType> domain_graphs;
	domain_graphs.reserve(support.size());
	std::size_t skipped_wheel = 0;
	for (const GraphType& graph : support) {
		if (graph == wheel) {
			++skipped_wheel;
			continue;
		}
		domain_graphs.push_back(graph);
	}

	CombinedSplitContractIndex<GCType> image_index;
	row_based_lil_matrix<fieldType> matrix;
	std::size_t split_entries = 0;
	std::size_t contraction_entries = 0;

	const double map_seconds = time_seconds([&]() {
		for (std::size_t column_index = 0; column_index < domain_graphs.size(); ++column_index) {
			GCType gamma(domain_graphs[column_index], AssumeBasisOrderTag{});
			const auto split_column = gamma.delta().data();
			const auto contraction_column = contraction_of(gamma).data();
			split_entries += static_cast<std::size_t>(split_column.size());
			contraction_entries += static_cast<std::size_t>(contraction_column.size());

			for (const auto& be : split_column) {
				matrix.add_element(
					image_index.split_index(be.getValue()),
					column_index,
					be.getCoefficient()
				);
			}
			for (const auto& be : contraction_column) {
				matrix.add_element(
					image_index.contraction_index(be.getValue()),
					column_index,
					be.getCoefficient()
				);
			}
		}
		matrix.sort_rows();
	});

	VectorSpace::LinComb<std::size_t, fieldType> target;
	const auto split_target = GCType(wheel, AssumeBasisOrderTag{}).delta().data();
	const std::size_t image_basis_before_target = image_index.size();
	for (const auto& be : split_target) {
		target.append_in_basis_order(
			image_index.split_index(be.getValue()),
			be.getCoefficient()
		);
	}
	target.sort_without_deduplicate();

	std::cout << "combined map domain = " << domain_graphs.size() << '\n';
	std::cout << "skipped wheel columns = " << skipped_wheel << '\n';
	std::cout << "combined image basis = " << image_index.size() << '\n';
	std::cout << "split image basis = " << image_index.split_basis().size() << '\n';
	std::cout << "contraction image basis = "
	          << image_index.contraction_basis().size() << '\n';
	std::cout << "raw split entries = " << split_entries << '\n';
	std::cout << "raw contraction entries = " << contraction_entries << '\n';
	std::cout << "matrix entries before solve = " << matrix.size() << '\n';
	std::cout << "target d_split(W) terms = " << split_target.size() << '\n';
	std::cout << "target d_contract terms = 0\n";
	std::cout << "map build time = " << map_seconds << " s\n";

	if (image_index.size() != image_basis_before_target) {
		std::cout << "target terms missing from combined image basis = "
		          << (image_index.size() - image_basis_before_target) << '\n';
		return std::nullopt;
	}

	std::optional<VectorSpace::LinComb<std::size_t, fieldType>> correction_coefficients;
	const double solve_seconds = time_seconds([&]() {
		correction_coefficients = matrix.solve_MX_equals_y(target);
	});
	std::cout << "solve time = " << solve_seconds << " s\n";

	if (!correction_coefficients.has_value()) {
		return std::nullopt;
	}

	GCType correction;
	for (const auto& be : *correction_coefficients) {
		GCType scaled(domain_graphs[be.getValue()], AssumeBasisOrderTag{});
		scaled.scalar_multiply(be.getCoefficient());
		correction += scaled;
	}
	correction.standardize_all();
	correction.sort_elements();

	GCType representative(wheel, AssumeBasisOrderTag{});
	GCType negative_correction = correction;
	negative_correction.scalar_multiply(fieldType{-1});
	representative += negative_correction;
	representative.standardize_all();
	representative.sort_elements();

	const auto correction_split = correction.delta();
	const auto correction_contract = contraction_of(correction);
	const auto representative_split = representative.delta();

	std::cout << "correction coefficients = "
	          << correction_coefficients->size() << '\n';
	std::cout << "correction graph terms = " << correction.data().size() << '\n';
	std::cout << "d_split(correction) terms = "
	          << correction_split.data().size() << '\n';
	std::cout << "d_contract(correction) terms = "
	          << correction_contract.data().size() << '\n';
	std::cout << "d_split(correction) == d_split(W): "
	          << ((correction_split.data() == split_target) ? "yes" : "no") << '\n';
	std::cout << "d_split(representative) terms = "
	          << representative_split.data().size() << '\n';

	if (!correction_contract.data().empty() || !representative_split.data().empty()) {
		if (representative_output_path != nullptr) {
			write_class_file(
				related_output_path<typename GCType::ContGC>(representative_output_path, "_correction_contract"),
				correction_contract
			);
			write_class_file(
				related_output_path<typename GCType::SplitGC>(representative_output_path, "_representative_split"),
				representative_split
			);
		}
		return std::nullopt;
	}

	return representative;
}

template <Int N>
int run_generated_constrained_case(
	int rounds,
	const std::filesystem::path* representatives_output_path,
	const std::filesystem::path* representative_output_path
) {
	using GraphType = OddGraphdegZero<N + 1>;

	std::vector<GraphType> graphs;
	const std::string path = cache_path<N>(rounds);
	if (load_graph_cache(path, graphs)) {
		std::cout << "loaded cache = " << path << '\n';
	} else {
		graphs = generate_same_degree_candidates<N>(rounds, representatives_output_path);
		save_graph_cache(path, graphs);
		std::cout << "saved cache = " << path << '\n';
	}

	std::cout << "wheel = W" << +N << '\n';
	std::cout << "split-contract rounds = " << rounds << '\n';
	std::cout << "same-degree candidates = " << graphs.size() << '\n';

	if (representatives_output_path != nullptr) {
		if (!write_representatives(*representatives_output_path, graphs)) {
			std::cerr << "failed to write " << representatives_output_path->string() << '\n';
			return EXIT_FAILURE;
		}
		std::cout << "saved representatives = " << representatives_output_path->string() << '\n';
	}

	auto representative = solve_generated_constrained_representative<N>(
		graphs,
		representative_output_path
	);
	if (!representative.has_value()) {
		std::cout << "generated constrained representative found = no\n";
		return EXIT_FAILURE;
	}

	std::cout << "generated constrained representative found = yes\n";
	if (representative_output_path != nullptr) {
		if (!write_class_file(*representative_output_path, *representative)) {
			std::cerr << "failed to write " << representative_output_path->string() << '\n';
			return EXIT_FAILURE;
		}
		std::cout << "saved representative = "
		          << representative_output_path->string() << '\n';
	}

	return EXIT_SUCCESS;
}

template <Int N>
int run_constrained_krylov_dual_case(
	int rounds,
	const std::filesystem::path* witness_output_path
) {
	using GCType = OddGCdegZero<N + 1>;
	using GraphType = typename GCType::GraphType;
	using SplitL = typename GCType::SplitL;
	using ContL = typename GCType::ContL;

	typename GraphType::Basis wheel_basis(wheel_graph<N>(), fieldType{1});
	GraphType::std(wheel_basis);
	const GraphType wheel = wheel_basis.getValue();

	std::vector<GCType> krylov_vectors;
	krylov_vectors.reserve(static_cast<std::size_t>(rounds) + 1);
	krylov_vectors.emplace_back(wheel, AssumeBasisOrderTag{});

	std::cout << "wheel = W" << +N << '\n';
	std::cout << "constrained krylov dual depth = " << rounds << '\n';
	std::cout << "krylov step 0 terms = "
	          << krylov_vectors.back().data().size() << '\n';
	for (int step = 1; step <= rounds; ++step) {
		const double seconds = time_seconds([&]() {
			krylov_vectors.push_back(split_contract_krylov_step(krylov_vectors.back()));
		});
		std::cout << "krylov step " << step
		          << " terms = " << krylov_vectors.back().data().size()
		          << ", time = " << seconds << " s\n";
	}

	CombinedSplitContractIndex<GCType> image_index;
	row_based_lil_matrix<fieldType> equations(static_cast<std::size_t>(rounds) + 1);
	for (std::size_t column_index = 1; column_index < krylov_vectors.size(); ++column_index) {
		const auto split_column = krylov_vectors[column_index].delta().data();
		const auto contraction_column = contraction_of(krylov_vectors[column_index]).data();
		std::cout << "column " << column_index
		          << ": d_split terms = " << split_column.size()
		          << ", d_contract terms = " << contraction_column.size() << '\n';

		for (const auto& be : split_column) {
			equations.add_element(
				column_index - 1,
				image_index.split_index(be.getValue()),
				be.getCoefficient()
			);
		}
		for (const auto& be : contraction_column) {
			equations.add_element(
				column_index - 1,
				image_index.contraction_index(be.getValue()),
				be.getCoefficient()
			);
		}
	}

	const std::size_t target_row = static_cast<std::size_t>(rounds);
	const auto split_target = GCType(wheel, AssumeBasisOrderTag{}).delta().data();
	for (const auto& be : split_target) {
		equations.add_element(
			target_row,
			image_index.split_index(be.getValue()),
			be.getCoefficient()
		);
	}
	equations.sort_rows();

	VectorSpace::LinComb<std::size_t, fieldType> rhs;
	rhs.append_in_basis_order(target_row, fieldType{1});

	std::optional<VectorSpace::LinComb<std::size_t, fieldType>> witness_indices;
	const double solve_seconds = time_seconds([&]() {
		witness_indices = equations.solve_MX_equals_y(rhs);
	});

	std::cout << "dual equations = " << (static_cast<std::size_t>(rounds) + 1) << '\n';
	std::cout << "combined image basis = " << image_index.size() << '\n';
	std::cout << "split image basis = " << image_index.split_basis().size() << '\n';
	std::cout << "contraction image basis = "
	          << image_index.contraction_basis().size() << '\n';
	std::cout << "dual matrix entries = " << equations.size() << '\n';
	std::cout << "dual solve time = " << solve_seconds << " s\n";

	if (!witness_indices.has_value()) {
		std::cout << "constrained dual witness found = no\n";
		return EXIT_FAILURE;
	}

	SplitL split_witness;
	ContL contraction_witness;
	for (const auto& be : *witness_indices) {
		const auto& ref = image_index.refs()[be.getValue()];
		if (ref.is_split) {
			split_witness.append_in_basis_order(
				image_index.split_basis()[ref.local_index],
				be.getCoefficient()
			);
		} else {
			contraction_witness.append_in_basis_order(
				image_index.contraction_basis()[ref.local_index],
				be.getCoefficient()
			);
		}
	}
	split_witness.sort_without_deduplicate();
	contraction_witness.sort_without_deduplicate();

	auto evaluate = [](const auto& functional, const auto& vector) {
		fieldType result{};
		std::size_t i = 0;
		std::size_t j = 0;
		while (i < functional.raw_elements().size() && j < vector.raw_elements().size()) {
			const auto& lhs = functional.raw_elements()[i];
			const auto& rhs_term = vector.raw_elements()[j];
			if (lhs.getValue() < rhs_term.getValue()) {
				++i;
			} else if (rhs_term.getValue() < lhs.getValue()) {
				++j;
			} else {
				result += lhs.getCoefficient() * rhs_term.getCoefficient();
				++i;
				++j;
			}
		}
		return result;
	};

	std::size_t nonzero_column_pairings = 0;
	for (std::size_t column_index = 1; column_index < krylov_vectors.size(); ++column_index) {
		const fieldType pairing =
			evaluate(split_witness, krylov_vectors[column_index].delta().data())
			+ evaluate(contraction_witness, contraction_of(krylov_vectors[column_index]).data());
		if (pairing != fieldType{}) {
			++nonzero_column_pairings;
		}
	}
	const fieldType target_pairing = evaluate(split_witness, split_target);

	std::cout << "split witness terms = " << split_witness.size() << '\n';
	std::cout << "contraction witness terms = " << contraction_witness.size() << '\n';
	std::cout << "nonzero pairings with constrained Krylov columns = "
	          << nonzero_column_pairings << '\n';
	std::cout << "pairing with (d_split(W), 0) = " << target_pairing << '\n';
	std::cout << "constrained dual witness found = "
	          << ((nonzero_column_pairings == 0 && target_pairing != fieldType{})
	              ? "yes"
	              : "no") << '\n';

	if (nonzero_column_pairings != 0 || target_pairing == fieldType{}) {
		return EXIT_FAILURE;
	}

	if (witness_output_path != nullptr) {
		typename GCType::SplitGC split_witness_gc(split_witness);
		typename GCType::ContGC contraction_witness_gc(contraction_witness);
		write_class_file(
			related_output_path<typename GCType::SplitGC>(witness_output_path, "_split"),
			split_witness_gc
		);
		write_class_file(
			related_output_path<typename GCType::ContGC>(witness_output_path, "_contract"),
			contraction_witness_gc
		);
		std::cout << "saved split witness = "
		          << related_output_path<typename GCType::SplitGC>(witness_output_path, "_split").string()
		          << '\n';
		std::cout << "saved contraction witness = "
		          << related_output_path<typename GCType::ContGC>(witness_output_path, "_contract").string()
		          << '\n';
	}

	return EXIT_SUCCESS;
}

void print_usage(const char* argv0) {
	std::cerr << "usage: " << argv0 << " wheel_N rounds [representatives_output_path] [cycle_output_path] [mode]\n";
	std::cerr << "  mode defaults to generated-support; use generated-constrained, krylov, krylov-dual, krylov-constrained, or krylov-constrained-dual\n";
	std::cerr << "  wheel_N must be one of 3,5,7,9,11,13,15,17,19,21,25,27\n";
}

} // namespace

int main(int argc, char** argv) {
	if (argc < 3 || argc > 6) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	const int wheel_n = std::stoi(argv[1]);
	const int rounds = std::stoi(argv[2]);
	if (rounds < 0) {
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

	std::filesystem::path representatives_output_path;
	const std::filesystem::path* representatives_output_path_ptr = nullptr;
	if (argc >= 4) {
		representatives_output_path = argv[3];
		representatives_output_path_ptr = &representatives_output_path;
	}

	std::filesystem::path cycle_output_path;
	const std::filesystem::path* cycle_output_path_ptr = nullptr;
	if (argc >= 5) {
		cycle_output_path = argv[4];
		cycle_output_path_ptr = &cycle_output_path;
	}

	ExperimentMode mode = ExperimentMode::GeneratedSupport;
	if (argc >= 6) {
		const std::string mode_arg = argv[5];
		if (mode_arg == "krylov") {
			mode = ExperimentMode::Krylov;
		} else if (mode_arg == "krylov-dual") {
			mode = ExperimentMode::KrylovDual;
		} else if (mode_arg == "krylov-constrained") {
			mode = ExperimentMode::KrylovConstrained;
		} else if (mode_arg == "krylov-constrained-dual") {
			mode = ExperimentMode::KrylovConstrainedDual;
		} else if (mode_arg == "generated-constrained") {
			mode = ExperimentMode::GeneratedConstrained;
		} else if (mode_arg == "generated-support") {
			mode = ExperimentMode::GeneratedSupport;
		} else {
			print_usage(argv[0]);
			return EXIT_FAILURE;
		}
	}

#define RUN_SELECTED_CASE(N_VALUE) \
	do { \
		if (mode == ExperimentMode::GeneratedConstrained) { \
			return run_generated_constrained_case<N_VALUE>(rounds, representatives_output_path_ptr, cycle_output_path_ptr); \
		} \
		if (mode == ExperimentMode::KrylovConstrainedDual) { \
			return run_constrained_krylov_dual_case<N_VALUE>(rounds, cycle_output_path_ptr); \
		} \
		if (mode == ExperimentMode::KrylovConstrained) { \
			return run_constrained_krylov_case<N_VALUE>(rounds, cycle_output_path_ptr); \
		} \
		if (mode == ExperimentMode::KrylovDual) { \
			return run_krylov_dual_case<N_VALUE>(rounds, cycle_output_path_ptr); \
		} \
		if (mode == ExperimentMode::Krylov) { \
			return run_krylov_case<N_VALUE>(rounds, cycle_output_path_ptr); \
		} \
		return run_case<N_VALUE>(rounds, representatives_output_path_ptr, cycle_output_path_ptr); \
	} while (false)

	switch (wheel_n) {
	case 3:
		RUN_SELECTED_CASE(3);
	case 5:
		RUN_SELECTED_CASE(5);
	case 7:
		RUN_SELECTED_CASE(7);
	case 9:
		RUN_SELECTED_CASE(9);
	case 11:
		RUN_SELECTED_CASE(11);
	case 13:
		RUN_SELECTED_CASE(13);
	case 15:
		RUN_SELECTED_CASE(15);
	case 17:
		RUN_SELECTED_CASE(17);
	case 19:
		RUN_SELECTED_CASE(19);
	case 21:
		RUN_SELECTED_CASE(21);
	case 25:
		RUN_SELECTED_CASE(25);
	case 27:
		RUN_SELECTED_CASE(27);
	default:
		print_usage(argv[0]);
		return EXIT_FAILURE;
	}

#undef RUN_SELECTED_CASE
}
