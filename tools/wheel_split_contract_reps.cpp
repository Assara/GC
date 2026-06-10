#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "examplegraphs.hpp"
#include "GraphStandardizer.hpp"
#include "VectorSpace/row_based_lil_matrix.hpp"

namespace {

enum class ExperimentMode {
	Enumerate,
	GeneratedSupport,
	GeneratedConstrained,
	SplitMap,
	SplitMapInfo,
	SplitMapSolve,
	SplitMapDual,
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
std::filesystem::path split_map_prefix(int rounds) {
	return std::filesystem::path("output")
		/ "split_contract"
		/ "maps"
		/ ("W" + std::to_string(N)
			+ "_rounds" + std::to_string(rounds)
			+ "_split_map");
}

template <typename GraphType>
bool load_graph_records(const std::string& path, std::vector<GraphType>& graphs) {
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
void save_graph_records(const std::string& path, const std::vector<GraphType>& graphs) {
	const std::filesystem::path output_file(path);
	if (output_file.has_parent_path()) {
		std::filesystem::create_directories(output_file.parent_path());
	}

	std::ofstream out(path, std::ios::binary);
	if (!out) {
		std::cerr << "warning: failed to write graph records " << path << '\n';
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

template <typename T>
void write_binary_value(std::ofstream& out, const T& value) {
	out.write(reinterpret_cast<const char*>(&value), sizeof(value));
}

template <typename T>
bool read_binary_value(std::ifstream& in, T& value) {
	in.read(reinterpret_cast<char*>(&value), sizeof(value));
	return static_cast<bool>(in);
}

template <typename GraphType>
void write_graph_record(std::ofstream& out, const GraphType& graph) {
	out.write(
		reinterpret_cast<const char*>(graph.half_edges.data()),
		static_cast<std::streamsize>(graph.half_edges.size() * sizeof(Int))
	);
}

template <typename GraphType>
bigInt automorphism_group_size4(const GraphType& input_graph) {
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

void write_field_vector(std::ofstream& out, const fieldType* values, std::size_t size) {
	for (std::size_t i = 0; i < size; ++i) {
		fieldType::write_value(out, values[i]);
	}
}

bool read_field_vector(std::ifstream& in, fieldType* values, std::size_t size) {
	for (std::size_t i = 0; i < size; ++i) {
		if (!fieldType::read_value(in, values[i])) {
			return false;
		}
	}
	return true;
}

struct SplitMapHeader {
	std::uint32_t wheel_n = 0;
	std::uint32_t rounds = 0;
	std::uint32_t int_size = 0;
	std::uint32_t coeff_size = 0;
	std::uint32_t domain_vertices = 0;
	std::uint32_t domain_edges = 0;
	std::uint32_t image_vertices = 0;
	std::uint32_t image_edges = 0;
	std::uint64_t total_graphs = 0;
	std::uint64_t domain_dim = 0;
	std::uint64_t image_dim = 0;
	std::uint64_t entries = 0;
};

std::optional<compressed_sparse_matrix<fieldType>>
load_split_map_columns_as_compressed_matrix(
	const std::filesystem::path& columns_path,
	SplitMapHeader* header_out = nullptr
) {
	std::ifstream in(columns_path, std::ios::binary);
	if (!in) {
		std::cerr << "failed to read " << columns_path.string() << '\n';
		return std::nullopt;
	}

	char magic[16]{};
	in.read(magic, sizeof(magic));
	const char expected_magic[16] = {
		'G','C','S','P','L','I','T','M','A','P','C','O','L','1','\0','\0'
	};
	if (!in || std::string_view(magic, sizeof(magic)) != std::string_view(expected_magic, sizeof(expected_magic))) {
		std::cerr << "invalid split-map column file magic: "
		          << columns_path.string() << '\n';
		return std::nullopt;
	}

	SplitMapHeader header;
	if (!read_binary_value(in, header.wheel_n)
		|| !read_binary_value(in, header.rounds)
		|| !read_binary_value(in, header.int_size)
		|| !read_binary_value(in, header.coeff_size)
		|| !read_binary_value(in, header.domain_vertices)
		|| !read_binary_value(in, header.domain_edges)
		|| !read_binary_value(in, header.image_vertices)
		|| !read_binary_value(in, header.image_edges)
		|| !read_binary_value(in, header.total_graphs)
		|| !read_binary_value(in, header.domain_dim)
		|| !read_binary_value(in, header.image_dim)
		|| !read_binary_value(in, header.entries)) {
		std::cerr << "failed to read split-map header: "
		          << columns_path.string() << '\n';
		return std::nullopt;
	}

	if (header.int_size != sizeof(Int)
		|| header.coeff_size != fieldType::serialized_value_size_hint()) {
		std::cerr << "split-map storage size mismatch in "
		          << columns_path.string() << '\n';
		return std::nullopt;
	}
	if (header.image_dim > std::numeric_limits<compressed_sparse_matrix<fieldType>::indexType>::max()
		|| header.domain_dim > std::numeric_limits<compressed_sparse_matrix<fieldType>::indexType>::max()) {
		std::cerr << "split-map dimensions exceed compressed_sparse_matrix index type\n";
		return std::nullopt;
	}

	compressed_sparse_matrix<fieldType> matrix(
		static_cast<compressed_sparse_matrix<fieldType>::indexType>(header.image_dim)
	);
	matrix.rows_and_coeffs_.reserve(static_cast<std::size_t>(header.entries));
	matrix.col_ptr_.reserve(static_cast<std::size_t>(header.domain_dim) + 1);

	for (std::uint64_t column = 0; column < header.domain_dim; ++column) {
		std::uint32_t nnz = 0;
		if (!read_binary_value(in, nnz)) {
			std::cerr << "failed to read column size at column " << column << '\n';
			return std::nullopt;
		}
		for (std::uint32_t entry = 0; entry < nnz; ++entry) {
			std::uint32_t row = 0;
			typename fieldType::serialized_value_type coefficient{};
			if (!read_binary_value(in, row) || !fieldType::read_serialized_value(in, coefficient)) {
				std::cerr << "failed to read map entry at column "
				          << column << '\n';
				return std::nullopt;
			}
			matrix.rows_and_coeffs_.emplace_back(row, fieldType{std::move(coefficient)});
		}
		matrix.col_ptr_.push_back(
			static_cast<compressed_sparse_matrix<fieldType>::offset_type>(
				matrix.rows_and_coeffs_.size()
			)
		);
	}

	if (matrix.rows_and_coeffs_.size() != static_cast<std::size_t>(header.entries)) {
		std::cerr << "split-map entry count mismatch: expected "
		          << header.entries << ", read "
		          << matrix.rows_and_coeffs_.size() << '\n';
		return std::nullopt;
	}

	if (header_out != nullptr) {
		*header_out = header;
	}
	return matrix;
}

compressed_sparse_matrix<fieldType> build_scaled_adjoint_matrix(
	const compressed_sparse_matrix<fieldType>& matrix,
	const std::vector<fieldType>& image_aut,
	const std::vector<fieldType>& domain_aut
) {
	using Matrix = compressed_sparse_matrix<fieldType>;
	using Basis = typename Matrix::Basis;

	std::vector<std::vector<Basis>> columns(static_cast<std::size_t>(matrix.image_dim()));
	for (std::uint32_t domain_col = 0; domain_col < matrix.domain_dim(); ++domain_col) {
		const fieldType domain_scale = domain_aut[domain_col];
		for (const auto& be : matrix.get_column(domain_col)) {
			const std::size_t image_row = static_cast<std::size_t>(be.getValue());
			columns[image_row].emplace_back(
				domain_col,
				be.getCoefficient() * image_aut[image_row] / domain_scale
			);
		}
	}

	Matrix adjoint(matrix.domain_dim());
	for (auto& column : columns) {
		adjoint.add_col(column);
	}
	return adjoint;
}

template <typename SparseColumn>
std::vector<typename compressed_sparse_matrix<fieldType>::Basis>
apply_compressed_matrix_to_sparse_column(
	const compressed_sparse_matrix<fieldType>& matrix,
	const SparseColumn& column
) {
	using Matrix = compressed_sparse_matrix<fieldType>;
	using Basis = typename Matrix::Basis;

	std::unordered_map<std::uint32_t, fieldType> accum;
	for (const auto& be : column) {
		const fieldType scale = be.getCoefficient();
		for (const auto& image_be : matrix.get_column(static_cast<std::uint32_t>(be.getValue()))) {
			accum[image_be.getValue()] += image_be.getCoefficient() * scale;
		}
	}

	std::vector<Basis> result;
	result.reserve(accum.size());
	for (const auto& [index, coefficient] : accum) {
		if (coefficient != fieldType{}) {
			result.emplace_back(index, coefficient);
		}
	}
	std::sort(result.begin(), result.end());
	return result;
}

template <typename GCType>
typename GCType::ContGC contraction_of(GCType gamma);

compressed_sparse_matrix<fieldType> transpose_compressed_matrix(
	const compressed_sparse_matrix<fieldType>& matrix
) {
	using Matrix = compressed_sparse_matrix<fieldType>;
	using Index = typename Matrix::indexType;
	using Offset = typename Matrix::offset_type;

	Matrix transpose(matrix.domain_dim());
	const std::size_t transpose_domain_dim = matrix.image_dim();
	transpose.col_ptr_.assign(transpose_domain_dim + 1, Offset{0});

	for (const auto& be : matrix.rows_and_coeffs_) {
		++transpose.col_ptr_[static_cast<std::size_t>(be.getValue()) + 1];
	}
	for (std::size_t i = 1; i < transpose.col_ptr_.size(); ++i) {
		transpose.col_ptr_[i] += transpose.col_ptr_[i - 1];
	}

	transpose.rows_and_coeffs_.resize(matrix.rows_and_coeffs_.size());
	std::vector<Offset> cursor = transpose.col_ptr_;
	for (Index col = 0; col < matrix.domain_dim(); ++col) {
		for (const auto& be : matrix.get_column(col)) {
			const std::size_t transpose_col = static_cast<std::size_t>(be.getValue());
			const Offset write = cursor[transpose_col]++;
			transpose.rows_and_coeffs_[write] = typename Matrix::Basis(col, be.getCoefficient());
		}
	}

	return transpose;
}

class compressed_wiedemann_solver {
public:
	using Matrix = compressed_sparse_matrix<fieldType>;
	using DenseImageVec = typename Matrix::DenseImageVec;
	using DenseDomainVec = typename Matrix::DenseDomainVec;

	explicit compressed_wiedemann_solver(Matrix matrix)
		: matrix_(std::move(matrix)),
		  transpose_matrix_(transpose_compressed_matrix(matrix_)),
		  random_vector_(make_random_image_vec()) {}

	compressed_wiedemann_solver(Matrix matrix, Matrix transpose_matrix)
		: matrix_(std::move(matrix)),
		  transpose_matrix_(std::move(transpose_matrix)),
		  random_vector_(make_random_image_vec()) {}

	std::optional<DenseDomainVec> solve_MX_equals_y(
		const DenseImageVec& y0,
		const std::filesystem::path* checkpoint_path = nullptr,
		std::size_t checkpoint_interval = 0
	) const {
		if (matrix_.image_dim() == 0) {
			return make_dense_domain_vec_zero();
		}

		const bool checkpointing_enabled =
			checkpoint_path != nullptr && checkpoint_interval > 0;

		std::cout << "using compressed split-map solver. image dim = "
		          << matrix_.image_dim()
		          << " domain dim = " << matrix_.domain_dim() << '\n';
		if (checkpointing_enabled) {
			std::cout << "checkpoint path = " << checkpoint_path->string() << '\n';
			std::cout << "checkpoint interval = " << checkpoint_interval << '\n';
		} else {
			std::cout << "checkpointing disabled\n";
		}

		bool is_zero = true;
		for (std::size_t i = 0; i < matrix_.image_dim(); ++i) {
			if (y0[i] != fieldType{}) {
				is_zero = false;
				break;
			}
		}
		if (is_zero) {
			return make_dense_domain_vec_zero();
		}

		const std::size_t slack = 8;
		std::vector<fieldType> signatures;
		signatures.reserve(2 * static_cast<std::size_t>(matrix_.image_dim()) + slack);

		DenseDomainVec mt_work = make_dense_domain_vec_zero();
		DenseImageVec yi = make_dense_image_vec_zero();
		DenseImageVec next_yi = make_dense_image_vec_zero();
		std::size_t last_checkpoint_iteration = 0;
		auto checkpoint_timer = std::chrono::steady_clock::now();

		if (checkpointing_enabled && load_checkpoint(*checkpoint_path, signatures, yi)) {
			last_checkpoint_iteration = latest_mmt_iteration(signatures);
			std::cout << "resumed Wiedemann checkpoint = "
			          << checkpoint_path->string()
			          << ", signatures = " << signatures.size()
			          << ", MMT iteration = " << latest_mmt_iteration(signatures)
			          << '\n';
		} else {
			signatures.emplace_back(get_signature(y0));
			M_MT_into(y0, mt_work, yi);
			signatures.emplace_back(get_signature(yi));
			if (checkpointing_enabled) {
				const double checkpoint_seconds = time_seconds([&]() {
					save_checkpoint(*checkpoint_path, signatures, yi);
				});
				last_checkpoint_iteration = latest_mmt_iteration(signatures);
				checkpoint_timer = std::chrono::steady_clock::now();
				std::cout << "checkpointed Wiedemann signatures = "
				          << signatures.size()
				          << ", MMT iteration = " << latest_mmt_iteration(signatures)
				          << ", recurrence order = initial"
				          << ", seconds since previous checkpoint = 0"
				          << ", checkpoint write seconds = " << checkpoint_seconds
				          << '\n';
			}
		}

		berlekamp_massey_state<fieldType> bm_state(signatures);
		bm_state.process_all_new();

		std::size_t more_needed = bm_state.more_needed(slack);
		while (more_needed > 0) {
			const std::size_t chunk = std::min<std::size_t>(
				more_needed,
				checkpointing_enabled ? checkpoint_interval : more_needed
			);
			for (std::size_t i = 0; i < chunk; ++i) {
				M_MT_into(yi, mt_work, next_yi);
				std::swap(yi, next_yi);
				signatures.emplace_back(get_signature(yi));
			}
			bm_state.process_all_new();
			const std::size_t current_iteration = latest_mmt_iteration(signatures);
			if (checkpointing_enabled
				&& current_iteration - last_checkpoint_iteration >= checkpoint_interval) {
				const auto now = std::chrono::steady_clock::now();
				const double seconds_since_checkpoint =
					std::chrono::duration<double>(now - checkpoint_timer).count();
				const double checkpoint_seconds = time_seconds([&]() {
					save_checkpoint(*checkpoint_path, signatures, yi);
				});
				last_checkpoint_iteration = current_iteration;
				checkpoint_timer = std::chrono::steady_clock::now();
				std::cout << "checkpointed Wiedemann signatures = "
				          << signatures.size()
				          << ", MMT iteration = " << current_iteration
				          << ", recurrence order = " << bm_state.order()
				          << ", seconds since previous checkpoint = "
				          << seconds_since_checkpoint
				          << ", checkpoint write seconds = " << checkpoint_seconds
				          << '\n';
			}
			more_needed = bm_state.more_needed(slack);
		}

		return recompute_solution(
			y0,
			bm_state.connection_poly(),
			checkpoint_path,
			checkpoint_interval
		);
	}

	const Matrix& matrix() const {
		return matrix_;
	}

private:
	Matrix matrix_;
	Matrix transpose_matrix_;
	DenseImageVec random_vector_;

	DenseImageVec make_dense_image_vec_zero() const {
		auto result = std::make_unique<fieldType[]>(static_cast<std::size_t>(matrix_.image_dim()));
		for (std::size_t i = 0; i < matrix_.image_dim(); ++i) {
			result[i] = fieldType{};
		}
		return result;
	}

	DenseDomainVec make_dense_domain_vec_zero() const {
		auto result = std::make_unique<fieldType[]>(static_cast<std::size_t>(matrix_.domain_dim()));
		for (std::size_t i = 0; i < matrix_.domain_dim(); ++i) {
			result[i] = fieldType{};
		}
		return result;
	}

	DenseImageVec make_random_image_vec() const {
		auto result = std::make_unique<fieldType[]>(static_cast<std::size_t>(matrix_.image_dim()));
		std::mt19937_64 rng(0xC0FFEEULL);
		for (std::size_t i = 0; i < matrix_.image_dim(); ++i) {
			result[i] = fieldType::sample(rng);
		}
		return result;
	}

	static std::size_t latest_mmt_iteration(const std::vector<fieldType>& signatures) {
		return signatures.empty() ? 0 : signatures.size() - 1;
	}

	void M_MT_into(const DenseImageVec& y, DenseDomainVec& v, DenseImageVec& result) const {
		matrix_.evaluate_transpose(y, v);
		transpose_matrix_.evaluate_transpose(v, result);
	}

	bool load_checkpoint(
		const std::filesystem::path& path,
		std::vector<fieldType>& signatures,
		DenseImageVec& yi
	) const {
		std::ifstream in(path, std::ios::binary);
		if (!in) {
			return false;
		}

		char magic[16]{};
		in.read(magic, sizeof(magic));
		const char expected_magic[16] = {
			'G','C','W','I','E','D','E','M','A','N','N','C','K','2','\0','\0'
		};
		if (!in || std::string_view(magic, sizeof(magic)) != std::string_view(expected_magic, sizeof(expected_magic))) {
			std::cerr << "ignoring invalid checkpoint = " << path.string() << '\n';
			return false;
		}

		std::uint64_t image_dim = 0;
		std::uint64_t domain_dim = 0;
		std::uint64_t signature_count = 0;
		if (!read_binary_value(in, image_dim)
			|| !read_binary_value(in, domain_dim)
			|| !read_binary_value(in, signature_count)) {
			std::cerr << "ignoring incomplete checkpoint header = " << path.string() << '\n';
			return false;
		}

		if (image_dim != matrix_.image_dim()
			|| domain_dim != matrix_.domain_dim()
			|| signature_count < 2) {
			std::cerr << "ignoring incompatible checkpoint = " << path.string() << '\n';
			return false;
		}

		signatures.assign(static_cast<std::size_t>(signature_count), fieldType{});
		DenseImageVec saved_random_vector = make_dense_image_vec_zero();
		if (!read_field_vector(in, signatures.data(), signatures.size())
			|| !read_field_vector(in, saved_random_vector.get(), static_cast<std::size_t>(matrix_.image_dim()))
			|| !read_field_vector(in, yi.get(), static_cast<std::size_t>(matrix_.image_dim()))) {
			std::cerr << "ignoring incomplete checkpoint payload = " << path.string() << '\n';
			signatures.clear();
			return false;
		}
		for (std::size_t i = 0; i < matrix_.image_dim(); ++i) {
			if (saved_random_vector[i] != random_vector_[i]) {
				std::cerr << "ignoring checkpoint with incompatible random vector = "
				          << path.string() << '\n';
				signatures.clear();
				return false;
			}
		}

		return true;
	}

	void save_checkpoint(
		const std::filesystem::path& path,
		const std::vector<fieldType>& signatures,
		const DenseImageVec& yi
	) const {
		if (path.has_parent_path()) {
			std::filesystem::create_directories(path.parent_path());
		}
		const std::filesystem::path tmp_path = path.string() + ".tmp";
		std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
		if (!out) {
			std::cerr << "warning: failed to write checkpoint = "
			          << tmp_path.string() << '\n';
			return;
		}

		const char magic[16] = {
			'G','C','W','I','E','D','E','M','A','N','N','C','K','2','\0','\0'
		};
		out.write(magic, sizeof(magic));
		const std::uint64_t image_dim = matrix_.image_dim();
		const std::uint64_t domain_dim = matrix_.domain_dim();
		const std::uint64_t signature_count = signatures.size();
		write_binary_value(out, image_dim);
		write_binary_value(out, domain_dim);
		write_binary_value(out, signature_count);
		write_field_vector(out, signatures.data(), signatures.size());
		write_field_vector(out, random_vector_.get(), static_cast<std::size_t>(matrix_.image_dim()));
		write_field_vector(out, yi.get(), static_cast<std::size_t>(matrix_.image_dim()));
		out.close();

		if (!out) {
			std::cerr << "warning: failed while writing checkpoint = "
			          << tmp_path.string() << '\n';
			return;
		}

		std::error_code ec;
		std::filesystem::rename(tmp_path, path, ec);
		if (ec) {
			std::filesystem::remove(path, ec);
			ec.clear();
			std::filesystem::rename(tmp_path, path, ec);
		}
		if (ec) {
			std::cerr << "warning: failed to install checkpoint = "
			          << path.string() << ": " << ec.message() << '\n';
		}
	}

	bool load_recompute_checkpoint(
		const std::filesystem::path& path,
		const std::vector<fieldType>& connection_poly,
		std::size_t& next_step,
		DenseImageVec& acc
	) const {
		std::ifstream in(path, std::ios::binary);
		if (!in) {
			return false;
		}

		char magic[16]{};
		in.read(magic, sizeof(magic));
		const char expected_magic[16] = {
			'G','C','W','I','E','D','R','E','C','O','M','P','2','\0','\0','\0'
		};
		if (!in || std::string_view(magic, sizeof(magic)) != std::string_view(expected_magic, sizeof(expected_magic))) {
			std::cerr << "ignoring invalid recomputation checkpoint = "
			          << path.string() << '\n';
			return false;
		}

		std::uint64_t image_dim = 0;
		std::uint64_t domain_dim = 0;
		std::uint64_t saved_next_step = 0;
		std::uint64_t poly_size = 0;
		if (!read_binary_value(in, image_dim)
			|| !read_binary_value(in, domain_dim)
			|| !read_binary_value(in, saved_next_step)
			|| !read_binary_value(in, poly_size)) {
			std::cerr << "ignoring incomplete recomputation checkpoint header = "
			          << path.string() << '\n';
			return false;
		}

		if (image_dim != matrix_.image_dim()
			|| domain_dim != matrix_.domain_dim()
			|| poly_size != connection_poly.size()
			|| saved_next_step == 0
			|| saved_next_step > connection_poly.size()) {
			std::cerr << "ignoring incompatible recomputation checkpoint = "
			          << path.string() << '\n';
			return false;
		}

		std::vector<fieldType> saved_poly(static_cast<std::size_t>(poly_size));
		if (!read_field_vector(in, saved_poly.data(), saved_poly.size())
			|| !read_field_vector(in, acc.get(), static_cast<std::size_t>(matrix_.image_dim()))) {
			std::cerr << "ignoring incomplete recomputation checkpoint payload = "
			          << path.string() << '\n';
			return false;
		}

		if (saved_poly != connection_poly) {
			std::cerr << "ignoring recomputation checkpoint with different connection polynomial = "
			          << path.string() << '\n';
			return false;
		}

		next_step = static_cast<std::size_t>(saved_next_step);
		return true;
	}

	void save_recompute_checkpoint(
		const std::filesystem::path& path,
		const std::vector<fieldType>& connection_poly,
		std::size_t next_step,
		const DenseImageVec& acc
	) const {
		if (path.has_parent_path()) {
			std::filesystem::create_directories(path.parent_path());
		}
		const std::filesystem::path tmp_path = path.string() + ".tmp";
		std::ofstream out(tmp_path, std::ios::binary | std::ios::trunc);
		if (!out) {
			std::cerr << "warning: failed to write recomputation checkpoint = "
			          << tmp_path.string() << '\n';
			return;
		}

		const char magic[16] = {
			'G','C','W','I','E','D','R','E','C','O','M','P','2','\0','\0','\0'
		};
		out.write(magic, sizeof(magic));
		const std::uint64_t image_dim = matrix_.image_dim();
		const std::uint64_t domain_dim = matrix_.domain_dim();
		const std::uint64_t saved_next_step = next_step;
		const std::uint64_t poly_size = connection_poly.size();
		write_binary_value(out, image_dim);
		write_binary_value(out, domain_dim);
		write_binary_value(out, saved_next_step);
		write_binary_value(out, poly_size);
		write_field_vector(out, connection_poly.data(), connection_poly.size());
		write_field_vector(out, acc.get(), static_cast<std::size_t>(matrix_.image_dim()));
		out.close();

		if (!out) {
			std::cerr << "warning: failed while writing recomputation checkpoint = "
			          << tmp_path.string() << '\n';
			return;
		}

		std::error_code ec;
		std::filesystem::rename(tmp_path, path, ec);
		if (ec) {
			std::filesystem::remove(path, ec);
			ec.clear();
			std::filesystem::rename(tmp_path, path, ec);
		}
		if (ec) {
			std::cerr << "warning: failed to install recomputation checkpoint = "
			          << path.string() << ": " << ec.message() << '\n';
		}
	}

	fieldType get_signature(const DenseImageVec& y) const {
		fieldType result{};
		for (std::size_t i = 0; i < matrix_.image_dim(); ++i) {
			result += y[i] * random_vector_[i];
		}
		return result;
	}

	static void add_scaled_inplace(
		DenseImageVec& dst,
		const DenseImageVec& src,
		std::size_t size,
		fieldType scalar
	) {
		for (std::size_t i = 0; i < size; ++i) {
			dst[i] += src[i] * scalar;
		}
	}

	template <typename DenseVec>
	static DenseVec scaled_copy(const DenseVec& vec, std::size_t size, fieldType scalar) {
		DenseVec result = std::make_unique<fieldType[]>(size);
		for (std::size_t i = 0; i < size; ++i) {
			result[i] = vec[i] * scalar;
		}
		return result;
	}

	std::optional<DenseDomainVec> recompute_solution(
		const DenseImageVec& y0,
		const std::vector<fieldType>& connection_poly,
		const std::filesystem::path* checkpoint_path,
		std::size_t checkpoint_interval
	) const {
		std::cout << "recomputing solution connection_poly.size() = "
		          << connection_poly.size() << '\n';

		if (connection_poly.size() < 2 || connection_poly.front() == fieldType{}) {
			return std::nullopt;
		}
		const std::size_t length = connection_poly.size() - 1;
		const fieldType constant = connection_poly.back();
		if (constant == fieldType{}) {
			std::cout << "singular\n";
			return std::nullopt;
		}

		const fieldType scale = -(constant.inv());
		DenseImageVec acc = scaled_copy(y0, matrix_.image_dim(), fieldType{1});
		DenseImageVec next_acc = make_dense_image_vec_zero();
		DenseDomainVec mt_work = make_dense_domain_vec_zero();
		std::size_t start_j = 1;

		const bool checkpointing_enabled =
			checkpoint_path != nullptr && checkpoint_interval > 0;
		const std::filesystem::path recompute_checkpoint_path =
			checkpointing_enabled
				? std::filesystem::path(checkpoint_path->string() + ".recompute")
				: std::filesystem::path{};
		if (checkpointing_enabled
			&& load_recompute_checkpoint(
				recompute_checkpoint_path,
				connection_poly,
				start_j,
				acc
			)) {
			std::cout << "resumed recomputation checkpoint = "
			          << recompute_checkpoint_path.string()
			          << ", next step = " << start_j << '\n';
		}
		std::size_t last_recompute_checkpoint_step = start_j > 0 ? start_j - 1 : 0;
		auto recompute_checkpoint_timer = std::chrono::steady_clock::now();

		for (std::size_t j = start_j; j <= length - 1; ++j) {
			M_MT_into(acc, mt_work, next_acc);
			std::swap(acc, next_acc);
			add_scaled_inplace(acc, y0, matrix_.image_dim(), connection_poly[j]);
			if (checkpointing_enabled
				&& (j - last_recompute_checkpoint_step >= checkpoint_interval
					|| j == length - 1)) {
				const auto now = std::chrono::steady_clock::now();
				const double seconds_since_checkpoint =
					std::chrono::duration<double>(now - recompute_checkpoint_timer).count();
				const double checkpoint_seconds = time_seconds([&]() {
					save_recompute_checkpoint(
						recompute_checkpoint_path,
						connection_poly,
						j + 1,
						acc
					);
				});
				last_recompute_checkpoint_step = j;
				recompute_checkpoint_timer = std::chrono::steady_clock::now();
				std::cout << "checkpointed recomputation next step = "
				          << (j + 1) << " / " << length
				          << ", recomputation MMT iterations completed = " << j
				          << ", seconds since previous checkpoint = "
				          << seconds_since_checkpoint
				          << ", checkpoint write seconds = " << checkpoint_seconds
				          << '\n';
			}
		}

		DenseImageVec y_minus_1 = scaled_copy(acc, matrix_.image_dim(), scale);
		DenseDomainVec result = make_dense_domain_vec_zero();
		matrix_.evaluate_transpose(y_minus_1, result);

		const fieldType lambda = find_scaling_and_verify(y0, result);
		if (lambda == fieldType{}) {
			return std::nullopt;
		}

		return scaled_copy(result, matrix_.domain_dim(), lambda);
	}

	fieldType find_scaling_and_verify(const DenseImageVec& y0, const DenseDomainVec& x) const {
		DenseImageVec y_found = make_dense_image_vec_zero();
		transpose_matrix_.evaluate_transpose(x, y_found);

		fieldType lambda{};
		std::size_t i = 0;
		for (; i < matrix_.image_dim(); ++i) {
			if (y0[i] != fieldType{} && y_found[i] != fieldType{}) {
				lambda = y0[i] / y_found[i];
				break;
			}
			if (y0[i] != fieldType{} || y_found[i] != fieldType{}) {
				std::cout << "false solution type 1\n";
				return fieldType{};
			}
		}

		for (; i < matrix_.image_dim(); ++i) {
			if (y0[i] != lambda * y_found[i]) {
				std::cout << "false solution type 2\n";
				return fieldType{};
			}
		}

		std::cout << "Solution validated! with lambda " << lambda << '\n';
		return lambda;
	}
};

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
	out << "field_type: " << fieldType::name() << "\n";
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
	out << "field_type: " << fieldType::name() << "\n";
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
	graphs = generate_same_degree_candidates<N>(rounds, representatives_output_path);

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

template <Int N>
int run_enumerate_case(
	int rounds,
	const std::filesystem::path* representatives_output_path
) {
	auto graphs = generate_same_degree_candidates<N>(rounds, representatives_output_path);

	std::cout << "wheel = W" << +N << '\n';
	std::cout << "split-contract rounds = " << rounds << '\n';
	std::cout << "same-degree candidates = " << graphs.size() << '\n';
	if (representatives_output_path != nullptr) {
		if (!write_representatives(*representatives_output_path, graphs)) {
			std::cerr << "failed to write " << representatives_output_path->string() << '\n';
			return EXIT_FAILURE;
		}
		std::cout << "saved representatives = "
		          << representatives_output_path->string() << '\n';
	}

	return EXIT_SUCCESS;
}

template <Int N>
int run_split_map_export_case(
	int rounds,
	const std::filesystem::path* representatives_output_path,
	const std::filesystem::path* map_output_prefix
) {
	using GCType = OddGCdegZero<N + 1>;
	using GraphType = typename GCType::GraphType;
	using SplitGraph = typename GCType::SplitGraphType;

	auto graphs = generate_same_degree_candidates<N>(rounds, representatives_output_path);

	typename GraphType::Basis wheel_basis(wheel_graph<N>(), fieldType{1});
	GraphType::std(wheel_basis);
	const GraphType wheel = wheel_basis.getValue();

	const std::filesystem::path prefix = map_output_prefix != nullptr
		? *map_output_prefix
		: split_map_prefix<N>(rounds);
	const std::filesystem::path columns_path = prefix.string() + "_columns.bin";
	const std::filesystem::path image_basis_path = prefix.string() + "_image_basis.bin";
	const std::filesystem::path metadata_path = prefix.string() + "_metadata.txt";

	if (columns_path.has_parent_path()) {
		std::filesystem::create_directories(columns_path.parent_path());
	}

	std::ofstream columns_out(columns_path, std::ios::binary | std::ios::trunc);
	if (!columns_out) {
		std::cerr << "failed to write " << columns_path.string() << '\n';
		return EXIT_FAILURE;
	}

	const char magic[16] = {
		'G','C','S','P','L','I','T','M','A','P','C','O','L','1','\0','\0'
	};
	columns_out.write(magic, sizeof(magic));
	const std::uint32_t wheel_n_u32 = static_cast<std::uint32_t>(N);
	const std::uint32_t rounds_u32 = static_cast<std::uint32_t>(rounds);
	const std::uint32_t int_size_u32 = static_cast<std::uint32_t>(sizeof(Int));
	const std::uint32_t coeff_size_u32 = fieldType::serialized_value_size_hint();
	const std::uint32_t domain_vertices_u32 = static_cast<std::uint32_t>(GraphType::N_VERTICES_);
	const std::uint32_t domain_edges_u32 = static_cast<std::uint32_t>(GraphType::N_EDGES_);
	const std::uint32_t image_vertices_u32 = static_cast<std::uint32_t>(SplitGraph::N_VERTICES_);
	const std::uint32_t image_edges_u32 = static_cast<std::uint32_t>(SplitGraph::N_EDGES_);
	const std::uint64_t total_graphs_u64 = static_cast<std::uint64_t>(graphs.size());
	std::uint64_t domain_dim_u64 = 0;
	std::uint64_t image_dim_u64 = 0;
	std::uint64_t entries_u64 = 0;

	write_binary_value(columns_out, wheel_n_u32);
	write_binary_value(columns_out, rounds_u32);
	write_binary_value(columns_out, int_size_u32);
	write_binary_value(columns_out, coeff_size_u32);
	write_binary_value(columns_out, domain_vertices_u32);
	write_binary_value(columns_out, domain_edges_u32);
	write_binary_value(columns_out, image_vertices_u32);
	write_binary_value(columns_out, image_edges_u32);
	write_binary_value(columns_out, total_graphs_u64);
	const std::streampos domain_dim_pos = columns_out.tellp();
	write_binary_value(columns_out, domain_dim_u64);
	const std::streampos image_dim_pos = columns_out.tellp();
	write_binary_value(columns_out, image_dim_u64);
	const std::streampos entries_pos = columns_out.tellp();
	write_binary_value(columns_out, entries_u64);

	std::unordered_map<SplitGraph, std::uint32_t> image_index;
	std::vector<SplitGraph> image_basis;
	image_index.reserve(1024);
	image_basis.reserve(1024);

	auto get_image_index = [&](const SplitGraph& graph) {
		auto [it, inserted] = image_index.try_emplace(
			graph,
			static_cast<std::uint32_t>(image_basis.size())
		);
		if (inserted) {
			image_basis.push_back(graph);
		}
		return it->second;
	};

	std::size_t skipped_wheel = 0;
	const double seconds = time_seconds([&]() {
		for (const GraphType& graph : graphs) {
			if (graph == wheel) {
				++skipped_wheel;
				continue;
			}

			const auto column = GCType(graph, AssumeBasisOrderTag{}).delta().data();
			const std::uint32_t nnz = static_cast<std::uint32_t>(column.size());
			write_binary_value(columns_out, nnz);
			for (const auto& be : column) {
				const std::uint32_t row = get_image_index(be.getValue());
				const auto coefficient = be.getCoefficient().value();
				write_binary_value(columns_out, row);
				fieldType::write_serialized_value(columns_out, coefficient);
			}

			++domain_dim_u64;
			entries_u64 += nnz;
			if (domain_dim_u64 % 100000 == 0) {
				std::cout << "split map columns written = " << domain_dim_u64
				          << ", image basis = " << image_basis.size()
				          << ", entries = " << entries_u64 << '\n';
			}
		}
	});

	image_dim_u64 = static_cast<std::uint64_t>(image_basis.size());
	columns_out.seekp(domain_dim_pos);
	write_binary_value(columns_out, domain_dim_u64);
	columns_out.seekp(image_dim_pos);
	write_binary_value(columns_out, image_dim_u64);
	columns_out.seekp(entries_pos);
	write_binary_value(columns_out, entries_u64);
	columns_out.close();

	save_graph_records(image_basis_path.string(), image_basis);

	std::ofstream metadata_out(metadata_path, std::ios::trunc);
	if (!metadata_out) {
		std::cerr << "failed to write " << metadata_path.string() << '\n';
		return EXIT_FAILURE;
	}
	metadata_out << "format: GC split map column binary v1\n";
	metadata_out << "wheel: W" << +N << "\n";
	metadata_out << "rounds: " << rounds << "\n";
	metadata_out << "columns_file: " << columns_path.string() << "\n";
	metadata_out << "image_basis_file: " << image_basis_path.string() << "\n";
	metadata_out << "domain_graphs_generated: " << graphs.size() << "\n";
	metadata_out << "domain_columns_excluding_wheel: " << domain_dim_u64 << "\n";
	metadata_out << "skipped_wheel_columns: " << skipped_wheel << "\n";
	metadata_out << "image_basis_size: " << image_dim_u64 << "\n";
	metadata_out << "matrix_entries: " << entries_u64 << "\n";
	metadata_out << "field_type: " << fieldType::name() << "\n";
	metadata_out << "coefficient_storage_size_hint: " << coeff_size_u32 << "\n";
	metadata_out << "coefficient_storage_type: "
	             << (coeff_size_u32 == 0 ? "structured" : "fixed-size binary") << "\n";
	metadata_out << "domain_graph_size: (" << +GraphType::N_VERTICES_
	             << "," << +GraphType::N_EDGES_ << ")\n";
	metadata_out << "image_graph_size: (" << +SplitGraph::N_VERTICES_
	             << "," << +SplitGraph::N_EDGES_ << ")\n";

	std::cout << "wheel = W" << +N << '\n';
	std::cout << "split-contract rounds = " << rounds << '\n';
	std::cout << "domain graphs generated = " << graphs.size() << '\n';
	std::cout << "split map columns = " << domain_dim_u64 << '\n';
	std::cout << "split map image basis = " << image_dim_u64 << '\n';
	std::cout << "split map entries = " << entries_u64 << '\n';
	std::cout << "split map export time = " << seconds << " s\n";
	std::cout << "saved columns = " << columns_path.string() << '\n';
	std::cout << "saved image basis = " << image_basis_path.string() << '\n';
	std::cout << "saved metadata = " << metadata_path.string() << '\n';

	return EXIT_SUCCESS;
}

template <Int N>
int run_split_map_info_case(
	int rounds,
	const std::filesystem::path* map_output_prefix
) {
	const std::filesystem::path prefix = map_output_prefix != nullptr
		? *map_output_prefix
		: split_map_prefix<N>(rounds);
	const std::filesystem::path columns_path = prefix.string() + "_columns.bin";

	SplitMapHeader header;
	auto matrix = load_split_map_columns_as_compressed_matrix(columns_path, &header);
	if (!matrix.has_value()) {
		return EXIT_FAILURE;
	}

	std::cout << "loaded split map = " << columns_path.string() << '\n';
	std::cout << "wheel = W" << header.wheel_n << '\n';
	std::cout << "rounds = " << header.rounds << '\n';
	std::cout << "domain columns = " << matrix->domain_dim() << '\n';
	std::cout << "image rows = " << matrix->image_dim() << '\n';
	std::cout << "entries = " << matrix->rows_and_coeffs_.size() << '\n';
	std::cout << "header domain columns = " << header.domain_dim << '\n';
	std::cout << "header image rows = " << header.image_dim << '\n';
	std::cout << "header entries = " << header.entries << '\n';

	return EXIT_SUCCESS;
}

template <Int N>
int run_split_map_solve_case(
	int rounds,
	const std::filesystem::path* map_output_prefix,
	std::size_t checkpoint_interval
) {
	using GCType = OddGCdegZero<N + 1>;
	using GraphType = typename GCType::GraphType;
	using SplitGraph = typename GCType::SplitGraphType;

	const std::filesystem::path prefix = map_output_prefix != nullptr
		? *map_output_prefix
		: split_map_prefix<N>(rounds);
	const std::filesystem::path columns_path = prefix.string() + "_columns.bin";
	const std::filesystem::path image_basis_path = prefix.string() + "_image_basis.bin";
	const std::filesystem::path solution_path = prefix.string() + "_solution.txt";
	const std::filesystem::path checkpoint_path = prefix.string() + "_wiedemann_checkpoint.bin";
	const std::filesystem::path* checkpoint_path_ptr =
		checkpoint_interval > 0 ? &checkpoint_path : nullptr;

	SplitMapHeader header;
	auto matrix = load_split_map_columns_as_compressed_matrix(columns_path, &header);
	if (!matrix.has_value()) {
		return EXIT_FAILURE;
	}

	std::vector<SplitGraph> image_basis;
	if (!load_graph_records(image_basis_path.string(), image_basis)) {
		std::cerr << "failed to load image basis = "
		          << image_basis_path.string() << '\n';
		return EXIT_FAILURE;
	}
	if (image_basis.size() != static_cast<std::size_t>(header.image_dim)) {
		std::cerr << "image basis size mismatch: file has "
		          << image_basis.size() << ", header has "
		          << header.image_dim << '\n';
		return EXIT_FAILURE;
	}

	std::unordered_map<SplitGraph, std::uint32_t> image_index;
	image_index.reserve(image_basis.size());
	for (std::uint32_t i = 0; i < image_basis.size(); ++i) {
		image_index.emplace(image_basis[i], i);
	}

	auto target = std::make_unique<fieldType[]>(
		static_cast<std::size_t>(matrix->image_dim())
	);
	for (std::size_t i = 0; i < matrix->image_dim(); ++i) {
		target[i] = fieldType{};
	}

	typename GraphType::Basis wheel_basis(wheel_graph<N>(), fieldType{1});
	GraphType::std(wheel_basis);
	const auto split_target = GCType(wheel_basis, AssumeBasisOrderTag{}).delta().data();

	std::size_t missing_target_terms = 0;
	for (const auto& be : split_target) {
		const auto it = image_index.find(be.getValue());
		if (it == image_index.end()) {
			++missing_target_terms;
			continue;
		}
		target[it->second] += be.getCoefficient();
	}
	if (missing_target_terms != 0) {
		std::cout << "target terms missing from image basis = "
		          << missing_target_terms << '\n';
		return EXIT_FAILURE;
	}

	std::cout << "loaded split map = " << columns_path.string() << '\n';
	std::cout << "loaded image basis = " << image_basis_path.string() << '\n';
	std::cout << "wheel = W" << +N << '\n';
	std::cout << "rounds = " << rounds << '\n';
	std::cout << "target delta(W) terms = " << split_target.size() << '\n';
	std::cout << "domain columns = " << matrix->domain_dim() << '\n';
	std::cout << "image rows = " << matrix->image_dim() << '\n';
	std::cout << "entries = " << matrix->rows_and_coeffs_.size() << '\n';

	compressed_wiedemann_solver solver(std::move(*matrix));
	std::optional<compressed_wiedemann_solver::DenseDomainVec> solution;
	const double solve_seconds = time_seconds([&]() {
		solution = solver.solve_MX_equals_y(
			target,
			checkpoint_path_ptr,
			checkpoint_interval
		);
	});
	std::cout << "solve time = " << solve_seconds << " s\n";

	if (!solution.has_value()) {
		std::cout << "solution found = no\n";
		return EXIT_FAILURE;
	}

	std::size_t nonzero_coefficients = 0;
	for (std::size_t i = 0; i < solver.matrix().domain_dim(); ++i) {
		if ((*solution)[i] != fieldType{}) {
			++nonzero_coefficients;
		}
	}
	std::cout << "solution found = yes\n";
	std::cout << "nonzero solution coefficients = "
	          << nonzero_coefficients << '\n';

	if (solution_path.has_parent_path()) {
		std::filesystem::create_directories(solution_path.parent_path());
	}
	std::ofstream out(solution_path, std::ios::trunc);
	if (!out) {
		std::cerr << "failed to write " << solution_path.string() << '\n';
		return EXIT_FAILURE;
	}
	out << "wheel: W" << +N << "\n";
	out << "rounds: " << rounds << "\n";
	out << "columns_file: " << columns_path.string() << "\n";
	out << "domain_columns_excluding_wheel: " << solver.matrix().domain_dim() << "\n";
	out << "nonzero_solution_coefficients: " << nonzero_coefficients << "\n";
	for (std::size_t i = 0; i < solver.matrix().domain_dim(); ++i) {
		if ((*solution)[i] != fieldType{}) {
			out << i << " " << (*solution)[i] << "\n";
		}
	}
	std::cout << "saved solution = " << solution_path.string() << '\n';

	return EXIT_SUCCESS;
}

template <Int N>
int run_split_map_dual_case(
	int rounds,
	const std::filesystem::path* representative_output_path,
	const std::filesystem::path* map_output_prefix,
	std::size_t checkpoint_interval
) {
	using GCType = OddGCdegZero<N + 1>;
	using GraphType = typename GCType::GraphType;
	using SplitGraph = typename GCType::SplitGraphType;

	const std::filesystem::path prefix = map_output_prefix != nullptr
		? *map_output_prefix
		: split_map_prefix<N>(rounds);
	const std::filesystem::path columns_path = prefix.string() + "_columns.bin";
	const std::filesystem::path image_basis_path = prefix.string() + "_image_basis.bin";
	const std::filesystem::path checkpoint_path = prefix.string() + "_dual_wiedemann_checkpoint.bin";
	const std::filesystem::path* checkpoint_path_ptr =
		checkpoint_interval > 0 ? &checkpoint_path : nullptr;

	SplitMapHeader header;
	auto split_matrix = load_split_map_columns_as_compressed_matrix(columns_path, &header);
	if (!split_matrix.has_value()) {
		return EXIT_FAILURE;
	}

	std::vector<GraphType> all_domain_graphs =
		generate_same_degree_candidates<N>(rounds, nullptr);

	typename GraphType::Basis wheel_basis(wheel_graph<N>(), fieldType{1});
	GraphType::std(wheel_basis);
	const GraphType wheel = wheel_basis.getValue();

	std::vector<GraphType> domain_graphs;
	domain_graphs.reserve(all_domain_graphs.size());
	for (const auto& graph : all_domain_graphs) {
		if (!(graph == wheel)) {
			domain_graphs.push_back(graph);
		}
	}
	if (domain_graphs.size() != static_cast<std::size_t>(header.domain_dim)) {
		if (all_domain_graphs.size() == static_cast<std::size_t>(header.domain_dim) + 1) {
			domain_graphs.assign(all_domain_graphs.begin() + 1, all_domain_graphs.end());
		}
	}
	if (domain_graphs.size() != static_cast<std::size_t>(header.domain_dim)) {
		std::cerr << "domain graph count mismatch: filtered cache has "
		          << domain_graphs.size() << ", header has "
		          << header.domain_dim << '\n';
		return EXIT_FAILURE;
	}

	std::vector<SplitGraph> image_basis;
	if (!load_graph_records(image_basis_path.string(), image_basis)) {
		std::cerr << "failed to load image basis = "
		          << image_basis_path.string() << '\n';
		return EXIT_FAILURE;
	}
	if (image_basis.size() != static_cast<std::size_t>(header.image_dim)) {
		std::cerr << "image basis size mismatch: file has "
		          << image_basis.size() << ", header has "
		          << header.image_dim << '\n';
		return EXIT_FAILURE;
	}

	std::unordered_map<SplitGraph, std::uint32_t> image_index;
	image_index.reserve(image_basis.size());
	for (std::uint32_t i = 0; i < image_basis.size(); ++i) {
		image_index.emplace(image_basis[i], i);
	}

	const auto split_target = GCType(wheel, AssumeBasisOrderTag{}).delta().data();
	std::size_t missing_target_terms = 0;
	auto dense_target = split_matrix->reserve_dense_image_vec();
	for (std::size_t i = 0; i < split_matrix->image_dim(); ++i) {
		dense_target[i] = fieldType{};
	}
	for (const auto& be : split_target) {
		const auto it = image_index.find(be.getValue());
		if (it == image_index.end()) {
			++missing_target_terms;
			continue;
		}
		dense_target[it->second] += be.getCoefficient();
	}
	if (missing_target_terms != 0) {
		std::cout << "target terms missing from image basis = "
		          << missing_target_terms << '\n';
		return EXIT_FAILURE;
	}

	std::cout << "loaded split map = " << columns_path.string() << '\n';
	std::cout << "loaded image basis = " << image_basis_path.string() << '\n';
	std::cout << "wheel = W" << +N << '\n';
	std::cout << "rounds = " << rounds << '\n';
	std::cout << "domain graphs = " << domain_graphs.size() << '\n';
	std::cout << "split image basis = " << image_basis.size() << '\n';
	std::cout << "split map entries = " << split_matrix->rows_and_coeffs_.size() << '\n';
	std::cout << "target d_split(W) terms = " << split_target.size() << '\n';
	std::cout << "target terms in image basis = " << split_target.size() - missing_target_terms << '\n';

	compressed_wiedemann_solver solver(std::move(*split_matrix));
	std::optional<compressed_wiedemann_solver::DenseDomainVec> correction_coefficients;
	const double solve_seconds = time_seconds([&]() {
		correction_coefficients = solver.solve_MX_equals_y(
			dense_target,
			checkpoint_path_ptr,
			checkpoint_interval
		);
	});
	std::cout << "solve time = " << solve_seconds << " s\n";

	if (!correction_coefficients.has_value()) {
		std::cout << "dual correction found = no\n";
		return EXIT_FAILURE;
	}

	GCType graph_correction;
	std::size_t nonzero_coefficients = 0;
	for (std::size_t i = 0; i < domain_graphs.size(); ++i) {
		const fieldType coefficient = (*correction_coefficients)[i];
		if (coefficient == fieldType{}) {
			continue;
		}
		++nonzero_coefficients;
		GCType scaled(domain_graphs[i], AssumeBasisOrderTag{});
		scaled.scalar_multiply(coefficient);
		graph_correction += scaled;
	}
	graph_correction.standardize_all();
	graph_correction.sort_elements();

	GCType representative(wheel, AssumeBasisOrderTag{});
	GCType negative_correction = graph_correction;
	negative_correction.scalar_multiply(fieldType{-1});
	representative += negative_correction;
	representative.standardize_all();
	representative.sort_elements();

	const auto correction_split = graph_correction.delta();
	const auto correction_contract = contraction_of(graph_correction);
	const auto representative_contract = contraction_of(representative);

	std::cout << "dual correction found = yes\n";
	std::cout << "graph correction coefficients = " << nonzero_coefficients << '\n';
	std::cout << "graph correction terms = " << graph_correction.data().size() << '\n';
	std::cout << "representative terms = " << representative.data().size() << '\n';
	std::cout << "d_split(graph correction) terms = "
	          << correction_split.data().size() << '\n';
	std::cout << "d_contraction(graph correction) terms = "
	          << correction_contract.data().size() << '\n';
	std::cout << "d_contraction(representative) terms = "
	          << representative_contract.data().size() << '\n';
	std::cout << "representative contraction-closed = "
	          << (representative_contract.data().empty() ? "yes" : "no") << '\n';

	if (representative_output_path != nullptr) {
		if (!write_class_file(*representative_output_path, representative)) {
			std::cerr << "failed to write " << representative_output_path->string() << '\n';
			return EXIT_FAILURE;
		}
		std::cout << "saved representative = "
		          << representative_output_path->string()
		          << '\n';
	}

	return representative_contract.data().empty() ? EXIT_SUCCESS : EXIT_FAILURE;
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
	lil_matrix<fieldType> matrix;
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
		matrix.sort_cols();
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

	auto compressed_matrix = matrix.to_compressed_sparse_matrix();
	std::cout << "created sparse solver: domain_dim = " << compressed_matrix.domain_dim() << '\n';
	std::cout << "image_space_dim = " << compressed_matrix.image_dim() << '\n';
	std::cout << "number of matrix entries =" << compressed_matrix.rows_and_coeffs_.size() << '\n';

	auto dense_target = compressed_matrix.reserve_dense_image_vec();
	for (std::size_t i = 0; i < compressed_matrix.image_dim(); ++i) {
		dense_target[i] = fieldType{};
	}
	for (const auto& be : target) {
		dense_target[be.getValue()] += be.getCoefficient();
	}

	compressed_wiedemann_solver solver(std::move(compressed_matrix));
	std::optional<compressed_wiedemann_solver::DenseDomainVec> correction_coefficients;
	const double solve_seconds = time_seconds([&]() {
		correction_coefficients = solver.solve_MX_equals_y(dense_target);
	});
	std::cout << "solve time = " << solve_seconds << " s\n";

	if (!correction_coefficients.has_value()) {
		std::cout << "compressed Wiedemann returned no solution\n";
		return std::nullopt;
	}

	GCType correction;
	std::size_t nonzero_correction_coefficients = 0;
	for (std::size_t i = 0; i < matrix.domain_dim(); ++i) {
		const fieldType coefficient = (*correction_coefficients)[i];
		if (coefficient == fieldType{}) {
			continue;
		}
		++nonzero_correction_coefficients;
		GCType scaled(domain_graphs[i], AssumeBasisOrderTag{});
		scaled.scalar_multiply(coefficient);
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
	          << nonzero_correction_coefficients << '\n';
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
	graphs = generate_same_degree_candidates<N>(rounds, representatives_output_path);

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
	std::cerr << "usage: " << argv0 << " wheel_N rounds [representatives_output_path] [cycle_output_path] [mode] [checkpoint_interval]\n";
	std::cerr << "  mode defaults to generated-support; use enumerate, split-map, split-map-info, split-map-solve, split-map-dual, generated-constrained, krylov, krylov-dual, krylov-constrained, or krylov-constrained-dual\n";
	std::cerr << "  wheel_N must be one of 3,5,7,9,11,13,15,17,19,21,25,27\n";
}

} // namespace

int main(int argc, char** argv) {
	if (argc < 3 || argc > 7) {
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
		if (mode_arg == "enumerate") {
			mode = ExperimentMode::Enumerate;
		} else if (mode_arg == "split-map") {
			mode = ExperimentMode::SplitMap;
		} else if (mode_arg == "split-map-info") {
			mode = ExperimentMode::SplitMapInfo;
		} else if (mode_arg == "split-map-solve") {
			mode = ExperimentMode::SplitMapSolve;
		} else if (mode_arg == "split-map-dual") {
			mode = ExperimentMode::SplitMapDual;
		} else if (mode_arg == "krylov") {
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

	std::size_t checkpoint_interval = 0;
	if (argc >= 7) {
		const int parsed_interval = std::stoi(argv[6]);
		if (parsed_interval < 0) {
			print_usage(argv[0]);
			return EXIT_FAILURE;
		}
		checkpoint_interval = static_cast<std::size_t>(parsed_interval);
	}

#define RUN_SELECTED_CASE(N_VALUE) \
	do { \
		if (mode == ExperimentMode::Enumerate) { \
			return run_enumerate_case<N_VALUE>(rounds, representatives_output_path_ptr); \
		} \
		if (mode == ExperimentMode::SplitMap) { \
			return run_split_map_export_case<N_VALUE>(rounds, representatives_output_path_ptr, cycle_output_path_ptr); \
		} \
		if (mode == ExperimentMode::SplitMapInfo) { \
			return run_split_map_info_case<N_VALUE>(rounds, cycle_output_path_ptr); \
		} \
		if (mode == ExperimentMode::SplitMapSolve) { \
			return run_split_map_solve_case<N_VALUE>(rounds, cycle_output_path_ptr, checkpoint_interval); \
		} \
		if (mode == ExperimentMode::SplitMapDual) { \
			return run_split_map_dual_case<N_VALUE>(rounds, representatives_output_path_ptr, cycle_output_path_ptr, checkpoint_interval); \
		} \
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
