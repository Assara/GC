#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

#include <optional>
#include <unordered_map>
#include <vector>

#include "GC.hpp"

namespace RationalCycleConverter {

	using BigInt = boost::multiprecision::cpp_int;
	using Rational = boost::rational<BigInt>;

	template <typename GraphType>
	using RationalGraph = typename GraphType::template RebindField<Rational>;

	template <typename GraphType>
	class RationalGC {
	public:
		using ThisGraph = GraphType;
		using ContGraphType = typename ThisGraph::ContGraph;
		using L = VectorSpace::LinComb<ThisGraph, Rational>;
		using ContL = VectorSpace::LinComb<ContGraphType, Rational>;
		using Base = BasisElement<ThisGraph, Rational>;
		using ContGC = RationalGC<ContGraphType>;

	private:
		L vec;

	public:
		RationalGC() = default;
		explicit RationalGC(const ThisGraph& g) : vec(g) {}
		explicit RationalGC(L lin_comb) : vec(std::move(lin_comb)) {}
		explicit RationalGC(std::vector<Base>&& elems) : vec(std::move(elems)) {}

		const L& data() const { return vec; }
		bigInt size() const { return vec.size(); }
		void standardize_all() { vec.standardize_all(); }
		void sort_elements() { vec.sort_elements(); }
		void standardize_and_sort() { vec.standardize_and_sort(); }

		ContGC d_contraction() const {
			std::vector<BasisElement<ContGraphType, Rational>> elems;
			elems.reserve(vec.raw_elements().size() * ThisGraph::N_EDGES_);

			for (const auto& be : vec.raw_elements()) {
				for (Int i = 0; i < ThisGraph::N_EDGES_; ++i) {
					auto contracted = ThisGraph::contract_edge(be, i);
					if (contracted.getCoefficient() != Rational{0}) {
						elems.push_back(std::move(contracted));
					}
				}
			}

			return ContGC(std::move(elems));
		}
	};

	template <typename SourceGC>
	using RationalGCOf = RationalGC<RationalGraph<typename SourceGC::GraphType>>;

	template <typename SourceGC>
	RationalGCOf<SourceGC> cast_gc_coefficients(const SourceGC& source) {
		using SrcGraph = typename SourceGC::GraphType;
		using DstGraph = RationalGraph<SrcGraph>;

		std::vector<BasisElement<DstGraph, Rational>> elems;
		elems.reserve(source.data().size());

		for (const auto& be : source.data()) {
			elems.emplace_back(
				be.getValue().template cast_field<Rational>(),
				Rational(be.getCoefficient().value())
			);
		}

		return RationalGCOf<SourceGC>(std::move(elems));
	}

	template <typename T>
	std::pair<std::vector<std::vector<T>>, std::vector<std::size_t>>
	rref(std::vector<std::vector<T>> matrix) {
		if (matrix.empty()) {
			return {std::move(matrix), {}};
		}

		const std::size_t m = matrix.size();
		const std::size_t n = matrix[0].size();

		std::size_t row = 0;
		std::vector<std::size_t> pivots;
		pivots.reserve(std::min(m, n));

		for (std::size_t col = 0; col < n && row < m; ++col) {
			std::size_t pivot = row;
			while (pivot < m && matrix[pivot][col] == T{0}) {
				++pivot;
			}
			if (pivot == m) {
				continue;
			}

			if (pivot != row) {
				std::swap(matrix[pivot], matrix[row]);
			}

			const T inv = T{1} / matrix[row][col];
			for (std::size_t j = col; j < n; ++j) {
				matrix[row][j] *= inv;
			}

			for (std::size_t i = 0; i < m; ++i) {
				if (i == row || matrix[i][col] == T{0}) {
					continue;
				}

				const T factor = matrix[i][col];
				for (std::size_t j = col; j < n; ++j) {
					matrix[i][j] -= factor * matrix[row][j];
				}
			}

			pivots.push_back(col);
			++row;
		}

		return {std::move(matrix), std::move(pivots)};
	}

	template <typename T>
	std::vector<std::vector<T>> nullspace_basis(std::vector<std::vector<T>> matrix) {
		if (matrix.empty()) {
			return {};
		}

		const std::size_t n = matrix[0].size();
		auto [reduced, pivots] = rref(std::move(matrix));

		std::vector<bool> is_pivot(n, false);
		for (std::size_t c : pivots) {
			is_pivot[c] = true;
		}

		std::vector<std::size_t> free_cols;
		for (std::size_t c = 0; c < n; ++c) {
			if (!is_pivot[c]) {
				free_cols.push_back(c);
			}
		}

		std::vector<std::vector<T>> basis;
		basis.reserve(free_cols.size());

		for (std::size_t free_col : free_cols) {
			std::vector<T> vec(n, T{0});
			vec[free_col] = T{1};
			for (std::size_t r = 0; r < pivots.size(); ++r) {
				vec[pivots[r]] = -reduced[r][free_col];
			}
			basis.push_back(std::move(vec));
		}

		return basis;
	}

	template <typename T>
	std::optional<std::vector<T>> find_all_nonzero_kernel_vector(
		const std::vector<std::vector<T>>& basis,
		std::size_t support_size,
		int max_t = 128
	) {
		if (basis.empty()) {
			return std::nullopt;
		}

		for (std::size_t j = 0; j < support_size; ++j) {
			bool seen_nonzero = false;
			for (const auto& vec : basis) {
				if (vec[j] != T{0}) {
					seen_nonzero = true;
					break;
				}
			}
			if (!seen_nonzero) {
				return std::nullopt;
			}
		}

		for (int t = 1; t <= max_t; ++t) {
			T t_power = T{1};
			std::vector<T> combined(support_size, T{0});

			for (const auto& vec : basis) {
				for (std::size_t j = 0; j < support_size; ++j) {
					combined[j] += t_power * vec[j];
				}
				t_power *= T{t};
			}

			bool all_nonzero = true;
			for (const auto& coeff : combined) {
				if (coeff == T{0}) {
					all_nonzero = false;
					break;
				}
			}

			if (all_nonzero) {
				return combined;
			}
		}

		return std::nullopt;
	}

	template <typename SourceGC>
	std::optional<RationalGCOf<SourceGC>> rationalize_cycle_on_same_support(const SourceGC& source) {
		using SrcGraph = typename SourceGC::GraphType;
		using RatGraph = RationalGraph<SrcGraph>;
		using RatContGraph = typename RatGraph::ContGraph;

		std::vector<RatGraph> support;
		support.reserve(source.data().size());

		for (const auto& be : source.data()) {
			support.push_back(be.getValue().template cast_field<Rational>());
		}

		std::vector<VectorSpace::LinComb<RatContGraph, Rational>> columns;
		columns.reserve(support.size());

		std::unordered_map<RatContGraph, std::size_t> row_index;
		std::vector<RatContGraph> rows;

		auto add_rows = [&](const auto& col) {
			for (const auto& be : col) {
				if (!row_index.contains(be.getValue())) {
					row_index.emplace(be.getValue(), rows.size());
					rows.push_back(be.getValue());
				}
			}
		};

		for (const auto& graph : support) {
			columns.push_back(graph.contraction_differential(Rational{1}));
			add_rows(columns.back());
		}

		std::vector<std::vector<Rational>> matrix(
			rows.size(),
			std::vector<Rational>(support.size(), Rational{0})
		);

		for (std::size_t j = 0; j < columns.size(); ++j) {
			for (const auto& be : columns[j]) {
				matrix[row_index.at(be.getValue())][j] += be.getCoefficient();
			}
		}

		auto basis = nullspace_basis(std::move(matrix));
		auto coeffs = find_all_nonzero_kernel_vector(basis, support.size());
		if (!coeffs.has_value()) {
			return std::nullopt;
		}

		std::vector<BasisElement<RatGraph, Rational>> elems;
		elems.reserve(support.size());
		for (std::size_t i = 0; i < support.size(); ++i) {
			elems.emplace_back(support[i], (*coeffs)[i]);
		}

		RationalGCOf<SourceGC> result(std::move(elems));
		result.standardize_and_sort();
		return result;
	}

	template <typename SourceGC>
	std::optional<RationalGCOf<SourceGC>> rationalize_split_cocycle_on_same_support(const SourceGC& source) {
		using SrcGraph = typename SourceGC::GraphType;
		using RatGraph = RationalGraph<SrcGraph>;
		using RatSplitGraph = typename RatGraph::SplitGraph;

		std::vector<RatGraph> support;
		support.reserve(source.data().size());

		for (const auto& be : source.data()) {
			support.push_back(be.getValue().template cast_field<Rational>());
		}

		std::vector<VectorSpace::LinComb<RatSplitGraph, Rational>> columns;
		columns.reserve(support.size());

		std::unordered_map<RatSplitGraph, std::size_t> row_index;
		std::vector<RatSplitGraph> rows;

		auto add_rows = [&](const auto& col) {
			for (const auto& be : col) {
				if (!row_index.contains(be.getValue())) {
					row_index.emplace(be.getValue(), rows.size());
					rows.push_back(be.getValue());
				}
			}
		};

		for (const auto& graph : support) {
			columns.push_back(graph.split_vertex_differential(Rational{1}));
			add_rows(columns.back());
		}

		std::vector<std::vector<Rational>> matrix(
			rows.size(),
			std::vector<Rational>(support.size(), Rational{0})
		);

		for (std::size_t j = 0; j < columns.size(); ++j) {
			for (const auto& be : columns[j]) {
				matrix[row_index.at(be.getValue())][j] += be.getCoefficient();
			}
		}

		auto basis = nullspace_basis(std::move(matrix));
		auto coeffs = find_all_nonzero_kernel_vector(basis, support.size());
		if (!coeffs.has_value()) {
			return std::nullopt;
		}

		std::vector<BasisElement<RatGraph, Rational>> elems;
		elems.reserve(support.size());
		for (std::size_t i = 0; i < support.size(); ++i) {
			elems.emplace_back(support[i], (*coeffs)[i]);
		}

		RationalGCOf<SourceGC> result(std::move(elems));
		result.standardize_and_sort();
		return result;
	}

} // namespace RationalCycleConverter
