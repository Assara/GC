#pragma once

#include <algorithm>
#include <optional>
#include <vector>

#include "OCGraph.hpp"
#include "GC.hpp"
#include "VectorSpace/BasisElement.hpp"
#include "VectorSpace/LinComb.hpp"

template <typename GraphType, Int N>
struct SplitGraphN {
	using type = typename SplitGraphN<typename GraphType::SplitGraph, N - 1>::type;
};

template <typename GraphType>
struct SplitGraphN<GraphType, 0> {
	using type = GraphType;
};

template <typename GraphType>
class OCGC {
	public:
		using ThisOCGC = OCGC<GraphType>;
		using OCGraphType = OCGraph<GraphType>;
		using SplitGraphType = typename GraphType::SplitGraph;
		using SplitOCGraphType = OCGraph<SplitGraphType>;
		using SplitOCGC = OCGC<SplitGraphType>;
		using L = VectorSpace::LinComb<OCGraphType, fieldType>;
		using Base = BasisElement<OCGraphType, fieldType>;

	private:
		L vec;

	public:
		OCGC() = default;

		explicit OCGC(const Base& basis_element)
			: vec(basis_element) {}

		explicit OCGC(const Base& basis_element, AssumeBasisOrderTag)
			: vec(basis_element, AssumeBasisOrderTag{}) {}

		explicit OCGC(const OCGraphType& graph)
			: vec(graph) {}

		explicit OCGC(const OCGraphType& graph, AssumeBasisOrderTag)
			: vec(graph, AssumeBasisOrderTag{}) {}

		explicit OCGC(L lincomb)
			: vec(std::move(lincomb)) {
			standardize_and_sort();
		}

		explicit OCGC(std::vector<Base>&& elems)
			: vec(std::move(elems), AssumeBasisOrderTag{}) {
			standardize_and_sort();
		}

		explicit OCGC(std::vector<Base>&& elems, AssumeBasisOrderTag)
			: vec(std::move(elems), AssumeBasisOrderTag{}) {}

		const L& data() const {
			return vec;
		}

		L& data_mutable() {
			return vec;
		}

		bigInt size() const {
			return vec.size();
		}

		bool empty() const {
			return vec.empty();
		}

		ThisOCGC& operator+=(const ThisOCGC& other) {
			vec += other.vec;
			return *this;
		}

		ThisOCGC operator+(const ThisOCGC& other) const {
			ThisOCGC result = *this;
			result += other;
			return result;
		}

		ThisOCGC scaled(fieldType scalar) const {
			std::vector<Base> elems;
			elems.reserve(vec.raw_elements().size());

			for (const Base& elem : vec.raw_elements()) {
				const fieldType coeff = elem.getCoefficient() * scalar;
				if (coeff != fieldType{0}) {
					elems.emplace_back(elem.getValue(), coeff);
				}
			}

			return ThisOCGC(std::move(elems));
		}

		void standardize_all() {
			auto& raw = vec.raw_elements_nonconst();
			for (Base& elem : raw) {
				elem.getValue().standardize_and_sort();
			}
			vec.remove_zeros();
		}

		void sort_elements() {
			vec.sort_elements();
		}

		void standardize_and_sort() {
			standardize_all();
			sort_elements();
		}

		OCGC delta_e() const {
			std::vector<Base> elems;
			elems.reserve(vec.raw_elements().size());

			for (const Base& elem : vec.raw_elements()) {
				OCGraphType image = elem.getValue().delta_e();
				if (image.size() == 0) {
					continue;
				}
				elems.emplace_back(std::move(image), elem.getCoefficient());
			}

			return OCGC(std::move(elems));
		}

		SplitOCGC delta_v() const {
			std::vector<typename SplitOCGC::Base> elems;

			for (const Base& elem : vec.raw_elements()) {
				const auto image = elem.getValue().delta_v();
				for (const auto& split_be : image.raw_elements()) {
					elems.emplace_back(
						split_be.getValue(),
						split_be.getCoefficient() * elem.getCoefficient()
					);
				}
			}

			return SplitOCGC(std::move(elems));
		}

		ThisOCGC delta_e_homotopy_correction() const {
			std::vector<Base> correction_terms;
			correction_terms.reserve(vec.raw_elements().size());

			for (const Base& elem : vec.raw_elements()) {
				OCGraphType correction = elem.getValue().homotopy_standardize();
				if (correction.size() != 0) {
					correction_terms.emplace_back(std::move(correction), -elem.getCoefficient());
				}
			}

			return ThisOCGC(std::move(correction_terms));
		}

		std::optional<ThisOCGC> greedy_delta_e_primitive() const {
			std::vector<Base> primitive_terms;
			primitive_terms.reserve(vec.raw_elements().size());

			for (const Base& elem : vec.raw_elements()) {
				auto primitive = elem.getValue().morse_delta_e_primitive();
				if (!primitive.has_value()) {
					return std::nullopt;
				}
				if (primitive->size() != 0) {
					primitive_terms.emplace_back(std::move(primitive.value()), elem.getCoefficient());
				}
			}

			return ThisOCGC(std::move(primitive_terms));
		}

		std::optional<ThisOCGC> exact_delta_e_primitive() const {
			std::vector<Base> primitive_terms;
			primitive_terms.reserve(vec.raw_elements().size());

			for (const Base& elem : vec.raw_elements()) {
				auto primitive = elem.getValue().exact_delta_e_primitive();
				if (!primitive.has_value()) {
					return std::nullopt;
				}
				if (primitive->size() != 0) {
					primitive_terms.emplace_back(std::move(primitive.value()), elem.getCoefficient());
				}
			}

			return ThisOCGC(std::move(primitive_terms));
		}

		std::optional<ThisOCGC> morse_delta_e_homotopy() const {
			std::vector<Base> homotopy_terms;
			homotopy_terms.reserve(vec.raw_elements().size());

			for (const Base& elem : vec.raw_elements()) {
				auto homotopy = elem.getValue().morse_delta_e_homotopy();
				if (!homotopy.has_value()) {
					return std::nullopt;
				}
				if (homotopy->size() != 0) {
					homotopy_terms.emplace_back(std::move(homotopy.value()), elem.getCoefficient());
				}
			}

			return ThisOCGC(std::move(homotopy_terms));
		}

		std::optional<ThisOCGC> morse_delta_e_primitive() const {
			std::vector<Base> primitive_terms;
			primitive_terms.reserve(vec.raw_elements().size());

			for (const Base& elem : vec.raw_elements()) {
				auto primitive = elem.getValue().morse_delta_e_primitive();
				if (!primitive.has_value()) {
					return std::nullopt;
				}
				if (primitive->size() != 0) {
					primitive_terms.emplace_back(std::move(primitive.value()), elem.getCoefficient());
				}
			}

			return ThisOCGC(std::move(primitive_terms));
		}
};

struct FullFGreedyReport {
	bool success = false;
	signedInt failed_n = -1;
	signedInt zero_n = -1;
	signedInt last_nonzero_n = -1;
	signedInt highest_valence_in_last_nonzero = 0;
};

template <typename GraphType>
Int max_vertex_valence_in_ocgc(const OCGC<GraphType>& ocgc) {
	Int best = 0;
	for (const auto& elem : ocgc.data().raw_elements()) {
		for (const Int valence : elem.getValue().graph.valence_array()) {
			best = std::max(best, valence);
		}
	}
	return best;
}

template <Int N, typename SourceGraphType>
auto OCG_Fn_greedy(const SourceGraphType& source);

template <Int N, typename SourceGraphType>
auto OCG_Fn_exact(const SourceGraphType& source);

template <Int N, typename SomeGCType>
auto OCG_Fn_greedy_on_GC(const SomeGCType& gc) {
	using SourceGraphType = typename SomeGCType::GraphType;
	using TargetGraphType = typename SplitGraphN<SourceGraphType, N>::type;
	using TargetOCGCType = OCGC<TargetGraphType>;

	std::vector<typename TargetOCGCType::Base> terms;

	for (const auto& be : gc.data().raw_elements()) {
		const auto image = OCG_Fn_greedy<N>(be.getValue());
		if (!image.has_value()) {
			return std::optional<TargetOCGCType>{};
		}
		for (const auto& image_be : image->data().raw_elements()) {
			terms.emplace_back(
				image_be.getValue(),
				image_be.getCoefficient() * be.getCoefficient()
			);
		}
	}

	return std::optional<TargetOCGCType>(TargetOCGCType(std::move(terms)));
}

template <Int N, typename SomeGCType>
auto OCG_Fn_exact_on_GC(const SomeGCType& gc) {
	using SourceGraphType = typename SomeGCType::GraphType;
	using TargetGraphType = typename SplitGraphN<SourceGraphType, N>::type;
	using TargetOCGCType = OCGC<TargetGraphType>;

	std::vector<typename TargetOCGCType::Base> terms;

	for (const auto& be : gc.data().raw_elements()) {
		const auto image = OCG_Fn_exact<N>(be.getValue());
		if (!image.has_value()) {
			return std::optional<TargetOCGCType>{};
		}
		for (const auto& image_be : image->data().raw_elements()) {
			terms.emplace_back(
				image_be.getValue(),
				image_be.getCoefficient() * be.getCoefficient()
			);
		}
	}

	return std::optional<TargetOCGCType>(TargetOCGCType(std::move(terms)));
}

template <Int N, typename SourceGraphType>
auto OCG_Fn_greedy(const SourceGraphType& source) {
	using TargetGraphType = typename SplitGraphN<SourceGraphType, N>::type;
	using TargetOCGCType = OCGC<TargetGraphType>;

	if constexpr (N == 0) {
		return std::optional<TargetOCGCType>(TargetOCGCType(OCGraph<SourceGraphType>::OCG_F0(source)));
	} else {
		using SourceGCType = GC<
			SourceGraphType::N_VERTICES_,
			SourceGraphType::N_EDGES_,
			SourceGraphType::N_OUT_HAIR_,
			SourceGraphType::N_IN_HAIR_,
			SourceGraphType::C_,
			SourceGraphType::D_
		>;

		SourceGCType gc(source);
		const auto lhs = OCG_Fn_greedy_on_GC<N - 1>(gc.delta());
		if (!lhs.has_value()) {
			return std::optional<TargetOCGCType>{};
		}

		const auto rhs = OCG_Fn_greedy<N - 1>(source);
		if (!rhs.has_value()) {
			return std::optional<TargetOCGCType>{};
		}

		TargetOCGCType boundary = lhs.value() + rhs->delta_v().scaled(fieldType{-1});
		boundary.standardize_and_sort();
		if (boundary.size() == 0) {
			return std::optional<TargetOCGCType>(TargetOCGCType{});
		}

		return boundary.morse_delta_e_primitive();
	}
}

template <Int N, typename SourceGraphType>
auto OCG_Fn_exact(const SourceGraphType& source) {
	using TargetGraphType = typename SplitGraphN<SourceGraphType, N>::type;
	using TargetOCGCType = OCGC<TargetGraphType>;

	if constexpr (N == 0) {
		return std::optional<TargetOCGCType>(TargetOCGCType(OCGraph<SourceGraphType>::OCG_F0(source)));
	} else {
		using SourceGCType = GC<
			SourceGraphType::N_VERTICES_,
			SourceGraphType::N_EDGES_,
			SourceGraphType::N_OUT_HAIR_,
			SourceGraphType::N_IN_HAIR_,
			SourceGraphType::C_,
			SourceGraphType::D_
		>;

		SourceGCType gc(source);
		const auto lhs = OCG_Fn_exact_on_GC<N - 1>(gc.delta());
		if (!lhs.has_value()) {
			return std::optional<TargetOCGCType>{};
		}

		const auto rhs = OCG_Fn_exact<N - 1>(source);
		if (!rhs.has_value()) {
			return std::optional<TargetOCGCType>{};
		}

		TargetOCGCType boundary = lhs.value() + rhs->delta_v().scaled(fieldType{-1});
		boundary.standardize_and_sort();
		if (boundary.size() == 0) {
			return std::optional<TargetOCGCType>(TargetOCGCType{});
		}

		return boundary.exact_delta_e_primitive();
	}
}

struct FullFExactReport {
	bool success = false;
	signedInt failed_n = -1;
	signedInt zero_n = -1;
	signedInt last_nonzero_n = -1;
	signedInt highest_valence_in_last_nonzero = 0;
};

template <Int MaxDepth, Int N, typename SourceGraphType>
void advance_full_F_exact_report(const SourceGraphType& source, FullFExactReport& report) {
	if constexpr (N > MaxDepth) {
		report.success = false;
		report.failed_n = MaxDepth + 1;
		return;
	} else {
		const auto fn = OCG_Fn_exact<N>(source);
		if (!fn.has_value()) {
			report.success = false;
			report.failed_n = N;
			return;
		}

		if (fn->size() == 0) {
			report.success = true;
			report.zero_n = N;
			return;
		}

		report.last_nonzero_n = N;
		report.highest_valence_in_last_nonzero = max_vertex_valence_in_ocgc(*fn);
		advance_full_F_exact_report<MaxDepth, N + 1>(source, report);
	}
}

template <typename SourceGraphType>
FullFExactReport build_full_F_exact_report(const SourceGraphType& source) {
	FullFExactReport report;
	advance_full_F_exact_report<SourceGraphType::N_EDGES_, 0>(source, report);
	return report;
}

template <Int MaxDepth, Int N, typename SourceGraphType>
void advance_full_F_greedy_report(const SourceGraphType& source, FullFGreedyReport& report) {
	if constexpr (N > MaxDepth) {
		report.success = false;
		report.failed_n = MaxDepth + 1;
		return;
	} else {
		const auto fn = OCG_Fn_greedy<N>(source);
		if (!fn.has_value()) {
			report.success = false;
			report.failed_n = N;
			return;
		}

		if (fn->size() == 0) {
			report.success = true;
			report.zero_n = N;
			return;
		}

		report.last_nonzero_n = N;
		report.highest_valence_in_last_nonzero = max_vertex_valence_in_ocgc(*fn);
		advance_full_F_greedy_report<MaxDepth, N + 1>(source, report);
	}
}

template <typename SourceGraphType>
FullFGreedyReport build_full_F_greedy_report(const SourceGraphType& source) {
	FullFGreedyReport report;
	advance_full_F_greedy_report<SourceGraphType::N_EDGES_, 0>(source, report);
	return report;
}

template <Int N, typename SourceGraphType>
auto OCG_Fn(const SourceGraphType& source);

template <Int N, typename SourceGraphType>
auto OCG_Fn_on_GC(
	const GC<
		SourceGraphType::N_VERTICES_,
		SourceGraphType::N_EDGES_,
		SourceGraphType::N_OUT_HAIR_,
		SourceGraphType::N_IN_HAIR_,
		SourceGraphType::C_,
		SourceGraphType::D_
	>& gc
) {
	using TargetGraphType = typename SplitGraphN<SourceGraphType, N>::type;
	using TargetOCGCType = OCGC<TargetGraphType>;

	std::vector<typename TargetOCGCType::Base> terms;

	for (const auto& be : gc.data().raw_elements()) {
		auto image = OCG_Fn<N>(be.getValue());
		for (const auto& image_be : image.data().raw_elements()) {
			terms.emplace_back(
				image_be.getValue(),
				image_be.getCoefficient() * be.getCoefficient()
			);
		}
	}

	return TargetOCGCType(std::move(terms));
}

template <Int N, typename SourceGraphType>
auto OCG_Fn(const SourceGraphType& source) {
	if constexpr (N == 0) {
		return OCGC<SourceGraphType>(OCGraph<SourceGraphType>::OCG_F0(source));
	} else {
		using TargetGraphType = typename SplitGraphN<SourceGraphType, N>::type;
		using TargetOCGCType = OCGC<TargetGraphType>;
		using SourceGCType = GC<
			SourceGraphType::N_VERTICES_,
			SourceGraphType::N_EDGES_,
			SourceGraphType::N_OUT_HAIR_,
			SourceGraphType::N_IN_HAIR_,
			SourceGraphType::C_,
			SourceGraphType::D_
		>;

		SourceGCType gc(source);
		auto lhs = OCG_Fn_on_GC<N - 1>(gc.delta());
		auto rhs = OCG_Fn<N - 1>(source).delta_v();
		TargetOCGCType boundary = lhs + rhs.scaled(fieldType{-1});
		return boundary.delta_e_homotopy_correction();
	}
}
