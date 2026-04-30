#pragma once

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
};

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
