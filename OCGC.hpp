#pragma once

#include <vector>

#include "OCGraph.hpp"
#include "VectorSpace/BasisElement.hpp"
#include "VectorSpace/LinComb.hpp"

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
				const OCGraphType& ocgraph = elem.getValue();
				const auto graph_splits = ocgraph.graph.split_vertex_differential(elem.getCoefficient());

				for (const auto& split_be : graph_splits.raw_elements()) {
					typename SplitOCGraphType::Base split_directions;

					for (const auto& direction_be : ocgraph.raw_elements()) {
						typename SplitOCGraphType::DirectionType left_right;
						left_right.fill(false);
						for (Int i = 0; i < GraphType::SIZE; ++i) {
							left_right[i] = direction_be.getValue()[i];
						}
						left_right[SplitGraphType::SIZE - 1] = true;
						split_directions.append_in_basis_order(left_right, direction_be.getCoefficient());

						typename SplitOCGraphType::DirectionType right_left;
						right_left.fill(false);
						for (Int i = 0; i < GraphType::SIZE; ++i) {
							right_left[i] = direction_be.getValue()[i];
						}
						right_left[SplitGraphType::SIZE - 2] = true;
						// TODO: the (true,false) copy should likely differ by a sign from (false,true).
						split_directions.append_in_basis_order(right_left, direction_be.getCoefficient());
					}

					SplitOCGraphType split_ocgraph(split_be.getValue(), std::move(split_directions));
					elems.emplace_back(std::move(split_ocgraph), split_be.getCoefficient());
				}
			}

			return SplitOCGC(std::move(elems));
		}
};
