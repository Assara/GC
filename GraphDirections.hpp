#pragma once

#include <array>
#include <algorithm>
#include <cstring>
#include <vector>

#include "VectorSpace/BasisElement.hpp"

template <typename GraphType>
class GraphDirections {
	public:
		using ThisType = GraphDirections<GraphType>;
		static constexpr Int SIZE = GraphType::SIZE;

		GraphDirections() = default;

		explicit GraphDirections(const std::array<bool, SIZE>& input)
			: values(input) {}

		constexpr bigInt size() const {
			return SIZE;
		}

		bool operator[](Int i) const {
			return values[i];
		}

		bool& operator[](Int i) {
			return values[i];
		}

		void fill(bool value) {
			values.fill(value);
		}

		const std::array<bool, SIZE>& raw() const {
			return values;
		}

		std::array<bool, SIZE>& raw() {
			return values;
		}

		bool operator==(const GraphDirections& other) const {
			return values == other.values;
		}

		bool operator!=(const GraphDirections& other) const {
			return !(*this == other);
		}

		signedInt compare(const GraphDirections& other) const {
			return -std::memcmp(values.data(), other.values.data(), SIZE * sizeof(bool));
		}

		bool operator<(const GraphDirections& other) const {
			return compare(other) < 0;
		}

		template <typename FieldType, typename IsoType>
		static BasisElement<ThisType, FieldType> standardize(
			const BasisElement<ThisType, FieldType>& input,
			const std::vector<IsoType>& automorphisms
		) {
			if (automorphisms.empty()) {
				return input;
			}

			ThisType best_value = automorphisms.front().permute(input.getValue());
			signedInt best_sign = automorphisms.front().direction_sign(input.getValue());
			bool contains_plus = best_sign > 0;
			bool contains_minus = best_sign < 0;

			for (bigInt i = 1; i < automorphisms.size(); ++i) {
				const ThisType candidate = automorphisms[i].permute(input.getValue());
				const signedInt candidate_sign = automorphisms[i].direction_sign(input.getValue());
				const signedInt comparison = candidate.compare(best_value);

				if (comparison < 0) {
					best_value = candidate;
					contains_plus = candidate_sign > 0;
					contains_minus = candidate_sign < 0;
					continue;
				}

				if (comparison == 0) {
					contains_plus = contains_plus || (candidate_sign > 0);
					contains_minus = contains_minus || (candidate_sign < 0);
				}
			}

			if (contains_plus && contains_minus) {
				return BasisElement<ThisType, FieldType>(best_value, FieldType{0});
			}

			const FieldType coeff = contains_plus ? input.getCoefficient() : -input.getCoefficient();
			return BasisElement<ThisType, FieldType>(best_value, coeff);
		}

	private:
		std::array<bool, SIZE> values{};
};
