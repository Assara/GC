#pragma once

#include <array>
#include <algorithm>
#include <cstring>

template <typename GraphType>
class GraphDirections {
	public:
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
			return std::memcmp(values.data(), other.values.data(), SIZE * sizeof(bool));
		}

		bool operator<(const GraphDirections& other) const {
			return compare(other) < 0;
		}

	private:
		std::array<bool, SIZE> values{};
};
