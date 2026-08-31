#pragma once

#include <bit>
#include <concepts>
#include <cstddef>
#include <cstdlib>
#include <limits>
#include <memory>
#include <new>
#include <type_traits>
#include <utility>

template <typename Value>
concept linear_probe_value = std::default_initializable<Value>
	&& std::equality_comparable<Value>
	&& std::is_nothrow_default_constructible_v<Value>
	&& std::is_nothrow_copy_constructible_v<Value>
	&& std::is_nothrow_move_constructible_v<Value>
	&& std::is_nothrow_copy_assignable_v<Value>
	&& std::is_nothrow_move_assignable_v<Value>
	&& std::is_nothrow_destructible_v<Value>
	&& requires(const Value& value) {
		{ value.hash() } noexcept -> std::convertible_to<std::size_t>;
		{ value.empty() } noexcept -> std::same_as<bool>;
		{ value == value } noexcept -> std::same_as<bool>;
	};

// Append-oriented open-addressing set with direct Value storage.
template <linear_probe_value Value>
class linear_probe_set {
	public:
		using value_type = Value;
		using size_type = std::size_t;

		static constexpr size_type MINIMUM_CAPACITY = 16;
		static constexpr size_type LOAD_NUMERATOR = 7;
		static constexpr size_type LOAD_DENOMINATOR = 8;

		linear_probe_set() = default;
		linear_probe_set(const linear_probe_set&) = delete;
		linear_probe_set& operator=(const linear_probe_set&) = delete;
		linear_probe_set(linear_probe_set&& other) noexcept
			: values_(std::move(other.values_)),
			  size_(std::exchange(other.size_, 0)),
			  capacity_(std::exchange(other.capacity_, 0)) {}

		linear_probe_set& operator=(linear_probe_set&& other) noexcept {
			if (this == &other) return *this;
			values_ = std::move(other.values_);
			size_ = std::exchange(other.size_, 0);
			capacity_ = std::exchange(other.capacity_, 0);
			return *this;
		}

		bool insert(const Value& value) noexcept {
			if (value.empty()) return false;

			const size_type old_capacity = capacity_;
			size_type insertion_index = 0;
			if (old_capacity != 0) {
				const size_type mask = old_capacity - 1;
				insertion_index = value.hash() & mask;
				while (!values_[insertion_index].empty()) {
					if (values_[insertion_index] == value) return false;
					insertion_index = (insertion_index + 1) & mask;
				}
			}

			const size_type new_capacity = capacity_for_size(size_ + 1);
			if (new_capacity != old_capacity) {
				reallocate(new_capacity, old_capacity);
				capacity_ = new_capacity;
				const size_type mask = new_capacity - 1;
				insertion_index = value.hash() & mask;
				while (!values_[insertion_index].empty()) {
					insertion_index = (insertion_index + 1) & mask;
				}
			}
			values_[insertion_index] = value;
			++size_;
			return true;
		}

		void reserve(size_type expected_elements) noexcept {
			const size_type old_capacity = capacity_;
			const size_type new_capacity = capacity_for_size(expected_elements);
			if (new_capacity <= old_capacity) return;
			reallocate(new_capacity, old_capacity);
			capacity_ = new_capacity;
		}

		bool contains(const Value& value) const noexcept {
			if (value.empty() || capacity_ == 0) return false;
			const size_type mask = capacity_ - 1;
			size_type index = value.hash() & mask;
			while (!values_[index].empty()) {
				if (values_[index] == value) return true;
				index = (index + 1) & mask;
			}
			return false;
		}

		size_type size() const noexcept { return size_; }

		size_type capacity() const noexcept { return capacity_; }

		Value* data() noexcept { return values_.get(); }
		const Value* data() const noexcept { return values_.get(); }

	private:
		static size_type capacity_for_size(size_type element_count) noexcept {
			if (element_count == 0) return 0;
			if (element_count > std::numeric_limits<size_type>::max()
				/ LOAD_DENOMINATOR) {
				std::abort();
			}
			const size_type required
				= (element_count * LOAD_DENOMINATOR + LOAD_NUMERATOR - 1)
				/ LOAD_NUMERATOR;
			const size_type requested
				= required < MINIMUM_CAPACITY ? MINIMUM_CAPACITY : required;
			const size_type maximum_size_capacity
				= std::bit_floor(std::numeric_limits<size_type>::max());
			const size_type maximum_array_capacity = std::bit_floor(
				std::numeric_limits<size_type>::max() / sizeof(Value)
			);
			const size_type maximum_capacity
				= maximum_size_capacity < maximum_array_capacity
				? maximum_size_capacity
				: maximum_array_capacity;
			if (requested > maximum_capacity) {
				std::abort();
			}
			return std::bit_ceil(requested);
		}

		void reallocate(size_type new_capacity, size_type old_capacity) noexcept {
			std::unique_ptr<Value[]> old_values = std::move(values_);
			if (new_capacity == 0) return;
			values_.reset(new (std::nothrow) Value[new_capacity]);
			if (values_ == nullptr) std::abort();
			if (!values_[0].empty()) {
				std::abort();
			}
			for (size_type index = 0; index < old_capacity; ++index) {
				if (old_values[index].empty()) continue;
				const size_type mask = new_capacity - 1;
				size_type new_index = old_values[index].hash() & mask;
				while (!values_[new_index].empty()) {
					new_index = (new_index + 1) & mask;
				}
				values_[new_index] = std::move(old_values[index]);
			}
		}

		std::unique_ptr<Value[]> values_;
		size_type size_ = 0;
		size_type capacity_ = 0;
};
