#pragma once
#include <algorithm>
#include <cstddef>
#include <memory>
#include <span>
#include <utility>

namespace VectorSpace {
// Fixed-size, move-only ownership. Moving also clears the source length.
template<class T>
class OwnedArray {
    std::unique_ptr<T[]> storage_;
    std::size_t size_ = 0;
public:
    OwnedArray() = default;
    explicit OwnedArray(std::size_t n)
        : storage_(n ? std::make_unique<T[]>(n) : nullptr), size_(n) {}
    explicit OwnedArray(std::span<const T> values)
        : storage_(values.empty() ? nullptr : std::make_unique_for_overwrite<T[]>(values.size())),
          size_(values.size()) {
        std::copy(values.begin(), values.end(), begin());
    }
    OwnedArray(const OwnedArray&) = delete;
    OwnedArray& operator=(const OwnedArray&) = delete;
    OwnedArray(OwnedArray&& other) noexcept
        : storage_(std::move(other.storage_)), size_(std::exchange(other.size_, 0)) {}
    OwnedArray& operator=(OwnedArray&& other) noexcept {
        if (this != &other) {
            storage_ = std::move(other.storage_);
            size_ = std::exchange(other.size_, 0);
        }
        return *this;
    }
    std::size_t size() const noexcept { return size_; }
    bool empty() const noexcept { return size_ == 0; }
    T* data() noexcept { return storage_.get(); }
    const T* data() const noexcept { return storage_.get(); }
    T* begin() noexcept { return data(); }
    const T* begin() const noexcept { return data(); }
    T* end() noexcept { return size_ ? data() + size_ : data(); }
    const T* end() const noexcept { return size_ ? data() + size_ : data(); }
    T& operator[](std::size_t i) noexcept { return data()[i]; }
    const T& operator[](std::size_t i) const noexcept { return data()[i]; }
    T& back() noexcept { return data()[size_ - 1]; }
    const T& back() const noexcept { return data()[size_ - 1]; }
};
}
