#pragma once
#include <cstdint>
#include <istream>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>

// Fixed-prime finite field with 32-bit storage.
// Prime p = 4294967291, the largest 32-bit prime.
// Multiplication uses a 64-bit intermediate exactly.

class Z4294967291 {
	static constexpr std::uint32_t MOD = 4294967291u;
	std::uint32_t v_;

	static Z4294967291 pow_mod(Z4294967291 base, std::uint64_t exponent) {
		Z4294967291 result{1};
		while (exponent != 0) {
			if ((exponent & 1u) != 0) {
				result *= base;
			}
			exponent >>= 1u;
			if (exponent != 0) {
				base *= base;
			}
		}
		return result;
	}

public:
	using serialized_value_type = std::uint32_t;

	constexpr Z4294967291() : v_(0) {}

	template <class T, std::enable_if_t<std::is_integral_v<T>, int> = 0>
	/* not explicit */ constexpr Z4294967291(T x) {
		long long r = static_cast<long long>(x) % static_cast<long long>(MOD);
		if (r < 0) {
			r += static_cast<long long>(MOD);
		}
		v_ = static_cast<std::uint32_t>(r);
	}

	explicit constexpr Z4294967291(serialized_value_type x) : v_(x % MOD) {}

	serialized_value_type value() const { return v_; }
	static constexpr std::uint32_t characteristic() { return MOD; }
	static constexpr std::uint32_t serialized_value_size_hint() {
		return static_cast<std::uint32_t>(sizeof(serialized_value_type));
	}
	static std::string name() { return "Z4294967291"; }
	static Z4294967291 sample(std::mt19937_64& rng) {
		std::uniform_int_distribution<std::uint64_t> dist(0, characteristic() - 1);
		return Z4294967291{dist(rng)};
	}

	explicit operator bool() const { return v_ != 0; }

	friend bool operator==(const Z4294967291& a, const Z4294967291& b) { return a.v_ == b.v_; }
	friend bool operator!=(const Z4294967291& a, const Z4294967291& b) { return a.v_ != b.v_; }

	Z4294967291 operator-() const {
		if (v_ == 0) {
			return Z4294967291{0};
		}
		Z4294967291 out;
		out.v_ = MOD - v_;
		return out;
	}

	Z4294967291 operator+(const Z4294967291& o) const {
		std::uint64_t s = static_cast<std::uint64_t>(v_) + static_cast<std::uint64_t>(o.v_);
		if (s >= MOD) {
			s -= MOD;
		}
		Z4294967291 out;
		out.v_ = static_cast<std::uint32_t>(s);
		return out;
	}

	Z4294967291& operator+=(const Z4294967291& o) {
		std::uint64_t s = static_cast<std::uint64_t>(v_) + static_cast<std::uint64_t>(o.v_);
		if (s >= MOD) {
			s -= MOD;
		}
		v_ = static_cast<std::uint32_t>(s);
		return *this;
	}

	Z4294967291 operator-(const Z4294967291& o) const {
		std::uint64_t d = static_cast<std::uint64_t>(v_) - static_cast<std::uint64_t>(o.v_);
		if (v_ < o.v_) {
			d += MOD;
		}
		Z4294967291 out;
		out.v_ = static_cast<std::uint32_t>(d);
		return out;
	}

	Z4294967291& operator-=(const Z4294967291& o) {
		if (v_ < o.v_) {
			v_ = static_cast<std::uint32_t>(
				static_cast<std::uint64_t>(v_) + MOD - static_cast<std::uint64_t>(o.v_)
			);
		} else {
			v_ = static_cast<std::uint32_t>(v_ - o.v_);
		}
		return *this;
	}

	Z4294967291 operator*(const Z4294967291& o) const {
		const std::uint64_t p = static_cast<std::uint64_t>(v_) * static_cast<std::uint64_t>(o.v_);
		Z4294967291 out;
		out.v_ = static_cast<std::uint32_t>(p % MOD);
		return out;
	}

	Z4294967291& operator*=(const Z4294967291& o) { return *this = *this * o; }

	Z4294967291 operator/(const Z4294967291& o) const {
		if (o.v_ == 0) {
			throw std::domain_error("division by zero");
		}
		return *this * o.inv();
	}

	Z4294967291& operator/=(const Z4294967291& o) {
		return *this = *this / o;
	}

	Z4294967291 inv() const {
		if (v_ == 0) {
			throw std::domain_error("division by zero");
		}
		return pow_mod(*this, MOD - 2u);
	}

	static void write_serialized_value(std::ostream& out, const serialized_value_type& value) {
		out.write(reinterpret_cast<const char*>(&value), sizeof(value));
	}

	static bool read_serialized_value(std::istream& in, serialized_value_type& value) {
		in.read(reinterpret_cast<char*>(&value), sizeof(value));
		return static_cast<bool>(in);
	}

	static void write_value(std::ostream& out, const Z4294967291& value) {
		write_serialized_value(out, value.value());
	}

	static bool read_value(std::istream& in, Z4294967291& value) {
		serialized_value_type storage{};
		if (!read_serialized_value(in, storage)) {
			return false;
		}
		value = Z4294967291{storage};
		return true;
	}

	std::int64_t signed_representative() const {
		return v_ <= MOD / 2
			? static_cast<std::int64_t>(v_)
			: static_cast<std::int64_t>(v_) - static_cast<std::int64_t>(MOD);
	}
};

inline std::ostream& operator<<(std::ostream& os, const Z4294967291& a) {
	return os << a.value();
}

inline std::istream& operator>>(std::istream& is, Z4294967291& a) {
	long long x;
	if (is >> x) {
		a = Z4294967291{x};
	}
	return is;
}
