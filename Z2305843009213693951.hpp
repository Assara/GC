#pragma once
#include <cstdint>
#include <istream>
#include <ostream>
#include <random>
#include <stdexcept>
#include <string>
#include <type_traits>

// Fixed-prime finite field with 64-bit storage.
// Prime p = 2305843009213693951 = 2^61 - 1.
// Uses 128-bit intermediate multiplication.

class Z2305843009213693951 {
	static constexpr std::uint64_t MOD = 2305843009213693951ull;
	std::uint64_t v_;

	static Z2305843009213693951 pow_mod(Z2305843009213693951 base, std::uint64_t exponent) {
		Z2305843009213693951 result{1};
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
	using serialized_value_type = std::uint64_t;

	constexpr Z2305843009213693951() : v_(0) {}

	template <class T, std::enable_if_t<std::is_integral_v<T>, int> = 0>
	/* not explicit */ constexpr Z2305843009213693951(T x) {
		long long r = static_cast<long long>(x) % static_cast<long long>(MOD);
		if (r < 0) {
			r += static_cast<long long>(MOD);
		}
		v_ = static_cast<std::uint64_t>(r);
	}

	explicit constexpr Z2305843009213693951(serialized_value_type x) : v_(x % MOD) {}

	serialized_value_type value() const { return v_; }
	static constexpr std::uint64_t characteristic() { return MOD; }
	static constexpr std::uint32_t serialized_value_size_hint() {
		return static_cast<std::uint32_t>(sizeof(serialized_value_type));
	}
	static std::string name() { return "Z2305843009213693951"; }
	static Z2305843009213693951 sample(std::mt19937_64& rng) {
		std::uniform_int_distribution<std::uint64_t> dist(0, characteristic() - 1);
		return Z2305843009213693951{dist(rng)};
	}

	explicit operator bool() const { return v_ != 0; }

	friend bool operator==(const Z2305843009213693951& a, const Z2305843009213693951& b) { return a.v_ == b.v_; }
	friend bool operator!=(const Z2305843009213693951& a, const Z2305843009213693951& b) { return a.v_ != b.v_; }

	Z2305843009213693951 operator-() const {
		if (v_ == 0) {
			return Z2305843009213693951{0};
		}
		Z2305843009213693951 out;
		out.v_ = MOD - v_;
		return out;
	}

	Z2305843009213693951 operator+(const Z2305843009213693951& o) const {
		std::uint64_t s = v_ + o.v_;
		if (s >= MOD) {
			s -= MOD;
		}
		Z2305843009213693951 out;
		out.v_ = s;
		return out;
	}

	Z2305843009213693951& operator+=(const Z2305843009213693951& o) {
		std::uint64_t s = v_ + o.v_;
		if (s >= MOD) {
			s -= MOD;
		}
		v_ = s;
		return *this;
	}

	Z2305843009213693951 operator-(const Z2305843009213693951& o) const {
		std::uint64_t d = v_ - o.v_;
		if (v_ < o.v_) {
			d += MOD;
		}
		Z2305843009213693951 out;
		out.v_ = d;
		return out;
	}

	Z2305843009213693951& operator-=(const Z2305843009213693951& o) {
		if (v_ < o.v_) {
			v_ = v_ + MOD - o.v_;
		} else {
			v_ = v_ - o.v_;
		}
		return *this;
	}

	Z2305843009213693951 operator*(const Z2305843009213693951& o) const {
		const __uint128_t p = static_cast<__uint128_t>(v_) * static_cast<__uint128_t>(o.v_);
		Z2305843009213693951 out;
		out.v_ = static_cast<std::uint64_t>(p % MOD);
		return out;
	}

	Z2305843009213693951& operator*=(const Z2305843009213693951& o) { return *this = *this * o; }

	Z2305843009213693951 operator/(const Z2305843009213693951& o) const {
		if (o.v_ == 0) {
			throw std::domain_error("division by zero");
		}
		return *this * o.inv();
	}

	Z2305843009213693951& operator/=(const Z2305843009213693951& o) {
		return *this = *this / o;
	}

	Z2305843009213693951 inv() const {
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

	static void write_value(std::ostream& out, const Z2305843009213693951& value) {
		write_serialized_value(out, value.value());
	}

	static bool read_value(std::istream& in, Z2305843009213693951& value) {
		serialized_value_type storage{};
		if (!read_serialized_value(in, storage)) {
			return false;
		}
		value = Z2305843009213693951{storage};
		return true;
	}

	std::int64_t signed_representative() const {
		return v_ <= MOD / 2
			? static_cast<std::int64_t>(v_)
			: static_cast<std::int64_t>(v_) - static_cast<std::int64_t>(MOD);
	}
};

inline std::ostream& operator<<(std::ostream& os, const Z2305843009213693951& a) {
	return os << a.value();
}

inline std::istream& operator>>(std::istream& is, Z2305843009213693951& a) {
	long long x;
	if (is >> x) {
		a = Z2305843009213693951{x};
	}
	return is;
}
