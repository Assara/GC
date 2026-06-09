#pragma once
#include <cstdint>
#include <stdexcept>
#include <ostream>
#include <istream>
#include <random>
#include <string>
#include <type_traits>

// Fixed-prime finite field with 32-bit storage.
// Prime p = 32783 (< 2^16). Implements +, -, *, +=, -= and inverse().
// No namespace, class name: Z32783.

class Z32783 {
	static constexpr std::uint32_t MOD = 32783u; // prime
	std::uint32_t v_;                            // always in [0, MOD-1]

	public:
	using serialized_value_type = std::uint32_t;
	// Constructors
	constexpr Z32783() : v_(0) {}

	template <class T, std::enable_if_t<std::is_integral_v<T>, int> = 0>
		/* not explicit */ constexpr Z32783(T x) {
			long long r = static_cast<long long>(x) % static_cast<long long>(MOD);
			if (r < 0) r += MOD;
			v_ = static_cast<std::uint32_t>(r);
	}

	// Accessors
	serialized_value_type value() const { return v_; }
	static constexpr std::uint32_t characteristic() { return MOD; }
	static constexpr std::uint32_t serialized_value_size_hint() {
		return static_cast<std::uint32_t>(sizeof(serialized_value_type));
	}
	static std::string name() { return "Z32783"; }
	static Z32783 sample(std::mt19937_64& rng) {
		std::uniform_int_distribution<std::uint64_t> dist(0, characteristic() - 1);
		return Z32783{dist(rng)};
	}

	// Truthiness (so you can write: if (x) ...)
	explicit operator bool() const { return v_ != 0; }

	// Comparisons (hidden friends) — symmetric and allow implicit conversions of both sides
	friend bool operator==(const Z32783& a, const Z32783& b) { return a.v_ == b.v_; }
	friend bool operator!=(const Z32783& a, const Z32783& b) { return a.v_ != b.v_; }

	// ---------------- Arithmetic ----------------
	// Unary minus
	Z32783 operator-() const {
		if (v_ == 0) return Z32783{0};
		Z32783 out; out.v_ = static_cast<std::uint32_t>(MOD - v_);
		return out;
	}
	// Addition (mod MOD)
	Z32783 operator+(const Z32783& o) const {
		std::uint32_t s = v_ + o.v_;                 // wrap is defined for unsigned
		if (s >= MOD) s -= MOD;                      // reduce
		Z32783 out; out.v_ = s; return out;
	}
	Z32783& operator+=(const Z32783& o) {
		v_ += o.v_;
		if (v_ >= MOD) v_ -= MOD;
		return *this;
	}

	// Subtraction (mod MOD)
	Z32783 operator-(const Z32783& o) const {
		std::uint32_t d = v_ - o.v_;
		if (v_ < o.v_) d += MOD;                     // borrow -> add MOD back
		Z32783 out; out.v_ = d; return out;
	}
	Z32783& operator-=(const Z32783& o) {
		if (v_ < o.v_) v_ = static_cast<std::uint32_t>(v_ + MOD - o.v_);
		else            v_ = static_cast<std::uint32_t>(v_ - o.v_);
		return *this;
	}

	// Multiplication (mod MOD)
	Z32783 operator*(const Z32783& o) const {
		// Safe with 32-bit: (MOD-1)^2 < 2^31
		std::uint32_t p = v_ * o.v_;
		p %= MOD;
		Z32783 out; out.v_ = p; return out;
	}
	Z32783& operator*=(const Z32783& o) { return *this = *this * o; }

	// Division (mod MOD): multiply by inverse; throws on division by zero
	Z32783 operator/(const Z32783& o) const {
		if (o.v_ == 0) throw std::domain_error("division by zero");
		return *this * o.inv();
	}
	Z32783& operator/=(const Z32783& o) {
		return *this = *this / o;
	}

	// Multiplicative inverse via Fermat: a^(MOD-2) mod MOD
	// Unrolled for (MOD-2) = 32781 = 0b1_0000_0000_0000_1101 (bits 0,2,3,15 set)
	Z32783 inv() const {
		if (v_ == 0) throw std::domain_error("division by zero");

		Z32783 a(v_);      // v^1
		a = a * a;              // v^2   (bit 1 = 0)
		a = a * a;              // v^4   (move to bit 2)

		Z32783 r = Z32783(v_) * a; // r = v^1 * v^4 = v^5 (bit 2 = 1)

		a = a * a;              // v^8   (bit 3)
		r = r * a;              // r = v^5 * v^8 = v^13 (bit 3 = 1)

		// bits 4..14 = 0 → 11 squarings to reach 2^15
		a = a * a; // v^16
		a = a * a; // v^32
		a = a * a; // v^64
		a = a * a; // v^128
		a = a * a; // v^256
		a = a * a; // v^512
		a = a * a; // v^1024
		a = a * a; // v^2048
		a = a * a; // v^4096
		a = a * a; // v^8192
		a = a * a; // v^16384

		// bit 15 = 1 (2^15 = 32768)
		a = a * a; // v^32768
		r = r * a; // r = v^13 * v^32768 = v^32781

		return r;  // v^(MOD-2) = v^{-1} mod MOD
	}

	static void write_serialized_value(std::ostream& out, const serialized_value_type& value) {
		out.write(reinterpret_cast<const char*>(&value), sizeof(value));
	}

	static bool read_serialized_value(std::istream& in, serialized_value_type& value) {
		in.read(reinterpret_cast<char*>(&value), sizeof(value));
		return static_cast<bool>(in);
	}

	static void write_value(std::ostream& out, const Z32783& value) {
		write_serialized_value(out, value.value());
	}

	static bool read_value(std::istream& in, Z32783& value) {
		serialized_value_type storage{};
		if (!read_serialized_value(in, storage)) {
			return false;
		}
		value = Z32783{storage};
		return true;
	}

	std::int64_t signed_representative() const {
		return v_ <= MOD / 2
			? static_cast<std::int64_t>(v_)
			: static_cast<std::int64_t>(v_) - static_cast<std::int64_t>(MOD);
	}
};

// Stream operators for easy I/O with iostreams
inline std::ostream& operator<<(std::ostream& os, const Z32783& a) {
	return os << a.value();
}
inline std::istream& operator>>(std::istream& is, Z32783& a) {
	long long x; // accepts negatives, etc.
	if (is >> x) a = Z32783{x};
	return is;
}
