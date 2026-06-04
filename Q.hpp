#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

#include <istream>
#include <ostream>
#include <stdexcept>
#include <string>
#include <type_traits>

class Q {
public:
	using BigInt = boost::multiprecision::cpp_int;
	using Storage = boost::rational<BigInt>;

private:
	Storage value_;

public:
	constexpr Q() : value_(0) {}

	template <class T, std::enable_if_t<std::is_integral_v<T>, int> = 0>
	/* not explicit */ Q(T x) : value_(static_cast<long long>(x)) {}

	Q(const BigInt& numerator) : value_(numerator) {}

	Q(const BigInt& numerator, const BigInt& denominator)
		: value_(numerator, denominator) {
		if (denominator == 0) {
			throw std::domain_error("division by zero");
		}
	}

	explicit Q(const Storage& x) : value_(x) {}

	const Storage& raw() const { return value_; }
	const BigInt& numerator() const { return value_.numerator(); }
	const BigInt& denominator() const { return value_.denominator(); }

	explicit operator bool() const { return value_ != 0; }

	friend bool operator==(const Q& a, const Q& b) { return a.value_ == b.value_; }
	friend bool operator!=(const Q& a, const Q& b) { return a.value_ != b.value_; }

	Q operator-() const { return Q(-value_); }

	Q operator+(const Q& o) const { return Q(value_ + o.value_); }
	Q& operator+=(const Q& o) {
		value_ += o.value_;
		return *this;
	}

	Q operator-(const Q& o) const { return Q(value_ - o.value_); }
	Q& operator-=(const Q& o) {
		value_ -= o.value_;
		return *this;
	}

	Q operator*(const Q& o) const { return Q(value_ * o.value_); }
	Q& operator*=(const Q& o) {
		value_ *= o.value_;
		return *this;
	}

	Q operator/(const Q& o) const {
		if (o.value_ == 0) {
			throw std::domain_error("division by zero");
		}
		return Q(value_ / o.value_);
	}
	Q& operator/=(const Q& o) {
		return *this = *this / o;
	}

	Q inv() const {
		if (value_ == 0) {
			throw std::domain_error("division by zero");
		}
		return Q(value_.denominator(), value_.numerator());
	}
};

inline std::ostream& operator<<(std::ostream& os, const Q& q) {
	os << q.numerator();
	if (q.denominator() != 1) {
		os << "/" << q.denominator();
	}
	return os;
}

inline std::istream& operator>>(std::istream& is, Q& q) {
	std::string token;
	if (!(is >> token)) {
		return is;
	}

	const std::size_t slash = token.find('/');
	if (slash == std::string::npos) {
		q = Q(Q::BigInt(token));
		return is;
	}

	const std::string numerator = token.substr(0, slash);
	const std::string denominator = token.substr(slash + 1);
	q = Q(Q::BigInt(numerator), Q::BigInt(denominator));
	return is;
}
