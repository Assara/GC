#pragma once

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/rational.hpp>

#include <istream>
#include <ostream>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>

class Q {
public:
	using arbInt = boost::multiprecision::cpp_int;
	using serialized_value_type = std::pair<arbInt, arbInt>;
	using Storage = boost::rational<arbInt>;

private:
	Storage value_;

public:
	constexpr Q() : value_(0) {}

	template <class T, std::enable_if_t<std::is_integral_v<T>, int> = 0>
	/* not explicit */ Q(T x) : value_(static_cast<long long>(x)) {}

	Q(const arbInt& numerator) : value_(numerator) {}

	Q(const arbInt& numerator, const arbInt& denominator)
		: value_(numerator, denominator) {
		if (denominator == 0) {
			throw std::domain_error("division by zero");
		}
	}

	Q(const serialized_value_type& serialized)
		: Q(serialized.first, serialized.second) {}

	explicit Q(const Storage& x) : value_(x) {}

	static constexpr std::uint32_t characteristic() { return 0; }
	static constexpr std::uint32_t serialized_value_size_hint() { return 0; }
	static std::string name() { return "Q"; }
	static Q sample(std::mt19937_64& rng) {
		std::uniform_int_distribution<int> dist(-1, 1);
		return Q{dist(rng)};
	}

	const Storage& raw() const { return value_; }
	serialized_value_type value() const {
		return {value_.numerator(), value_.denominator()};
	}
	const arbInt& numerator() const { return value_.numerator(); }
	const arbInt& denominator() const { return value_.denominator(); }

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

	static void write_serialized_value(std::ostream& out, const serialized_value_type& value) {
		write_arb_int(out, value.first);
		write_arb_int(out, value.second);
	}

	static bool read_serialized_value(std::istream& in, serialized_value_type& value) {
		return read_arb_int(in, value.first) && read_arb_int(in, value.second);
	}

	static void write_value(std::ostream& out, const Q& value) {
		write_serialized_value(out, value.value());
	}

	static bool read_value(std::istream& in, Q& value) {
		serialized_value_type storage{};
		if (!read_serialized_value(in, storage)) {
			return false;
		}
		value = Q{storage};
		return true;
	}

private:
	static void write_arb_int(std::ostream& out, const arbInt& value) {
		const std::string text = value.convert_to<std::string>();
		const std::uint64_t size = static_cast<std::uint64_t>(text.size());
		out.write(reinterpret_cast<const char*>(&size), sizeof(size));
		out.write(text.data(), static_cast<std::streamsize>(text.size()));
	}

	static bool read_arb_int(std::istream& in, arbInt& value) {
		std::uint64_t size = 0;
		in.read(reinterpret_cast<char*>(&size), sizeof(size));
		if (!in) {
			return false;
		}
		std::string text(size, '\0');
		in.read(text.data(), static_cast<std::streamsize>(text.size()));
		if (!in) {
			return false;
		}
		std::istringstream stream(text);
		stream >> value;
		return static_cast<bool>(stream);
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
		q = Q(Q::arbInt(token));
		return is;
	}

	const std::string numerator = token.substr(0, slash);
	const std::string denominator = token.substr(slash + 1);
	q = Q(Q::arbInt(numerator), Q::arbInt(denominator));
	return is;
}
