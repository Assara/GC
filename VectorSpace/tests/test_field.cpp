#include <iostream>
#include <stdexcept>
#include "VectorSpace/Field.hpp"
#if __has_include(<boost/multiprecision/cpp_int.hpp>)
#include "Q.hpp"
#define GC_TEST_RATIONAL_FIELD
#endif

static_assert(std::same_as<SmallSignedInt, std::int8_t>);
static_assert(sizeof(SmallSignedInt) == 1);
static_assert(!VectorSpace::Field<int>);
static_assert(!VectorSpace::Field<std::string>);

void check(bool passed) {
    if (!passed) throw std::runtime_error("small signed field multiplication mismatch");
}

template <VectorSpace::Field F>
void check_values(const F& value) {
    for (int i = -128; i <= 127; ++i) {
        const auto small = static_cast<SmallSignedInt>(i);
        const F expected = value * F{i};
        check(value * small == expected);
        check(small * value == expected);
        F inplace = value;
        check(&(inplace *= small) == &inplace);
        check(inplace == expected);
    }
    // Ordinary integers must not select a narrowing signed-byte overload.
    for (int i : {-1000, -129, 128, 256, 1000}) {
        check(value * i == value * F{i});
        F inplace = value;
        inplace *= i;
        check(inplace == value * F{i});
    }
}

template <VectorSpace::Field F>
void finite_field() {
    for (const auto& x : {F{0},F{1},F{-1},F{2},F{-2},F{F::characteristic()/2}})
        check_values(x);
    std::cout << F::name() << ": all signed-byte multipliers passed\n";
}
int main() {
    finite_field<Z32783>();
    finite_field<Z2179564669>();
    finite_field<Z2305843009213693951>();
    finite_field<Z34821139123>();
    finite_field<Z4294967291>();
#ifdef GC_TEST_RATIONAL_FIELD
    check_values(Q{0}); check_values(Q{1}); check_values(Q{-1});
    check_values(Q{2}/Q{3}); check_values(Q{-7}/Q{11});
    std::cout << "Q: all signed-byte multipliers passed\n";
#else
    std::cout << "Q checks skipped: Boost headers unavailable\n";
#endif
}
