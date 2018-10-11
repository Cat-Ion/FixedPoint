/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file at the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#include <iostream>
#include <stdint.h>
#include <UnitTest++.h>
#include "FixedPoint/FixedPoint.hpp"

SUITE(MultiwordInteger) {
  TEST(zero_construction) {
    MultiwordInteger<4, uint16_t> a(int64_t(0));
    CHECK_EQUAL(double(0), a.operator double());
  }

  TEST(lsb_construction) {
    MultiwordInteger<4, uint16_t> a(int64_t(1));
    CHECK_EQUAL(double(1), a.operator double());
  }

  TEST(msb_construction) {
    MultiwordInteger<4, uint16_t> a(int64_t(-1));
    CHECK_EQUAL(double(-1), a.operator double());
  }

  TEST(storagetype_construction) {
    MultiwordInteger<2, uint16_t> a;

    // Positive values
    a = uint16_t(5);
    CHECK_EQUAL(double(a), 5.);

    // Negative values get sign-extended
    a = uint16_t(-5);
    CHECK_EQUAL(double(a), -5.);
  }

  TEST(double_construction) {
    MultiwordInteger<20, uint16_t> a;
    for (unsigned i = 0; i < 200; i++) {
      double b = pow(2., i);
      a = b;
      CHECK_EQUAL(b, a.operator double());
    }
  }

  TEST(other_construction) {
    MultiwordInteger<4, uint16_t> a(int64_t(0x7FFF7FFF7FFF7FFF));
    MultiwordInteger<2, uint16_t> b(a);
    // Assignment should truncate
    CHECK_EQUAL(b.operator double(), double(0x7FFF7FFF));

    a = int64_t(0x7FFF7FFF80000000);
    b = a;
    // This should be negative
    CHECK(b.is_negative());
    CHECK_EQUAL(double(b), -double(0x80000000));

    b = int64_t(0x7FFFFFFF);
    a = b;
    // Assignment from a smaller type should be exact
    CHECK_EQUAL(double(a), double(b));

    // Assignment from a type that doesn't fill the storage completely
    MultiwordInteger<3, uint8_t> c(int64_t(0xC0FFEE));
    a = c;
    CHECK_EQUAL(double(a), double(c));

    MultiwordInteger<1, uint32_t> d = MultiwordInteger<2, uint8_t>(-5.);
    CHECK_EQUAL(double(d), -5.);
  }

  TEST(limits) {
    MultiwordInteger<2, uint16_t> maxv(MultiwordInteger<2, uint16_t>::maxVal);
    MultiwordInteger<2, uint16_t> minv(MultiwordInteger<2, uint16_t>::minVal);
    MultiwordInteger<2, uint16_t> a(int64_t(0x80000000));

    CHECK_EQUAL(double(a), double(minv));

    a=int64_t(0x7FFFFFFF);
    CHECK_EQUAL(double(a), double(maxv));
  }

  template<typename T>
  void addition_perform() {
      T a(0.), b(1.), c(-1.), d(2.);

      // Adding zero
      CHECK_EQUAL(double(a+b), double(b));

      // Basic addition
      CHECK_EQUAL(double(b+b), double(d));

      // Adding with overflow
      CHECK_EQUAL(double(a+b+c), double(a));

      // Commutativity
      CHECK_EQUAL(double(b+d), double(d+b));

      b = 0.;
      CHECK_EQUAL(double(a++), double(b));
      CHECK_EQUAL(double(++a), double(d));
  }

  TEST(addition) {
      addition_perform<MultiwordInteger<2, uint16_t>>();
      addition_perform<MultiwordInteger<1, uint16_t>>();
  }

  template<typename T>
  void subtraction_perform() {
      T a(0.), b(1.), c(-1.), d(3.);

      // Subtracting from zero should yield negative
      CHECK_EQUAL(double(a-b), double(c));

      // Subtracting itself should yield zero
      CHECK_EQUAL(double(b-b), double(a));

      // Basic subtraction
      CHECK_EQUAL(double(d-b-b), double(b));

      // Anticommutativity
      CHECK_EQUAL(double(a-b), double(a-(b-a)));

      d = 2.;
      c = 2.;
      CHECK_EQUAL(double(d--), double(c));
      CHECK_EQUAL(double(--d), double(a));
  }

  TEST(subtraction) {
      subtraction_perform<MultiwordInteger<2, uint16_t>>();
      subtraction_perform<MultiwordInteger<1, uint16_t>>();
  }

  template<typename T>
  void negation_perform() {
    T a(1.), b(-1.), c(-1.);

    CHECK_EQUAL(double(a), double(-b));

    c.negate();
    CHECK_EQUAL(double(a), double(c));
  }

  TEST(negation) {
      negation_perform<MultiwordInteger<2, uint16_t>>();
      negation_perform<MultiwordInteger<1, uint16_t>>();
  }

  template<typename T>
  void multiplication_perform() {
    T a(0.), b(1.), c(-1.), d(2.);

    // Multiplying with zero yields zero
    CHECK_EQUAL(double(a*d), double(a));

    // Multiplying with 1 yields the other operand
    CHECK_EQUAL(double(b*d), double(d));

    // Multiplying with -1 yields -d
    CHECK_EQUAL(double(d*c), double(-d));

    // Commutativity
    CHECK_EQUAL(double(c*d), double(d*c));

    // Results out of range of the smaller type should still work in the resulting type
    b = pow(2., T::numWords * T::storageSize - 2.);
    MultiwordInteger<T::numWords*2, typename T::storageType> e(pow(2., 2*T::numWords*T::storageSize-4.));
    CHECK_EQUAL(double(b*b), double(e));

    // But overflow when assigned to a smaller type
    T f(b);
    f *= b;
    CHECK_EQUAL(double(f), double(a));
  }

  TEST(multiplication) {
      multiplication_perform<MultiwordInteger<2, uint16_t>>();
      multiplication_perform<MultiwordInteger<1, uint16_t>>();
  }

  template<typename T>
  void division_perform() {
    T a(1.), b(2.), c, d, zero(0.), unity(1.);

    // Dividing by larger values yields zero for positive values, -1 for negative values
    CHECK_EQUAL(double(a/b), double(zero));
    CHECK_EQUAL(double(a/(-b)), double(-unity));

    // Dividing by 1 yields the numerator
    CHECK_EQUAL(double(b/unity), double(b));

    // x/x == 1
    CHECK_EQUAL(double(b/b), double(unity));

    if (T::numWords > 1) {
        // Test again, for larger values
        a = pow(2., double(T::storageSize))*5;
        b = pow(2., double(T::storageSize));
        CHECK_EQUAL(double(a/a), double(unity));
    }

    // If there is no remainder, a/b == c <-> a/c == b
    a = 5.;
    b = 1.;
    c = 5.;
    CHECK_EQUAL(double(a/b), double(c));
    CHECK_EQUAL(double(a/c), double(b));
    CHECK_EQUAL(double((-a)/b), double(-c));
    CHECK_EQUAL(double(a/(-b)), double(-c));
    CHECK_EQUAL(double((-a)/(-b)), double(c));

    // Test with remainder
    a = 5.;
    b = 2.;
    c = 3.;
    CHECK_EQUAL(double(a/-b), double(-c));
    CHECK_EQUAL(double(a/-c), double(-b));
    CHECK_EQUAL(double(-a/b), double(-c));
    CHECK_EQUAL(double(-a/c), double(-b));

    if (T::numWords > 1) {
        // And for large values
        a = 5*pow(2., T::storageSize);
        b = 3*pow(2., T::storageSize);
        c = 1.;
        d = 2.;
        CHECK_EQUAL(double(a/-b), double(-d));
        CHECK_EQUAL(double((-a)/b), double(-d));
        CHECK_EQUAL(double((-a)/(-b)), double(c));
    }

    // Division by zero
    a = 1.;
    b = 0.;
    c = T::_maxVal();
    d = T::_minVal();
    CHECK_EQUAL(double(a/b), double(c));
    CHECK_EQUAL(double((-a)/b), double(d));
  }

  TEST(division) {
      division_perform<MultiwordInteger<2, uint16_t>>();
      division_perform<MultiwordInteger<1, uint32_t>>();

      MultiwordInteger<1, uint32_t> a;
      MultiwordInteger<2, uint32_t> b;

      a = 1.;
      b = 2.;
      CHECK_EQUAL(double(a/b), 0.);
      CHECK_EQUAL(double((-a)/b), -1);
      CHECK_EQUAL(double(a/(-b)), 0.);
  }

  template<typename T>
  void modulo_perform() {
    T a, b, c, d, e;

    // Without remainder, positive and negative
    a = 5.;
    b = 1.;
    c = 0.;
    d = a/b * b;
    e = a%b;
    CHECK_EQUAL(double(e), double(c));
    CHECK_EQUAL(double(d + e), double(a));

    d = a / (-b) * (-b);
    e = a % (-b);
    CHECK_EQUAL(double(e), double(c));
    CHECK_EQUAL(double(d + e), double(a));

    d = (-a) / b * b;
    e = (-a) % b;
    CHECK_EQUAL(double(e), double(c));
    CHECK_EQUAL(double(d + e), double(-a));

    d = (-a) / (-b) * (-b);
    e = (-a) % (-b);
    CHECK_EQUAL(double(e), double(-c));
    CHECK_EQUAL(double(d + e), double(-a));

    // With remainder
    a = 7.;
    b = 3.;
    c = 1.;

    d = a/b*b;
    e = a%b;
    CHECK_EQUAL(double(e), double(c));
    CHECK_EQUAL(double(d + e), double(a));

    d = (-a) / (-b) * (-b);
    e = (-a) % (-b);
    CHECK_EQUAL(double(e), double(-c));
    CHECK_EQUAL(double(d + e), double((-a)));

    c = 2.;
    d = a/(-b)*(-b);
    e = a%(-b);
    CHECK_EQUAL(double(e), double(-c));
    CHECK_EQUAL(double(d + e), double(a));

    d = (-a)/b*b;
    e = (-a)%b;
    CHECK_EQUAL(double(e), double(c));
    CHECK_EQUAL(double(d + e), double(-a));

    // Large divisor
    a = 65536.*2 - 1;
    b = 65536.;

    c = 65535.;
    d = a/b*b;
    e = a%b;
    CHECK_EQUAL(double(e), double(c));
    CHECK_EQUAL(double(d + e), double(a));

    d = (-a) / (-b) * (-b);
    e = (-a) % (-b);
    CHECK_EQUAL(double(e), double(-c));
    CHECK_EQUAL(double(d + e), double(-a));

    c = 1.;
    d = (-a) / b * b;
    e = (-a) % b;
    CHECK_EQUAL(double(e), double(c));
    CHECK_EQUAL(double(d + e), double(-a));

    d = a / (-b) * (-b);
    e = a % (-b);
    CHECK_EQUAL(double(e), double(-c));
    CHECK_EQUAL(double(d + e), double(a));

    // Division of zero returns 0
    a = 0.;
    b = 1.;
    d = a / b * b;
    e = a % b;
    CHECK_EQUAL(double(d), double(a));
    CHECK_EQUAL(double(e), double(a));

    // Division by zero returns 0
    a = 1.;
    b = 0.;
    e = a % b;
    CHECK_EQUAL(double(e), double(b));
  }

  TEST(modulo) {
      modulo_perform<MultiwordInteger<2, uint16_t>>();
      modulo_perform<MultiwordInteger<1, uint32_t>>();

      MultiwordInteger<1, uint32_t> a;
      MultiwordInteger<2, uint32_t> b;

      a = 1.;
      b = 2.;
      CHECK_EQUAL(double(a%b), 1.);
      CHECK_EQUAL(double((-a)%b), 1.);
      CHECK_EQUAL(double(a%(-b)), 1.);
  }

  TEST(comparison) {
    MultiwordInteger<2, uint16_t> a, b;

    a = 1.;
    b = 0.;

    CHECK(a > b);
    CHECK(a >= b);
    CHECK(b < a);
    CHECK(b <= a);

    CHECK(!(a < a));

    CHECK(a != b);
    CHECK_EQUAL(double(a), double(a));

    CHECK(!(a == b));
    CHECK(!(a != a));

    b = -1.;
    CHECK(b < a);
  }

  template<typename T>
  void shifts_perform() {
    T a = 65536., b = 0.;

    // Test shifting all bits out, when positive
    a >>= 17;
    CHECK_EQUAL(double(a), double(b));

    // Test shifting 0
    a <<= 1;
    CHECK_EQUAL(double(a), double(b));

    // Huge shift
    a = 1.;
    a <<= 100;
    CHECK_EQUAL(double(a), double(b));

    a = 1.;
    // Test normal shift
    b = 4.;
    CHECK_EQUAL(double(a<<1), double(b>>1));
    a <<= 1;
    b >>= 1;
    CHECK_EQUAL(double(a), double(b));

    // Shift over boundary
    a <<= 17;
    b = 65536.*4;
    CHECK_EQUAL(double(a), double(b));

    a <<= 12;

    // Test right shift
    a >>= 1;
    b = 536870912.;
    CHECK_EQUAL(double(a), double(b));

    // Test large right shift
    a >>= 29;
    b = 1.;
    CHECK_EQUAL(double(a), double(b));

    // Result should be negative
    a <<= 31;
    b = T::minVal;
    CHECK_EQUAL(double(a), double(b));

    // Result should be sign extended, and negative
    a >>= 30;
    b = -2.;
    CHECK_EQUAL(double(a), double(b));

    // Too large shift should still sign extend
    a >>= 300;
    b = -1.;
    CHECK_EQUAL(double(a), double(b));

    // Huge shift for positive values
    a = 555.;
    a >>= 300;
    b = 0.;
    CHECK_EQUAL(double(a), double(b));

    // Shift by storageSize
    a = 65536.;
    a >>= 16;
    b = 1.;
    CHECK_EQUAL(double(a), double(b));

    // Shift by storageSize, negative number
    a = -65536.;
    a >>= 16;
    b = -1.;
    CHECK_EQUAL(double(a), double(b));
  }
  TEST(shifts) {
      shifts_perform<MultiwordInteger<2, uint16_t>>();
      shifts_perform<MultiwordInteger<1, uint32_t>>();
  }

  template<typename T>
  void bitwise_perform() {
      T a, b, c;
      a = int64_t(0xFFFFFFFF);
      b = int64_t(0xFFFFFFFF);
      c = int64_t(0);

      CHECK_EQUAL(double(a^a), double(c));

      CHECK_EQUAL(double(a&a), double(a));
      CHECK_EQUAL(double(a&c), double(c));

      CHECK_EQUAL(double(a|c), double(a));
      CHECK_EQUAL(double(c|c), double(c));

      CHECK_EQUAL(double(~a), double(c));
      CHECK_EQUAL(double(~c), double(a));
  }

  TEST(bitwise) {
      bitwise_perform<MultiwordInteger<1, uint32_t>>();
      bitwise_perform<MultiwordInteger<2, uint16_t>>();
  }

  TEST(conversions) {
    MultiwordInteger<4, uint16_t> a;

    MultiwordInteger<1, uint32_t> b; // Larger backing type, less bits
    MultiwordInteger<2, uint32_t> c; // Larger backing type, same number of bits
    MultiwordInteger<4, uint32_t> d; // Larger backing type, more bits

    MultiwordInteger<2, uint16_t> e; // Same backing type, less bits
    MultiwordInteger<4, uint16_t> f; // Same backing type, same number of bits
    MultiwordInteger<4, uint16_t> g; // Same backing type, more bits

    MultiwordInteger<4, uint8_t> h; // Smaller backing type, less bits
    MultiwordInteger<8, uint8_t> i; // Smaller backing type, same number of bits
    MultiwordInteger<16, uint8_t> j; // Smaller backing type, more bits

    // Small value conversion
    a = 5.;
    b = a;
    c = a;
    d = a;
    e = a;
    f = a;
    g = a;
    h = a;
    i = a;
    j = a;
    CHECK_EQUAL(double(b), double(a));
    CHECK_EQUAL(double(c), double(a));
    CHECK_EQUAL(double(d), double(a));
    CHECK_EQUAL(double(e), double(a));
    CHECK_EQUAL(double(f), double(a));
    CHECK_EQUAL(double(g), double(a));
    CHECK_EQUAL(double(h), double(a));
    CHECK_EQUAL(double(i), double(a));
    CHECK_EQUAL(double(j), double(a));

    // Larger value conversion
    a = 65536.;
    b = a; c = a; d = a; e = a; f = a; g = a; h = a; i = a; j = a;
    CHECK_EQUAL(double(b), double(a));
    CHECK_EQUAL(double(c), double(a));
    CHECK_EQUAL(double(d), double(a));
    CHECK_EQUAL(double(e), double(a));
    CHECK_EQUAL(double(f), double(a));
    CHECK_EQUAL(double(g), double(a));
    CHECK_EQUAL(double(h), double(a));
    CHECK_EQUAL(double(i), double(a));
    CHECK_EQUAL(double(j), double(a));

    // Negative conversion
    a = -1.;
    b = a;
    c = a;
    d = a;
    e = a;
    f = a;
    g = a;
    h = a;
    i = a;
    j = a;
    CHECK_EQUAL(double(b), double(a));
    CHECK_EQUAL(double(c), double(a));
    CHECK_EQUAL(double(d), double(a));
    CHECK_EQUAL(double(e), double(a));
    CHECK_EQUAL(double(f), double(a));
    CHECK_EQUAL(double(g), double(a));
    CHECK_EQUAL(double(h), double(a));
    CHECK_EQUAL(double(i), double(a));
    CHECK_EQUAL(double(j), double(a));

    // Overflow
    a = pow(2., 32.);
    b = a; e = a; h = a;
    CHECK(!b);
    CHECK(!e);
    CHECK(!h);
  }

  TEST(casts) {
      MultiwordInteger<8, uint8_t> a(int64_t(0x1122334455667788));
      MultiwordInteger<1, uint32_t> b(int64_t(0x55667788));

      CHECK_EQUAL(double(a), double(int64_t(0x1122334455667788UL)));
      CHECK_EQUAL(int64_t(a), int64_t(0x1122334455667788UL));
      CHECK_EQUAL(int32_t(a), 0x55667788);
      CHECK_EQUAL(int16_t(a), 0x7788);
      CHECK_EQUAL(int8_t(a),  int8_t(0x88));

      CHECK_EQUAL(double(b), double(int64_t(0x55667788UL)));
      CHECK_EQUAL(int64_t(b), int64_t(0x55667788UL));
      CHECK_EQUAL(int32_t(b), 0x55667788);
      CHECK_EQUAL(int16_t(b), 0x7788);
      CHECK_EQUAL(int8_t(b),  int8_t(0x88));
  }

  TEST(raw_data) {
      MultiwordInteger<1, uint32_t> a = 5. + 6.*256 + 7.*256*256 + 8.*256*256*256;
      uint8_t data[4];

      a.get_raw(data);
      CHECK_EQUAL(data[0], 5);
      CHECK_EQUAL(data[1], 6);
      CHECK_EQUAL(data[2], 7);
      CHECK_EQUAL(data[3], 8);

      MultiwordInteger<2, uint16_t> b = 5. + 6.*256 + 7.*256*256 + 8.*256*256*256;
      b.get_raw(data);
      CHECK_EQUAL(data[0], 5);
      CHECK_EQUAL(data[1], 6);
      CHECK_EQUAL(data[2], 7);
      CHECK_EQUAL(data[3], 8);
  }
 
  TEST(helpers) {
      MultiwordInteger<2, uint32_t> a(int64_t(1));
      CHECK_EQUAL(a.leading_zeros(), 63);
      CHECK_EQUAL(FixedPointHelpers::nlz_constexpr(int32_t(1)), 31);
      CHECK_EQUAL(FixedPointHelpers::nlz_constexpr(int32_t(0)), 32);
      CHECK_EQUAL(FixedPointHelpers::nlz(uint64_t(1)), 63);
      CHECK_EQUAL(FixedPointHelpers::nlz(uint32_t(1)), 31);
      CHECK_EQUAL(FixedPointHelpers::nlz(uint16_t(1)), 15);
      CHECK_EQUAL(FixedPointHelpers::nlz(uint8_t(1)), 7);
      CHECK_EQUAL(FixedPointHelpers::ilogb(0.), INT_MIN);
      CHECK_EQUAL(FixedPointHelpers::ilogb(0.5), -1);
      CHECK_EQUAL(FixedPointHelpers::ilogb(2.), 1);
      CHECK_EQUAL(FixedPointHelpers::ilogb(-2.), 1);
      CHECK_EQUAL(FixedPointHelpers::ilogb(NAN), 0);
      CHECK_EQUAL(FixedPointHelpers::ilogb(INFINITY), INT_MAX);
      CHECK_EQUAL(FixedPointHelpers::ilogb(0.f), INT_MIN);
      CHECK_EQUAL(FixedPointHelpers::ilogb(0.5f), -1);
      CHECK_EQUAL(FixedPointHelpers::ilogb(2.f), 1);
      CHECK_EQUAL(FixedPointHelpers::ilogb(-2.f), 1);
      CHECK_EQUAL(FixedPointHelpers::ilogb(float(NAN)), 0);
      CHECK_EQUAL(FixedPointHelpers::ilogb(float(INFINITY)), INT_MAX);
  }
}

SUITE(FixedPoint) {
  TEST(casts) {
    FixedPoint<15, 16, uint32_t> a;
    a = -0.5;
    CHECK_EQUAL(int32_t(a), -1);
  }

  TEST(construction) {
    FixedPoint<3, 28, uint16_t> a, b;
    FixedPoint<5, 26, uint32_t> c;
    FixedPoint<5, 26, uint32_t> d;

    FixedPoint<31, 0>::StorageType s(int64_t(1));

    a = 5.;
    CHECK_EQUAL(double(a), 5.);

    a = int(5);
    CHECK_EQUAL(double(a), 5.);

    b = a;
    CHECK_EQUAL(double(b), 5.);

    c = a;
    CHECK_EQUAL(double(c), 5.);

    d = a;
    CHECK_EQUAL(double(d), 5.);

    FixedPoint<31, 0> e((s));
    CHECK_EQUAL(double(e), 1.);

    FixedPoint<30, 1> f((s));
    CHECK_EQUAL(double(f), 0.5);

    FixedPoint<32, -1> g(a);
    CHECK_EQUAL(double(g), 6); // Round up

    a = 100.;
    CHECK_EQUAL(double(a), double(FixedPoint<3, 28, uint16_t>::maxVal));

    a = -100.;
    CHECK_EQUAL(double(a), double(FixedPoint<3, 28, uint16_t>::minVal));

    a = int(100);
    CHECK_EQUAL(double(a), double(FixedPoint<3, 28, uint16_t>::maxVal));

    a = int(-100);
    CHECK_EQUAL(double(a), double(FixedPoint<3, 28, uint16_t>::minVal));
  }

  TEST(limits) {
    FixedPoint<5, 10, uint16_t> a;
    double max = pow(2., 5) - pow(2., -(double)a.fractionalWidth);
    double min = -pow(2., 5.);
    double smallest = pow(2., -(double)a.fractionalWidth);

    a = FixedPoint<5, 10, uint16_t>::maxVal;
    CHECK_EQUAL(double(a), max);

    a = FixedPoint<5, 10, uint16_t>::minVal;
    CHECK_EQUAL(double(a), min);

    a = FixedPoint<5, 10, uint16_t>::smallestVal;
    CHECK_EQUAL(double(a), smallest);
  }

  TEST(addition) {
    FixedPoint<5, 10, uint16_t> a, b, c;

    a = 1.;
    b = pow(2., -10.);
    c = (1. + pow(2., -10));
    CHECK_EQUAL(double((a + b)), double(c));
    CHECK_EQUAL(double((b + a)), double(c));

    a = -1.;
    b = -2.;
    c = -3.;
    CHECK_EQUAL(double((a+b)), double(c));
  }

  TEST(subtraction) {
    FixedPoint<5, 10, uint16_t> a, b, c;

    a = 1.;
    b = pow(2., -10.);
    c = (1. - pow(2., -10));
    CHECK_EQUAL(double((a - b)), double(c));
    CHECK_EQUAL(double((b - a)), double(-c));
  }

  TEST(multiplication) {
    FixedPoint<5, 10, uint16_t> a, b, c;

    a = 1.;
    b = pow(2., -10.);
    CHECK_EQUAL(double(a * b), double(b));

    a = 5.;
    b = 0.;
    CHECK_EQUAL(double(a * b), double(b));

    b = 2.;
    c = 10.;
    CHECK_EQUAL(double(a * b), double(c));
    CHECK_EQUAL(double(b * a), double(c));

    CHECK_EQUAL(double(b * 5), double(c));
  }

  TEST(division) {
    FixedPoint<5, 10, uint16_t> a, b, c;

    a = 1.;
    b = pow(2., -4.);
    c = pow(2., 4.);
    CHECK_EQUAL(double((a / b)), double(c));
    CHECK_EQUAL(double((a / c)), double(b));
    CHECK_EQUAL(double((a / (-b))), double(-c));
    CHECK_EQUAL(double((a / (-c))), double(-b));
    CHECK_EQUAL(double(((-a) / (-b))), double(c));
    CHECK_EQUAL(double(((-a) / (-c))), double(b));
    CHECK_EQUAL(double(a / 16), double(b));
    CHECK_EQUAL(double(b % a), double(b));

    FixedPoint<4, 91, uint32_t> y32(-0.13779029068463858), x32(-0.99046142569665119);
    FixedPoint<4, 91, uint32_t> d32 = y32/x32;

    FixedPoint<4, 91, uint16_t> y16(y32), x16(x32), d16(y16/x16);
    FixedPoint<4, 91, uint8_t> y8(y32), x8(x32), d8(y16/x16);
    FixedPoint<4, 91> d8_32(d8), d16_32(d16);

    CHECK_EQUAL(double(double(y32)), double(double(y16)));
    CHECK_EQUAL(double(double(x32)), double(double(x16)));
    CHECK_EQUAL(double(d8_32), double(d32));
    CHECK_EQUAL(double(d16_32), double(d32));
  }

  template <typename T>
  void test_comparison_perform() {
    T a, b;

    a = 1.;
    b = 0.;

    CHECK(a > b);
    CHECK(a >= b);
    CHECK(b < a);
    CHECK(b <= a);

    CHECK(a != b);
    CHECK_EQUAL(double(a), double(a));

    CHECK(!(a == b));
    CHECK(!(a != a));

    b = -1.;
    a = -0.5;
    CHECK(b < a);
    CHECK(a > b);
    CHECK(b <= a);
    CHECK(a >= b);

    CHECK(T(1.).is_positive());
    CHECK(! T(0.).is_positive());
    CHECK(! T(-1.).is_positive());
  }

  TEST(comparison) {
    test_comparison_perform<FixedPoint<3, 28, uint32_t>>();
    test_comparison_perform<FixedPoint<3, 28, uint16_t>>();
    test_comparison_perform<FixedPoint<3, 28, uint8_t>>();
  }

  template<typename T, int bits>
  void test_atan2_perform() {
    T x, y, a;

    for (unsigned i = 0; i < 1000; i++) {
      double angle = i * 2*M_PI / 1000.;
      double dy = sin(angle);
      double dx = cos(angle);
      double da = atan2(dy, dx);

      y = dy;
      x = dx;
      a = std::atan2(y, x);

      CHECK_CLOSE(da, double(a), 0.5*1e-4 * 2*M_PI);
    }
  }

  template<typename T, int bits>
  void test_sin_perform() {
    for (int i = -1000; i < 2000; i++) {
      double angle = i * 2*M_PI / 1000.;
      T a = angle;

      double floats = sin(angle);
      T b = std::sin(a);

      CHECK_CLOSE(floats, double(b), 1.8e-4 * 2*M_PI);
    }
  }

  template<typename T, int bits>
  void test_cos_perform() {
      for (int i = -1000; i < 2000; i++) {
        double angle = i * 2*M_PI / 1000.;
        T a = angle;

        double floats = cos(angle);
        T b = std::cos(a);

        CHECK_CLOSE(floats, double(b), 1.8e-4 * 2*M_PI);
      }
  }

  TEST(atan2) {
    test_atan2_perform<FixedPoint<4,91,uint32_t>, 32>();
    test_atan2_perform<FixedPoint<4,91,uint16_t>, 16>();
    test_atan2_perform<FixedPoint<4,91,uint8_t>, 8>();
  }

  TEST(sin) {
    test_sin_perform<FixedPoint<4,91,uint32_t>, 32>();
    test_sin_perform<FixedPoint<4,91,uint16_t>, 16>();
    test_sin_perform<FixedPoint<4,91,uint8_t>, 8>();
  }

  TEST(cos) {
    test_cos_perform<FixedPoint<4,91,uint32_t>, 32>();
    test_cos_perform<FixedPoint<4,91,uint16_t>, 16>();
    test_cos_perform<FixedPoint<4,91,uint8_t>, 8>();
  }

  TEST(abs) {
    FixedPoint<4, 91, uint32_t> a, b;
    FixedPoint<4, 91, uint16_t> c, d;
    FixedPoint<4, 91, uint8_t> e, f;

    a = c = e = -1.;
    b = d = f = 1.;

    CHECK_EQUAL(double(std::abs(a)), double(b));
    CHECK_EQUAL(double(std::abs(c)), double(d));
    CHECK_EQUAL(double(std::abs(e)), double(f));
    CHECK_EQUAL(double(std::abs(b)), double(b));
    CHECK_EQUAL(double(std::abs(d)), double(d));
    CHECK_EQUAL(double(std::abs(f)), double(f));

    a = c = e = 0.;
    CHECK_EQUAL(double(std::abs(a)), double(a));
    CHECK_EQUAL(double(std::abs(c)), double(c));
    CHECK_EQUAL(double(std::abs(e)), double(e));

    a = FixedPoint<4, 91, uint32_t>::maxVal;
    CHECK_EQUAL(double(std::abs(a)), double(a));
    CHECK_EQUAL(double(std::abs(-a)), double(a));

    a = FixedPoint<4, 91, uint32_t>::minVal;
    b = FixedPoint<4, 91, uint32_t>::maxVal;
    CHECK_EQUAL(double(std::abs(a)), double(b));
    CHECK_EQUAL(double(std::abs(b)), double(b));
  }

  TEST(ln2) {
    using Type = FixedPoint<3, 60, uint32_t>;
    Type b = std::ln2<3,60,uint32_t>();
    CHECK_CLOSE(std::log(2), double(b), 1e-14);
  } 
  
  TEST(log) {
    using Type = FixedPoint<3, 28+32, uint32_t>;
    CHECK_EQUAL(double(std::log(Type(-1))), double(Type::minVal));
    CHECK_CLOSE(double(std::log(Type(std::exp(1)))), 1, 1e-14);
    CHECK_CLOSE(double(std::log(Type(std::exp(2)))), 2, 1e-14);

    CHECK_CLOSE(200, double(std::log(FixedPoint<319, 64>(std::exp(200)))), 1e-14);
  } 
  
  TEST(exp) {
    using Type = FixedPoint<319, 64, uint32_t>;
    Type b = 5;
    CHECK_CLOSE(double(std::exp(b)), std::exp(double(b)), 1e-14);
  }

  TEST(pow) {
    using Type = FixedPoint<5, 26+32>;
    CHECK_CLOSE(double(std::pow(Type(5), -1)), 1./5, 1e-14);
    CHECK_CLOSE(double(std::pow(Type(5), -1)), 1./5, 1e-14);
    CHECK_CLOSE(double(std::pow(Type(25), Type(0.5))), 5, 1e-14);
  }

  TEST(raw_data) {
      FixedPoint<4, 11, uint16_t> a = 5.;
      uint8_t data[2];

      a.get_raw(data);
      CHECK_EQUAL(data[0], 0);
      CHECK_EQUAL(data[1], 5 << 3);
  }
}

int main()
{
  return UnitTest::RunAllTests();
}
