#include <iostream>
#include <stdint.h>
#include <UnitTest++.h>
#include "FixedPoint.hpp"
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
    MultiwordInteger<4, uint16_t> a(int64_t(1)<<63);
    CHECK_EQUAL(double(1L<<63), a.operator double());
  }

  TEST(storagetype_construction) {
    MultiwordInteger<2, uint16_t> a;

    // Positive values
    a = uint16_t(5);
    CHECK(a.operator double() == 5.);

    // Negative values get sign-extended
    a = uint16_t(-5);
    CHECK(a.operator double() == -5.);
  }

  TEST(double_construction) {
    MultiwordInteger<4, uint16_t> a;
    for (unsigned i = 0; i < 63; i++) {
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
    CHECK(double(a) == double(b));

    // Assignment from a type that doesn't fill the storage completely
    MultiwordInteger<3, uint8_t> c(int64_t(0xC0FFEE));
    a = c;
    CHECK(double(a) == double(c));
  }

  TEST(limits) {
    MultiwordInteger<2, uint16_t> maxv(MultiwordInteger<2, uint16_t>::maxVal());
    MultiwordInteger<2, uint16_t> minv(MultiwordInteger<2, uint16_t>::minVal());
    MultiwordInteger<2, uint16_t> a(int64_t(0x80000000));

    CHECK(a == minv);

    a=int64_t(0x7FFFFFFF);
    CHECK(a == maxv);
  }

  TEST(addition) {
    MultiwordInteger<2, uint16_t> a(0.), b(1.), c(-1.), d(2.);

    // Adding zero
    CHECK(a+b == b);

    // Basic addition
    CHECK(b+b == d);

    // Adding with overflow
    CHECK(a+b+c == a);

    // Commutativity
    CHECK(b+d == d+b);
  }

  TEST(subtraction) {
    MultiwordInteger<2, uint16_t> a(0.), b(1.), c(-1.), d(3.);

    // Subtracting from zero should yield negative
    CHECK(a-b == c);

    // Subtracting itself should yield zero
    CHECK(b-b == a);

    // Basic subtraction
    CHECK(d-b-b == b);

    // Anticommutativity
    CHECK(a-b == a-(b-a));
  }

  TEST(negation) {
    MultiwordInteger<2, uint16_t> a(1.), b(-1.), c(-1.);

    CHECK(a == -b);

    c.negate();
    CHECK(a == c);
  }

  TEST(multiplication) {
    MultiwordInteger<2, uint16_t> a(0.), b(1.), c(-1.), d(2.);

    // Multiplying with zero yields zero
    CHECK(a*d == a);

    // Multiplying with 1 yields the other operand
    CHECK(b*d == d);

    // Multiplying with -1 yields -d
    CHECK(d*c == -d);

    // Commutativity
    CHECK(c*d == d*c);

    // Results out of range of the smaller type should still work in the resulting type
    b = double(0x10000);
    MultiwordInteger<4, uint16_t> e(double(0x100000000));
    CHECK(b*b == e);

    // But overflow when assigned to a smaller type
    MultiwordInteger<2, uint16_t> f(b);
    f *= b;
    CHECK(f == a);
  }

  TEST(division) {
    MultiwordInteger<2, uint16_t> a(1.), b(2.), c, d, zero(0.), unity(1.);

    // Dividing by larger values yields zero for positive values, -1 for negative values
    CHECK(a/b == zero);
    CHECK(a/(-b) == -unity);

    // Dividing by 1 yields the numerator
    CHECK(b/unity == b);

    // x/x == 1
    CHECK(b/b == unity);

    // Test again, for larger values
    a = 65536.*5;
    b = 65536.;
    CHECK(a/a == unity);

    // If there is no remainder, a/b == c <-> a/c == b
    c = 5.;
    CHECK(a/b == c);
    CHECK(a/c == b);
    CHECK((-a)/b == -c);
    CHECK(a/(-b) == -c);
    CHECK((-a)/(-b) == c);

    // Test with remainder
    a = 5.;
    b = 2.;
    c = 3.;
    CHECK(a/-b == -c);
    CHECK(a/-c == -b);
    CHECK(-a/b == -c);
    CHECK(-a/c == -b);

    // And for large values
    a = 5*65536.;
    b = 3*65536.;
    c = 1.;
    d = 2.;
    CHECK(a/-b == -d);
    CHECK((-a)/b == -d);
    CHECK((-a)/(-b) == c);

    // Division by zero
    a = 1.;
    b = 0.;
    c = MultiwordInteger<2, uint16_t>::maxVal();
    d = MultiwordInteger<2, uint16_t>::minVal();
    CHECK(a/b == c);
    CHECK((-a)/b == d);
  }

  TEST(modulo) {
    MultiwordInteger<2, uint16_t> a, b, c, d, e;

    // Without remainder, positive and negative
    a = 5.;
    b = 1.;
    c = 0.;
    d = a/b * b;
    e = a%b;
    CHECK(e == c);
    CHECK(d + e == a);

    d = a / (-b) * (-b);
    e = a % (-b);
    CHECK(e == c);
    CHECK(d + e == a);

    d = (-a) / b * b;
    e = (-a) % b;
    CHECK(e == c);
    CHECK(d + e == -a);

    d = (-a) / (-b) * (-b);
    e = (-a) % (-b);
    CHECK(e == -c);
    CHECK(d + e == -a);

    // With remainder
    a = 7.;
    b = 3.;
    c = 1.;

    d = a/b*b;
    e = a%b;
    CHECK(e == c);
    CHECK(d + e == a);

    d = (-a) / (-b) * (-b);
    e = (-a) % (-b);
    CHECK(e == -c);
    CHECK(d + e == (-a));

    c = 2.;
    d = a/(-b)*(-b);
    e = a%(-b);
    CHECK(e == -c);
    CHECK(d + e == a);

    d = (-a)/b*b;
    e = (-a)%b;
    CHECK(e == c);
    CHECK(d + e == -a);

    // Large divisor
    a = 65536.*2 - 1;
    b = 65536.;

    c = 65535.;
    d = a/b*b;
    e = a%b;
    CHECK(e == c);
    CHECK(d + e == a);

    d = (-a) / (-b) * (-b);
    e = (-a) % (-b);
    CHECK(e == -c);
    CHECK(d + e == -a);

    c = 1.;
    d = (-a) / b * b;
    e = (-a) % b;
    CHECK(e == c);
    CHECK(d + e == -a);

    d = a / (-b) * (-b);
    e = a % (-b);
    CHECK(e == -c);
    CHECK(d + e == a);

    // Division of zero returns 0
    a = 0.;
    b = 1.;
    d = a / b * b;
    e = a % b;
    CHECK(d == a);
    CHECK(e == a);

    // Division by zero returns 0
    a = 1.;
    b = 0.;
    e = a % b;
    CHECK(e == b);
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
    CHECK(a == a);

    CHECK(!(a == b));
    CHECK(!(a != a));

    b = -1.;
    CHECK(b < a);
  }

  TEST(shifts) {
    MultiwordInteger<2, uint16_t> a = 65536., b = 0., c;

    // Test shifting all bits out, when positive
    a >>= 17;
    CHECK(a == b);

    // Test shifting 0
    a <<= 1;
    CHECK(a == b);

    // Huge shift
    a = 1.;
    a <<= 100;
    CHECK(a == b);

    a = 1.;
    // Test normal shift
    a <<= 1;
    b = 2.;
    CHECK(a == b);

    // Shift over boundary
    a <<= 17;
    b = 65536.*4;
    CHECK(a == b);

    a <<= 12;

    // Test right shift
    a >>= 1;
    b = 536870912.;
    CHECK(a == b);

    // Test large right shift
    a >>= 29;
    b = 1.;
    CHECK(a == b);

    // Result should be negative
    a <<= 31;
    b = MultiwordInteger<2, uint16_t>::minVal();
    CHECK(a == b);

    // Result should be sign extended, and negative
    a >>= 30;
    b = -2.;
    CHECK(a == b);

    // Too large shift should still sign extend
    a >>= 300;
    b = -1.;
    CHECK(a == b);

    // Huge shift for positive values
    a = 555.;
    a >>= 300;
    b = 0.;
    CHECK(a == b);

    // Shift by storageSize
    a = 65536.;
    a >>= 16;
    b = 1.;
    CHECK(a == b);

    // Shift by storageSize, negative number
    a = -65536.;
    a >>= 16;
    b = -1.;
    CHECK(a == b);
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
    b = a; c = a; d = a; e = a; f = a; g = a; h = a; i = a; j = a;
    CHECK(double(b) == double(a));
    CHECK(double(c) == double(a));
    CHECK(double(d) == double(a));
    CHECK(double(e) == double(a));
    CHECK(double(f) == double(a));
    CHECK(double(g) == double(a));
    CHECK(double(h) == double(a));
    CHECK(double(i) == double(a));
    CHECK(double(j) == double(a));

    // Larger value conversion
    a = 65536.;
    b = a; c = a; d = a; e = a; f = a; g = a; h = a; i = a; j = a;
    CHECK(double(b) == double(a));
    CHECK(double(c) == double(a));
    CHECK(double(d) == double(a));
    CHECK(double(e) == double(a));
    CHECK(double(f) == double(a));
    CHECK(double(g) == double(a));
    CHECK(double(h) == double(a));
    CHECK(double(i) == double(a));
    CHECK(double(j) == double(a));

    // Negative conversion
    a = -1.;
    b = a; c = a; d = a; e = a; f = a; g = a; h = a; i = a; j = a;
    CHECK(double(b) == double(a));
    CHECK(double(c) == double(a));
    CHECK(double(d) == double(a));
    CHECK(double(e) == double(a));
    CHECK(double(f) == double(a));
    CHECK(double(g) == double(a));
    CHECK(double(h) == double(a));
    CHECK(double(i) == double(a));
    CHECK(double(j) == double(a));

    // Overflow
    a = pow(2., 32.);
    b = a; e = a; h = a;
    CHECK(!b);
    CHECK(!e);
    CHECK(!h);
  }
}

SUITE(FixedPoint) {
  TEST(construction) {
    FixedPoint<3, 15, uint16_t> a(5.), b;
    FixedPoint<5, 4, uint32_t> c;
    FixedPoint<5, 18, uint32_t> d;

    FixedPoint<31, 0>::StorageType s(int64_t(1));

    CHECK(double(a) == 5.);

    a = int(5);
    CHECK(double(a) == 5.);

    b = a;
    CHECK(double(b) == 5.);

    c = a;
    CHECK(double(c) == 5.);

    d = a;
    CHECK(double(d) == 5.);

    FixedPoint<31, 0> e(s);
    CHECK(double(e) == 1.);

    FixedPoint<30, 1> f(s);
    CHECK(double(f) == 0.5);

    a = 100.;
    CHECK(double(a) == (FixedPoint<3, 15, uint16_t>::maxVal<double>()));

    a = -100.;
    CHECK(double(a) == (FixedPoint<3, 15, uint16_t>::minVal<double>()));

    a = int(100);
    CHECK(double(a) == (FixedPoint<3, 15, uint16_t>::maxVal<double>()));

    a = int(-100);
    CHECK(double(a) == (FixedPoint<3, 15, uint16_t>::minVal<double>()));
  }

  TEST(limits) {
    FixedPoint<5, 10, uint16_t> a;
    double max = pow(2., 5) - pow(2., -(double)a.fractionalWidth());
    double min = -pow(2., 5.);
    double smallest = pow(2., -(double)a.fractionalWidth());

    a = FixedPoint<5, 10, uint16_t>::template maxVal<FixedPoint<5,10,uint16_t>>();
    assert(double(a) == max);

    a = FixedPoint<5, 10, uint16_t>::template minVal<FixedPoint<5,10,uint16_t>>();
    assert(double(a) == min);

    a = FixedPoint<5, 10, uint16_t>::template smallestVal<FixedPoint<5,10,uint16_t>>();
    assert(double(a) == smallest);
  }

  TEST(addition) {
    FixedPoint<5, 10, uint16_t> a, b, c;

    a = 1.;
    b = pow(2., -10.);
    c = (1. + pow(2., -10));
    CHECK((a + b) == c);
    CHECK((b + a) == c);

    a = -1.;
    b = -2.;
    c = -3.;
    CHECK((a+b) == c);
  }

  TEST(subtraction) {
    FixedPoint<5, 10, uint16_t> a, b, c;

    a = 1.;
    b = pow(2., -10.);
    c = (1. - pow(2., -10));
    CHECK((a - b) == c);
    CHECK((b - a) == -c);
  }

  TEST(multiplication) {
    FixedPoint<5, 10, uint16_t> a, b, c;

    a = 1.;
    b = pow(2., -10.);
    CHECK((a * b) == b);

    a = 5.;
    b = 0.;
    CHECK(a * b == b);

    b = 2.;
    c = 10.;
    CHECK(a * b == c);
    CHECK(b * a == c);
  }

  TEST(division) {
    FixedPoint<5, 10, uint16_t> a, b, c;

    a = 1.;
    b = pow(2., -4.);
    c = pow(2., 4.);
    CHECK((a / b) == c);
    CHECK((a / c) == b);
    CHECK((a / (-b)) == -c);
    CHECK((a / (-c)) == -b);
    CHECK(((-a) / (-b)) == c);
    CHECK(((-a) / (-c)) == b);

    FixedPoint<4, 91, uint32_t> y32(-0.13779029068463858), x32(-0.99046142569665119);
    FixedPoint<4, 91, uint32_t> d32 = y32/x32;

    FixedPoint<4, 91, uint16_t> y16(y32), x16(x32), d16(y16/x16);
    FixedPoint<4, 91, uint8_t> y8(y32), x8(x32), d8(y16/x16);
    FixedPoint<4, 91> d8_32(d8), d16_32(d16);

    CHECK(double(y32) == double(y16));
    CHECK(double(x32) == double(x16));
    CHECK(d8_32 == d32);
    CHECK(d16_32 == d32);
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
    CHECK(a == a);

    CHECK(!(a == b));
    CHECK(!(a != a));

    b = -1.;
    a = -0.5;
    CHECK(b < a);
    CHECK(a > b);
    CHECK(b <= a);
    CHECK(a >= b);
  }

  TEST(comparison) {
    test_comparison_perform<FixedPoint<3, 4, uint32_t>>();
    test_comparison_perform<FixedPoint<3, 4, uint16_t>>();
    test_comparison_perform<FixedPoint<3, 4, uint8_t>>();
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

  TEST(atan2) {
    test_atan2_perform<FixedPoint<4,80,uint32_t>, 32>();
    test_atan2_perform<FixedPoint<4,80,uint16_t>, 16>();
    test_atan2_perform<FixedPoint<4,80,uint8_t>, 8>();
  }

  TEST(abs) {
    FixedPoint<4, 91, uint32_t> a, b;
    FixedPoint<4, 91, uint16_t> c, d;
    FixedPoint<4, 91, uint8_t> e, f;

    a = c = e = -1.;
    b = d = f = 1.;

    CHECK(std::abs(a) == b);
    CHECK(std::abs(c) == d);
    CHECK(std::abs(e) == f);
    CHECK(std::abs(b) == b);
    CHECK(std::abs(d) == d);
    CHECK(std::abs(f) == f);

    a = c = e = 0.;
    CHECK(std::abs(a) == a);
    CHECK(std::abs(c) == c);
    CHECK(std::abs(e) == e);

    a = FixedPoint<4, 91, uint32_t>::template maxVal<FixedPoint<4, 91, uint32_t>>();
    CHECK(std::abs(a) == a);
    CHECK(std::abs(-a) == a);

    a = FixedPoint<4, 91, uint32_t>::template minVal<FixedPoint<4, 91, uint32_t>>();
    b = FixedPoint<4, 91, uint32_t>::template maxVal<FixedPoint<4, 91, uint32_t>>();
    CHECK(std::abs(a) == b);
    CHECK(std::abs(b) == b);
  }
}

int main()
{
  return UnitTest::RunAllTests();
}
