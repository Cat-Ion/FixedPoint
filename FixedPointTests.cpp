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
    CHECK_EQUAL(double(int64_t(1)<<63), a.operator double());
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
    CHECK_EQUAL(double(a), double(b));

    // Assignment from a type that doesn't fill the storage completely
    MultiwordInteger<3, uint8_t> c(int64_t(0xC0FFEE));
    a = c;
    CHECK_EQUAL(double(a), double(c));
  }

  TEST(limits) {
    MultiwordInteger<2, uint16_t> maxv(MultiwordInteger<2, uint16_t>::maxVal);
    MultiwordInteger<2, uint16_t> minv(MultiwordInteger<2, uint16_t>::minVal);
    MultiwordInteger<2, uint16_t> a(int64_t(0x80000000));

    CHECK_EQUAL(double(a), double(minv));

    a=int64_t(0x7FFFFFFF);
    CHECK_EQUAL(double(a), double(maxv));
  }

  TEST(addition) {
    MultiwordInteger<2, uint16_t> a(0.), b(1.), c(-1.), d(2.);

    // Adding zero
    CHECK_EQUAL(double(a+b), double(b));

    // Basic addition
    CHECK_EQUAL(double(b+b), double(d));

    // Adding with overflow
    CHECK_EQUAL(double(a+b+c), double(a));

    // Commutativity
    CHECK_EQUAL(double(b+d), double(d+b));
  }

  TEST(subtraction) {
    MultiwordInteger<2, uint16_t> a(0.), b(1.), c(-1.), d(3.);

    // Subtracting from zero should yield negative
    CHECK_EQUAL(double(a-b), double(c));

    // Subtracting itself should yield zero
    CHECK_EQUAL(double(b-b), double(a));

    // Basic subtraction
    CHECK_EQUAL(double(d-b-b), double(b));

    // Anticommutativity
    CHECK_EQUAL(double(a-b), double(a-(b-a)));
  }

  TEST(negation) {
    MultiwordInteger<2, uint16_t> a(1.), b(-1.), c(-1.);

    CHECK_EQUAL(double(a), double(-b));

    c.negate();
    CHECK_EQUAL(double(a), double(c));
  }

  TEST(multiplication) {
    MultiwordInteger<2, uint16_t> a(0.), b(1.), c(-1.), d(2.);

    // Multiplying with zero yields zero
    CHECK_EQUAL(double(a*d), double(a));

    // Multiplying with 1 yields the other operand
    CHECK_EQUAL(double(b*d), double(d));

    // Multiplying with -1 yields -d
    CHECK_EQUAL(double(d*c), double(-d));

    // Commutativity
    CHECK_EQUAL(double(c*d), double(d*c));

    // Results out of range of the smaller type should still work in the resulting type
    b = double(0x10000);
    MultiwordInteger<4, uint16_t> e(double(0x100000000));
    CHECK_EQUAL(double(b*b), double(e));

    // But overflow when assigned to a smaller type
    MultiwordInteger<2, uint16_t> f(b);
    f *= b;
    CHECK_EQUAL(double(f), double(a));
  }

  TEST(division) {
    MultiwordInteger<2, uint16_t> a(1.), b(2.), c, d, zero(0.), unity(1.);

    // Dividing by larger values yields zero for positive values, -1 for negative values
    CHECK_EQUAL(double(a/b), double(zero));
    CHECK_EQUAL(double(a/(-b)), double(-unity));

    // Dividing by 1 yields the numerator
    CHECK_EQUAL(double(b/unity), double(b));

    // x/x == 1
    CHECK_EQUAL(double(b/b), double(unity));

    // Test again, for larger values
    a = 65536.*5;
    b = 65536.;
    CHECK_EQUAL(double(a/a), double(unity));

    // If there is no remainder, a/b == c <-> a/c == b
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

    // And for large values
    a = 5*65536.;
    b = 3*65536.;
    c = 1.;
    d = 2.;
    CHECK_EQUAL(double(a/-b), double(-d));
    CHECK_EQUAL(double((-a)/b), double(-d));
    CHECK_EQUAL(double((-a)/(-b)), double(c));

    // Division by zero
    a = 1.;
    b = 0.;
    c = MultiwordInteger<2, uint16_t>::maxVal;
    d = MultiwordInteger<2, uint16_t>::minVal;
    CHECK_EQUAL(double(a/b), double(c));
    CHECK_EQUAL(double((-a)/b), double(d));
  }

  TEST(modulo) {
    MultiwordInteger<2, uint16_t> a, b, c, d, e;

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

  TEST(shifts) {
    MultiwordInteger<2, uint16_t> a = 65536., b = 0., c;

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
    a <<= 1;
    b = 2.;
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
    b = MultiwordInteger<2, uint16_t>::minVal;
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

  TEST(bitwise) {
      MultiwordInteger<4, uint8_t> a, b, c;
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
    CHECK_EQUAL(double(double(b)), double(double(a)));
    CHECK_EQUAL(double(double(c)), double(double(a)));
    CHECK_EQUAL(double(double(d)), double(double(a)));
    CHECK_EQUAL(double(double(e)), double(double(a)));
    CHECK_EQUAL(double(double(f)), double(double(a)));
    CHECK_EQUAL(double(double(g)), double(double(a)));
    CHECK_EQUAL(double(double(h)), double(double(a)));
    CHECK_EQUAL(double(double(i)), double(double(a)));
    CHECK_EQUAL(double(double(j)), double(double(a)));

    // Larger value conversion
    a = 65536.;
    b = a; c = a; d = a; e = a; f = a; g = a; h = a; i = a; j = a;
    CHECK_EQUAL(double(double(b)), double(double(a)));
    CHECK_EQUAL(double(double(c)), double(double(a)));
    CHECK_EQUAL(double(double(d)), double(double(a)));
    CHECK_EQUAL(double(double(e)), double(double(a)));
    CHECK_EQUAL(double(double(f)), double(double(a)));
    CHECK_EQUAL(double(double(g)), double(double(a)));
    CHECK_EQUAL(double(double(h)), double(double(a)));
    CHECK_EQUAL(double(double(i)), double(double(a)));
    CHECK_EQUAL(double(double(j)), double(double(a)));

    // Negative conversion
    a = -1.;
    b = a; c = a; d = a; e = a; f = a; g = a; h = a; i = a; j = a;
    CHECK_EQUAL(double(double(b)), double(double(a)));
    CHECK_EQUAL(double(double(c)), double(double(a)));
    CHECK_EQUAL(double(double(d)), double(double(a)));
    CHECK_EQUAL(double(double(e)), double(double(a)));
    CHECK_EQUAL(double(double(f)), double(double(a)));
    CHECK_EQUAL(double(double(g)), double(double(a)));
    CHECK_EQUAL(double(double(h)), double(double(a)));
    CHECK_EQUAL(double(double(i)), double(double(a)));
    CHECK_EQUAL(double(double(j)), double(double(a)));

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
    FixedPoint<3, 28, uint16_t> a(5.), b;
    FixedPoint<5, 26, uint32_t> c;
    FixedPoint<5, 26, uint32_t> d;

    FixedPoint<31, 0>::StorageType s(int64_t(1));

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
}

int main()
{
  return UnitTest::RunAllTests();
}
