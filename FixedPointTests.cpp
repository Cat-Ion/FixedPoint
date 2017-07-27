#include <iostream>
#include <cpptest.h>
#include <stdint.h>
#include "FixedPoint.hpp"

class MultiwordIntegerTestSuite : public Test::Suite
{
public:
    MultiwordIntegerTestSuite()
    {
        TEST_ADD(MultiwordIntegerTestSuite::test_zero_construction);
        TEST_ADD(MultiwordIntegerTestSuite::test_lsb_construction);
        TEST_ADD(MultiwordIntegerTestSuite::test_msb_construction);
        TEST_ADD(MultiwordIntegerTestSuite::test_storagetype_construction);
        TEST_ADD(MultiwordIntegerTestSuite::test_double_construction);
        TEST_ADD(MultiwordIntegerTestSuite::test_limits);
        TEST_ADD(MultiwordIntegerTestSuite::test_addition);
        TEST_ADD(MultiwordIntegerTestSuite::test_subtraction);
        TEST_ADD(MultiwordIntegerTestSuite::test_negation);
        TEST_ADD(MultiwordIntegerTestSuite::test_multiplication);
        TEST_ADD(MultiwordIntegerTestSuite::test_division);
        TEST_ADD(MultiwordIntegerTestSuite::test_modulo);
        TEST_ADD(MultiwordIntegerTestSuite::test_comparison);
        TEST_ADD(MultiwordIntegerTestSuite::test_shifts);
        TEST_ADD(MultiwordIntegerTestSuite::test_conversions);
    }
protected:
    void test_zero_construction() {
        MultiwordInteger<4, uint16_t> a(int64_t(0));
        TEST_ASSERT_EQUALS(double(0), a.operator double());
    }

    void test_lsb_construction() {
        MultiwordInteger<4, uint16_t> a(int64_t(1));
        TEST_ASSERT_EQUALS(double(1), a.operator double());
    }

    void test_msb_construction() {
        MultiwordInteger<4, uint16_t> a(int64_t(1)<<63);
        TEST_ASSERT_EQUALS(double(1L<<63), a.operator double());
    }

    void test_storagetype_construction() {
        MultiwordInteger<2, uint16_t> a;

        // Positive values
        a = uint16_t(5);
        TEST_ASSERT(a.operator double() == 5.);

        // Negative values get sign-extended
        a = uint16_t(-5);
        TEST_ASSERT(a.operator double() == -5.);
    }

    void test_double_construction() {
        MultiwordInteger<4, uint16_t> a;
        for (unsigned i = 0; i < 63; i++) {
            double b = pow(2., i);
            a = b;
            TEST_ASSERT_EQUALS(b, a.operator double());
        }
    }

    void test_other_construction() {
        MultiwordInteger<4, uint16_t> a(int64_t(0x7FFF7FFF7FFF7FFF));
        MultiwordInteger<2, uint16_t> b(a);
        // Assignment should truncate
        TEST_ASSERT_EQUALS(b.operator double(), double(0x7FFF7FFF));

        a = int64_t(0x7FFF7FFF80000000);
        b = a;
        // This should be negative
        TEST_ASSERT(b.is_negative());
        TEST_ASSERT_EQUALS(b.operator double(), double(0x80000000));

        b = int64_t(0x7FFFFFFF);
        a = b;
        // Assignment from a smaller type should be exact
        TEST_ASSERT(a == b);
    }

    void test_limits() {
        MultiwordInteger<2, uint16_t> maxv(MultiwordInteger<2, uint16_t>::maxVal());
        MultiwordInteger<2, uint16_t> minv(MultiwordInteger<2, uint16_t>::minVal());
        MultiwordInteger<2, uint16_t> a(int64_t(0x80000000));

        TEST_ASSERT(a == minv);

        a=int64_t(0x7FFFFFFF);
        TEST_ASSERT(a == maxv);
    }

    void test_addition() {
        MultiwordInteger<2, uint16_t> a(0.), b(1.), c(-1.), d(2.);

        // Adding zero
        TEST_ASSERT(a+b == b);

        // Basic addition
        TEST_ASSERT(b+b == d);

        // Adding with overflow
        TEST_ASSERT(a+b+c == a);

        // Commutativity
        TEST_ASSERT(b+d == d+b);
    }

    void test_subtraction() {
        MultiwordInteger<2, uint16_t> a(0.), b(1.), c(-1.), d(3.);

        // Subtracting from zero should yield negative
        TEST_ASSERT(a-b == c);

        // Subtracting itself should yield zero
        TEST_ASSERT(b-b == a);

        // Basic subtraction
        TEST_ASSERT(d-b-b == b);

        // Anticommutativity
        TEST_ASSERT(a-b == a-(b-a));
    }

    void test_negation() {
        MultiwordInteger<2, uint16_t> a(1.), b(-1.), c(-1.);

        TEST_ASSERT(a == -b);

        c.negate();
        TEST_ASSERT(a == c);
    }

    void test_multiplication() {
        MultiwordInteger<2, uint16_t> a(0.), b(1.), c(-1.), d(2.);

        // Multiplying with zero yields zero
        TEST_ASSERT(a*d == a);

        // Multiplying with 1 yields the other operand
        TEST_ASSERT(b*d == d);

        // Multiplying with -1 yields -d
        TEST_ASSERT(d*c == -d);

        // Commutativity
        TEST_ASSERT(c*d == d*c);

        // Results out of range of the smaller type should still work in the resulting type
        b = double(0x10000);
        MultiwordInteger<4, uint16_t> e(double(0x100000000));
        TEST_ASSERT(b*b == e);

        // But overflow when assigned to a smaller type
        MultiwordInteger<2, uint16_t> f(b);
        f *= b;
        TEST_ASSERT(f == a);
    }

    void test_division() {
        MultiwordInteger<2, uint16_t> a(1.), b(2.), c, d, zero(0.), unity(1.);

        // Dividing by larger values yields zero for positive values, -1 for negative values
        TEST_ASSERT(a/b == zero);
        TEST_ASSERT(a/(-b) == -unity);

        // Dividing by 1 yields the numerator
        TEST_ASSERT(b/unity == b);

        // x/x == 1
        TEST_ASSERT(b/b == unity);

        // Test again, for larger values
        a = 65536.*5;
        b = 65536.;
        TEST_ASSERT(a/a == unity);

        // If there is no remainder, a/b == c <-> a/c == b
        c = 5.;
        TEST_ASSERT(a/b == c);
        TEST_ASSERT(a/c == b);
        TEST_ASSERT((-a)/b == -c);
        TEST_ASSERT(a/(-b) == -c);
        TEST_ASSERT((-a)/(-b) == c);

        // Test with remainder
        a = 5.;
        b = 2.;
        c = 3.;
        TEST_ASSERT(a/-b == -c);
        TEST_ASSERT(a/-c == -b);
        TEST_ASSERT(-a/b == -c);
        TEST_ASSERT(-a/c == -b);

        // And for large values
        a = 5*65536.;
        b = 3*65536.;
        c = 1.;
        d = 2.;
        TEST_ASSERT(a/-b == -d);
        TEST_ASSERT((-a)/b == -d);
        TEST_ASSERT((-a)/(-b) == c);
    }

    void test_modulo() {
        MultiwordInteger<2, uint16_t> a, b, c, d, e;

        // Without remainder, positive and negative
        a = 5.;
        b = 1.;
        c = 0.;
        d = a/b * b;
        e = a%b;
        TEST_ASSERT(e == c);
        TEST_ASSERT(d + e == a);

        d = a / (-b) * (-b);
        e = a % (-b);
        TEST_ASSERT(e == c);
        TEST_ASSERT(d + e == a);

        d = (-a) / b * b;
        e = (-a) % b;
        TEST_ASSERT(e == c);
        TEST_ASSERT(d + e == -a);

        d = (-a) / (-b) * (-b);
        e = (-a) % (-b);
        TEST_ASSERT(e == -c);
        TEST_ASSERT(d + e == -a);

        // With remainder
        a = 7.;
        b = 3.;
        c = 1.;

        d = a/b*b;
        e = a%b;
        TEST_ASSERT(e == c);
        TEST_ASSERT(d + e == a);

        d = (-a) / (-b) * (-b);
        e = (-a) % (-b);
        TEST_ASSERT(e == -c);
        TEST_ASSERT(d + e == (-a));

        c = 2.;
        d = a/(-b)*(-b);
        e = a%(-b);
        TEST_ASSERT(e == -c);
        TEST_ASSERT(d + e == a);

        d = (-a)/b*b;
        e = (-a)%b;
        TEST_ASSERT(e == c);
        TEST_ASSERT(d + e == -a);

        // Large divisor
        a = 65536.*2 - 1;
        b = 65536.;

        c = 65535.;
        d = a/b*b;
        e = a%b;
        TEST_ASSERT(e == c);
        TEST_ASSERT(d + e == a);

        d = (-a) / (-b) * (-b);
        e = (-a) % (-b);
        TEST_ASSERT(e == -c);
        TEST_ASSERT(d + e == -a);

        c = 1.;
        d = (-a) / b * b;
        e = (-a) % b;
        TEST_ASSERT(e == c);
        TEST_ASSERT(d + e == -a);

        d = a / (-b) * (-b);
        e = a % (-b);
        TEST_ASSERT(e == -c);
        TEST_ASSERT(d + e == a);
    }

    void test_comparison() {
        MultiwordInteger<2, uint16_t> a, b;

        a = 1.;
        b = 0.;

        TEST_ASSERT(a > b);
        TEST_ASSERT(a >= b);
        TEST_ASSERT(b < a);
        TEST_ASSERT(b <= a);

        TEST_ASSERT(a != b);
        TEST_ASSERT(a == a);

        TEST_ASSERT(!(a == b));
        TEST_ASSERT(!(a != a));

        b = -1.;
        TEST_ASSERT(b < a);
    }

    void test_shifts() {
        MultiwordInteger<2, uint16_t> a = 65536., b = 0., c;

        // Test shifting all bits out, when positive
        a >>= 17;
        TEST_ASSERT(a == b);

        // Test shifting 0
        a <<= 1;
        TEST_ASSERT(a == b);

        // Huge shift
        a = 1.;
        a <<= 100;
        TEST_ASSERT(a == b);

        a = 1.;
        // Test normal shift
        a <<= 1;
        b = 2.;
        TEST_ASSERT(a == b);

        // Shift over boundary
        a <<= 17;
        b = 65536.*4;
        TEST_ASSERT(a == b);

        a <<= 12;

        // Test right shift
        a >>= 1;
        b = 536870912.;
        TEST_ASSERT(a == b);

        // Test large right shift
        a >>= 29;
        b = 1.;
        TEST_ASSERT(a == b);

        // Result should be negative
        a <<= 31;
        b = MultiwordInteger<2, uint16_t>::minVal();
        TEST_ASSERT(a == b);

        // Result should be sign extended, and negative
        a >>= 30;
        b = -2.;
        TEST_ASSERT(a == b);

        // Too large shift should still sign extend
        a >>= 300;
        b = -1.;
        TEST_ASSERT(a == b);
    }

    void test_conversions() {
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
        TEST_ASSERT(double(b) == double(a));
        TEST_ASSERT(double(c) == double(a));
        TEST_ASSERT(double(d) == double(a));
        TEST_ASSERT(double(e) == double(a));
        TEST_ASSERT(double(f) == double(a));
        TEST_ASSERT(double(g) == double(a));
        TEST_ASSERT(double(h) == double(a));
        TEST_ASSERT(double(i) == double(a));
        TEST_ASSERT(double(j) == double(a));

        // Larger value conversion
        a = 65536.;
        b = a; c = a; d = a; e = a; f = a; g = a; h = a; i = a; j = a;
        TEST_ASSERT(double(b) == double(a));
        TEST_ASSERT(double(c) == double(a));
        TEST_ASSERT(double(d) == double(a));
        TEST_ASSERT(double(e) == double(a));
        TEST_ASSERT(double(f) == double(a));
        TEST_ASSERT(double(g) == double(a));
        TEST_ASSERT(double(h) == double(a));
        TEST_ASSERT(double(i) == double(a));
        TEST_ASSERT(double(j) == double(a));

        // Negative conversion
        a = -1.;
        b = a; c = a; d = a; e = a; f = a; g = a; h = a; i = a; j = a;
        TEST_ASSERT(double(b) == double(a));
        TEST_ASSERT(double(c) == double(a));
        TEST_ASSERT(double(d) == double(a));
        TEST_ASSERT(double(e) == double(a));
        TEST_ASSERT(double(f) == double(a));
        TEST_ASSERT(double(g) == double(a));
        TEST_ASSERT(double(h) == double(a));
        TEST_ASSERT(double(i) == double(a));
        TEST_ASSERT(double(j) == double(a));

        // Overflow
        a = pow(2., 32.);
        b = a; e = a; h = a;
        TEST_ASSERT(!b);
        TEST_ASSERT(!e);
        TEST_ASSERT(!h);
    }
};

class FixedPointTestSuite : public Test::Suite {
public:
    FixedPointTestSuite() {
        TEST_ADD(FixedPointTestSuite::test_construction);
        TEST_ADD(FixedPointTestSuite::test_limits);
        TEST_ADD(FixedPointTestSuite::test_addition);
        TEST_ADD(FixedPointTestSuite::test_subtraction);
        TEST_ADD(FixedPointTestSuite::test_multiplication);
        TEST_ADD(FixedPointTestSuite::test_division);
        TEST_ADD(FixedPointTestSuite::test_comparison);
        TEST_ADD(FixedPointTestSuite::test_atan2);
        TEST_ADD(FixedPointTestSuite::test_abs);
    }
protected:
    void test_construction() {
        FixedPoint<3, 15, uint16_t> a(5.), b;
        FixedPoint<5, 4, uint32_t> c;
        FixedPoint<5, 18, uint32_t> d;

        FixedPoint<31, 0>::StorageType s(int64_t(1));

        TEST_ASSERT(double(a) == 5.);

        a = int(5);
        TEST_ASSERT(double(a) == 5.);

        b = a;
        TEST_ASSERT(double(b) == 5.);

        c = a;
        TEST_ASSERT(double(c) == 5.);

        d = a;
        TEST_ASSERT(double(d) == 5.);

        FixedPoint<31, 0> e(s);
        TEST_ASSERT(double(e) == 1.);

        FixedPoint<30, 1> f(s);
        TEST_ASSERT(double(f) == 0.5);
    }

    void test_limits() {
        FixedPoint<5, 31> a;

        a = FixedPoint<5, 31>::template maxVal<FixedPoint<5,31>>();
        assert(double(a) == (pow(2, 5) - pow(2, -31)));

        a = FixedPoint<5, 31>::template minVal<FixedPoint<5,31>>();
        assert(double(a) == -pow(2, 5));

        a = FixedPoint<5, 31>::template smallestVal<FixedPoint<5,31>>();
        assert(double(a) == pow(2, -31));
    }

    void test_addition() {
        FixedPoint<5, 10, uint16_t> a, b, c;

        a = 1.;
        b = pow(2., -10.);
        c = (1. + pow(2., -10));
        TEST_ASSERT((a + b) == c);
        TEST_ASSERT((b + a) == c);
    }

    void test_subtraction() {
        FixedPoint<5, 10, uint16_t> a, b, c;

        a = 1.;
        b = pow(2., -10.);
        c = (1. - pow(2., -10));
        TEST_ASSERT((a - b) == c);
        TEST_ASSERT((b - a) == -c);
    }

    void test_multiplication() {
        FixedPoint<5, 10, uint16_t> a, b, c;

        a = 1.;
        b = pow(2., -10.);
        TEST_ASSERT((a * b) == b);

        a = 5.;
        b = 0.;
        TEST_ASSERT(a * b == b);

        b = 2.;
        c = 10.;
        TEST_ASSERT(a * b == c);
        TEST_ASSERT(b * a == c);
    }

    void test_division() {
        FixedPoint<5, 10, uint16_t> a, b, c;

        a = 1.;
        b = pow(2., -4.);
        c = pow(2., 4.);
        TEST_ASSERT((a / b) == c);
        TEST_ASSERT((a / c) == b);
        TEST_ASSERT((a / (-b)) == -c);
        TEST_ASSERT((a / (-c)) == -b);
        TEST_ASSERT(((-a) / (-b)) == c);
        TEST_ASSERT(((-a) / (-c)) == b);

        FixedPoint<4, 91, uint32_t> y32(-0.13779029068463858), x32(-0.99046142569665119);
        FixedPoint<4, 91, uint32_t> d32 = y32/x32;

        FixedPoint<4, 91, uint16_t> y16(y32), x16(x32), d16(y16/x16);
        FixedPoint<4, 91, uint8_t> y8(y32), x8(x32), d8(y16/x16);
        FixedPoint<4, 91> d8_32(d8), d16_32(d16);

        TEST_ASSERT(double(y32) == double(y16));
        TEST_ASSERT(double(x32) == double(x16));
        TEST_ASSERT(d8_32 == d32);
        TEST_ASSERT(d16_32 == d32);
    }

    void test_comparison() {
        test_comparison_perform<FixedPoint<3, 4, uint32_t>>();
        test_comparison_perform<FixedPoint<3, 4, uint16_t>>();
        test_comparison_perform<FixedPoint<3, 4, uint8_t>>();
    }

    template <typename T>
    void test_comparison_perform() {
        T a, b;

        a = 1.;
        b = 0.;

        TEST_ASSERT(a > b);
        TEST_ASSERT(a >= b);
        TEST_ASSERT(b < a);
        TEST_ASSERT(b <= a);

        TEST_ASSERT(a != b);
        TEST_ASSERT(a == a);

        TEST_ASSERT(!(a == b));
        TEST_ASSERT(!(a != a));

        b = -1.;
        a = -0.5;
        TEST_ASSERT(b < a);
        TEST_ASSERT(a > b);
        TEST_ASSERT(b <= a);
        TEST_ASSERT(a >= b);
    }

    void test_atan2() {
        test_atan2_perform<FixedPoint<4,80,uint32_t>, 32>();
        test_atan2_perform<FixedPoint<4,80,uint16_t>, 16>();
        test_atan2_perform<FixedPoint<4,80,uint8_t>, 8>();
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

            std::string msg = std::string("Bits: ") + std::to_string(bits) + ": Iteration " + std::to_string(i) + ", angle: " + std::to_string(atan2(dy,dx)/M_PI*180);
            TEST_ASSERT_DELTA_MSG(double(a), da, 0.5*1e-4 * 2*M_PI, msg.c_str());
        }
    }

    void test_abs() {
        FixedPoint<4, 91, uint32_t> a, b;
        FixedPoint<4, 91, uint16_t> c, d;
        FixedPoint<4, 91, uint8_t> e, f;

        a = c = e = -1.;
        b = d = f = 1.;

        TEST_ASSERT(std::abs(a) == b);
        TEST_ASSERT(std::abs(c) == d);
        TEST_ASSERT(std::abs(e) == f);
        TEST_ASSERT(std::abs(b) == b);
        TEST_ASSERT(std::abs(d) == d);
        TEST_ASSERT(std::abs(f) == f);

        a = c = e = 0.;
        TEST_ASSERT(std::abs(a) == a);
        TEST_ASSERT(std::abs(c) == c);
        TEST_ASSERT(std::abs(e) == e);

        a = FixedPoint<4, 91, uint32_t>::template maxVal<FixedPoint<4, 91, uint32_t>>();
        TEST_ASSERT(std::abs(a) == a);
        TEST_ASSERT(std::abs(-a) == a);

        a = FixedPoint<4, 91, uint32_t>::template minVal<FixedPoint<4, 91, uint32_t>>();
        b = FixedPoint<4, 91, uint32_t>::template maxVal<FixedPoint<4, 91, uint32_t>>();
        TEST_ASSERT(std::abs(a) == b);
        TEST_ASSERT(std::abs(b) == b);
    }
};

int main()
{
    Test::Suite ts;
    ts.add(std::auto_ptr<Test::Suite>(new MultiwordIntegerTestSuite));
    ts.add(std::auto_ptr<Test::Suite>(new FixedPointTestSuite));

    Test::TextOutput output(Test::TextOutput::Verbose);

    return !ts.run(output);
}
