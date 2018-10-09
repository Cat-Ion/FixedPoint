FixedPoint
==========

To use, get a compiler that supports C++14 (gcc and clang are supported, Visual
Studio will currently not compile the code, and other compilers have not been
tested), and

    #include <FixedPoint/FixedPoint.hpp>

Then declare some variables with an integer width, a fractional width, and
optionally an unsigned storage type, which defaults to 32-bit ints. The integer
and fractional width need to add up to one less than a multiple of the size of
the storage type.

The variables

    FixedPoint<5, 26> a;
    FixedPoint<5, 26, uint32_t> b;
    FixedPoint<5, 26, uint16_t> c;

have the same integer precision of 5 bits, and the same fractional precision
of 26 bits. `a` and `b` both use 32-bit numbers for storage, while `c` uses 16-bit
numbers.

Arithmetics work as expected:

    FixedPoint<5, 26> a, b, c;
    a = 5.;
    b = 6.;

    c = a/b;

    printf("%g / %g = %g\n", double(a), double(b), double(c));

will print

    5 / 6 = 0.833333

Other than the standard arithmetic (`+`, `-`, `*`, `/`, `%` and their
assignment equivalents) and comparison operators, the following functions
are defined as well:

    std::sin
    std::cos
    std::atan2
    std::abs
    std::exp
    std::pow
    std::log

For the trigonometric functions, the results are at the moment only tested
with an integral width of 3 and will probably return nonsensical results
with smaller widths.

For more examples, look at the provided tests in FixedPointTests.cpp

