FixedPoint
==========

To use:

    #include <FixedPoint.hpp>

Then declare some variables with an integer width, a (minimum) fractional
width, and optionally an unsigned storage type, which defaults to 32-bit ints.

The variables

    FixedPoint<5, 32> a;
    FixedPoint<5, 35, uint32_t> b;
    FixedPoint<5, 58, uint32_t> c;

have the same integer precision of 5 bits, and the same fractional precision
of

    32 * ceil(5+35)/32 - 5 - 1 = 58

bits - the minimal storage size is rounded up to a multiple of the size of the
storage type; then the bits needed for the integral part and sign bit are subtracted.

Arithmetics work as expected:

    a = 5.;
    b = 6.;

    c = a/b;

    printf("%g / %g = %g\n", double(a), double(b), double(c));

will print

    5 / 6 = 0.833333

For more examples, look at the provided tests in FixedPointTests.cpp

