#ifndef FIXEDPOINTHELPERS_HPP
#define FIXEDPOINTHELPERS_HPP
namespace FixedPointHelpers {
    constexpr double dipow(double base, int exponent) {
        return exponent == 0 ? 1. :
               exponent  < 0 ? dipow(base, exponent + 1) / base :
                               dipow(base, exponent - 1) * base;
    }
}
#endif // FIXEDPOINTHELPERS_HPP
