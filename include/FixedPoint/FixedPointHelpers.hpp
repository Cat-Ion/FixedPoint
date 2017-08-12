#ifndef FIXEDPOINTHELPERS_HPP
#define FIXEDPOINTHELPERS_HPP
#include <cmath>
#include <limits.h>

namespace FixedPointHelpers {
    constexpr double dipow(double base, int exponent) {
        double ret = 1.;
        if (exponent < 0) {
            base = 1./base;
            exponent = -exponent;
        }
        while (exponent > 0) {
            if (exponent & 1) {
                ret *= base;
            }
            exponent /= 2;
            base *= base;
        }
        return ret;
    }

    constexpr int ilogb(double v) {
        return v < 0         ? ilogb(-v) :
               std::isnan(v) ? 0 :
               std::isinf(v) ? INT_MAX :
               v == 0        ? INT_MIN :
               v < 1         ? ilogb(v*2)-1 :
               v >= 2        ? ilogb(v/2) + 1 :
                               0;
    }
}
#endif // FIXEDPOINTHELPERS_HPP
