#ifndef FIXEDPOINTHELPERS_HPP
#define FIXEDPOINTHELPERS_HPP
#include <cmath>
#include <limits.h>

namespace FixedPointHelpers {
    constexpr double dipow_accum(double base, int exponent, double accum) {
        if (exponent == 0) {
            return accum;
        }
        if (exponent < 0) {
            return dipow_accum(base, exponent + 1, accum / base);
        }

        return dipow_accum(base, exponent - 1, accum * base);
    }

    constexpr double dipow(double base, int exponent) {
        return dipow_accum(base, exponent, 1.);
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
