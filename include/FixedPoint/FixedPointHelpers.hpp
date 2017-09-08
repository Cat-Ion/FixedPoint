/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef FIXEDPOINTHELPERS_HPP
#define FIXEDPOINTHELPERS_HPP
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdint.h>

namespace FixedPointHelpers {
    template<typename T> class make_bigger { };
    template<> class make_bigger<uint8_t > { public: typedef uint16_t type; };
    template<> class make_bigger<uint16_t> { public: typedef uint32_t type; };
    template<> class make_bigger<uint32_t> { public: typedef uint64_t type; };

    template<typename T>
    constexpr unsigned nlz_constexpr(T x_) {
        typedef typename std::make_unsigned<T>::type UT;
        UT x = x_;
        unsigned rv = 0;
        for (UT i = UT(1)<<(sizeof(UT)*8-1); i; i >>= 1) {
            if (x & i) {
                return rv;
            }
            rv++;
        }
        return rv;
    }

    template<typename T> inline constexpr unsigned nlz(T x) {
        return nlz_constexpr(x);
    }

    template<> inline constexpr unsigned nlz(unsigned long long x) {
        return __builtin_constant_p(x) ? nlz_constexpr(x) : __builtin_clzll(x);
    }
    template<> inline constexpr unsigned nlz(unsigned long x) {
        return __builtin_constant_p(x) ? nlz_constexpr(x) : __builtin_clzl(x);
    }
    template<> inline constexpr unsigned nlz(unsigned int x) {
        return __builtin_constant_p(x) ? nlz_constexpr(x) : __builtin_clz(x);
    }

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
               !(v == v)     ? 0 : // NAN check
               v > std::numeric_limits<double>::max() ? std::numeric_limits<int>::max() :
               v == 0        ? std::numeric_limits<int>::min() :
               v < 1         ? ilogb(v*2)-1 :
               v >= 2        ? ilogb(v/2) + 1 :
                               0;
    }
}
#endif // FIXEDPOINTHELPERS_HPP
