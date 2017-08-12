#ifndef FIXEDPOINTSTD_HPP
#define FIXEDPOINTSTD_HPP
#include "FixedPoint.hpp"
#define _USE_MATH_DEFINES
#include <cmath>

namespace std {
    template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, mfw, backingStorageType>
    abs(FixedPoint<iw,mfw,backingStorageType> const &y) {
        if (y.is_negative()) {
            if (y == FixedPoint<iw, mfw,backingStorageType>::_minVal()) {
                return FixedPoint<iw, mfw, backingStorageType>::_maxVal();
            } else {
                return -y;
            }
        } else {
            return y;
        }
    }

    template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw,mfw, backingStorageType>
    atan2(FixedPoint<iw,mfw,backingStorageType> const &y, FixedPoint<iw,mfw,backingStorageType> const &x)
    {
        FixedPoint<iw,mfw,backingStorageType> ya = y.nabs(), xa = x.nabs();
        auto mxmn = std::minmax(ya, xa);
        constexpr FixedPoint<iw, mfw, backingStorageType>
                c1(0.15931422),
                c2(-0.327622764),
                pi2(M_PI_2),
                pi(M_PI),
                rstart(-0.0464964749);
        FixedPoint<iw,mfw,backingStorageType>
                mx = mxmn.first,
                mn = mxmn.second,
                a = mn / mx,
                s = a*a;

        FixedPoint<iw,mfw,backingStorageType> r = rstart * s + c1;
        r = r * s + c2;
        r = r * s * a + a;

        if (ya < xa) {
            r = pi2 - r;
        }
        if (x.is_negative()) {
            r = pi - r;
        }
        if (y.is_negative()) {
            r = -r;
        }

        return r;
    }

    template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, mfw, backingStorageType>
    sin(FixedPoint<iw, mfw, backingStorageType> x) {
        if (x.is_negative()) {
            x %= FixedPoint<iw, mfw, backingStorageType>(-2*M_PI);
        } else {
            x %= FixedPoint<iw, mfw, backingStorageType>(2*M_PI);
        }

        if (x > FixedPoint<iw, mfw, backingStorageType>(M_PI)) {
            x -= FixedPoint<iw, mfw, backingStorageType>(2*M_PI);
        } else if (x < FixedPoint<iw, mfw, backingStorageType>(-M_PI)) {
            x += FixedPoint<iw, mfw, backingStorageType>(2*M_PI);
        }

        constexpr FixedPoint<iw, mfw, backingStorageType> c1(.405284735), c2(1.27323954), c3(0.225);
        FixedPoint<iw, mfw, backingStorageType> ret;
        ret = c1 * x.nabs();
        ret += c2;
        ret *= x;

        ret += c3 * (ret * std::abs(ret) - ret);

        return ret;
    }

    template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, mfw, backingStorageType>
    cos(FixedPoint<iw, mfw, backingStorageType> x) {
        if (x.is_negative()) {
            x %= FixedPoint<iw, mfw, backingStorageType>(-2*M_PI);
        } else {
            x %= FixedPoint<iw, mfw, backingStorageType>(2*M_PI);
        }

        x += FixedPoint<iw, mfw, backingStorageType>(M_PI_2);
        if (x > FixedPoint<iw, mfw, backingStorageType>(M_PI)) {
            x -= FixedPoint<iw, mfw, backingStorageType>(2*M_PI);
        } else if (x < FixedPoint<iw, mfw, backingStorageType>(-M_PI)) {
            x += FixedPoint<iw, mfw, backingStorageType>(2*M_PI);
        }

        constexpr FixedPoint<iw, mfw, backingStorageType> c1(.405284735), c2(1.27323954), c3(0.225);
        FixedPoint<iw, mfw, backingStorageType> ret;
        ret = c1 * x.nabs();
        ret += c2;
        ret *= x;

        ret += c3 * (ret * std::abs(ret) - ret);

        return ret;
    }
}
#endif // FIXEDPOINTSTD_HPP
