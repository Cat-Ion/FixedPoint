#ifndef FIXEDPOINTSTD_HPP
#define FIXEDPOINTSTD_HPP
#include "FixedPoint.hpp"
#define _USE_MATH_DEFINES
#include <cmath>

namespace std {
    template<int iw, int fw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, fw, backingStorageType>
    abs(FixedPoint<iw,fw,backingStorageType> const &y) {
        if (y.is_negative()) {
            if (y == FixedPoint<iw, fw,backingStorageType>::_minVal()) {
                return FixedPoint<iw, fw, backingStorageType>::_maxVal();
            } else {
                return -y;
            }
        } else {
            return y;
        }
    }

    template<int iw, int fw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw,fw, backingStorageType>
    atan2(FixedPoint<iw,fw,backingStorageType> const &y, FixedPoint<iw,fw,backingStorageType> const &x)
    {
        FixedPoint<iw,fw,backingStorageType> ya = y.nabs(), xa = x.nabs();
        auto mxmn = std::minmax(ya, xa);
        constexpr FixedPoint<iw, fw, backingStorageType>
                c1(0.15931422),
                c2(-0.327622764),
                pi2(M_PI_2),
                pi(M_PI),
                rstart(-0.0464964749);
        FixedPoint<iw,fw,backingStorageType>
                mx = mxmn.first,
                mn = mxmn.second,
                a = mn / mx,
                s = a*a;

        FixedPoint<iw,fw,backingStorageType> r = rstart * s + c1;
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

    template<int iw, int fw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, fw, backingStorageType>
    sin(FixedPoint<iw, fw, backingStorageType> x) {
        if (x.is_negative()) {
            x %= FixedPoint<iw, fw, backingStorageType>(-2*M_PI);
        } else {
            x %= FixedPoint<iw, fw, backingStorageType>(2*M_PI);
        }

        if (x > FixedPoint<iw, fw, backingStorageType>(M_PI)) {
            x -= FixedPoint<iw, fw, backingStorageType>(2*M_PI);
        } else if (x < FixedPoint<iw, fw, backingStorageType>(-M_PI)) {
            x += FixedPoint<iw, fw, backingStorageType>(2*M_PI);
        }

        constexpr FixedPoint<iw, fw, backingStorageType> c1(.405284735), c2(1.27323954), c3(0.225);
        FixedPoint<iw, fw, backingStorageType> ret;
        ret = c1 * x.nabs();
        ret += c2;
        ret *= x;

        ret += c3 * (ret * std::abs(ret) - ret);

        return ret;
    }

    template<int iw, int fw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, fw, backingStorageType>
    cos(FixedPoint<iw, fw, backingStorageType> x) {
        if (x.is_negative()) {
            x %= FixedPoint<iw, fw, backingStorageType>(-2*M_PI);
        } else {
            x %= FixedPoint<iw, fw, backingStorageType>(2*M_PI);
        }

        x += FixedPoint<iw, fw, backingStorageType>(M_PI_2);
        if (x > FixedPoint<iw, fw, backingStorageType>(M_PI)) {
            x -= FixedPoint<iw, fw, backingStorageType>(2*M_PI);
        } else if (x < FixedPoint<iw, fw, backingStorageType>(-M_PI)) {
            x += FixedPoint<iw, fw, backingStorageType>(2*M_PI);
        }

        constexpr FixedPoint<iw, fw, backingStorageType> c1(.405284735), c2(1.27323954), c3(0.225);
        FixedPoint<iw, fw, backingStorageType> ret;
        ret = c1 * x.nabs();
        ret += c2;
        ret *= x;

        ret += c3 * (ret * std::abs(ret) - ret);

        return ret;
    }
}
#endif // FIXEDPOINTSTD_HPP
