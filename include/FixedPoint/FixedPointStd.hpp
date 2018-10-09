/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
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

    template<int iw, int fw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, fw, backingStorageType>
    exp(FixedPoint<iw, fw, backingStorageType> x) {
        using Type = FixedPoint<iw, fw, backingStorageType>;
        Type val = x;
        Type series = Type(1) + val;
        for(int i = 2; val != Type(0); i++) {
            val = val * x / i;
            series += val;
        }
        return series;
    }

    template<int iw, int fw, typename backingStorageType = uint32_t>
    FixedPoint<iw, fw, backingStorageType>
    ln2() {
        using Type = FixedPoint<iw, fw + 32, backingStorageType>;
        Type retval = 0.;
        Type div = 1;
        Type el = 1;
        for (int n = 1; el != Type(0); n++) {
            div.v >>= 1;
            el = div / n;
            retval += el;
        }
        return retval;
    }
    
    template<int iw, int fw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, fw, backingStorageType>
    log(FixedPoint<iw, fw, backingStorageType> x) {
        using Type = FixedPoint<iw, fw, backingStorageType>;
        Type retval = 0;
        if (x < Type(0)) {
            return FixedPoint<iw, fw, backingStorageType>::minVal;
        }

        while (x > Type(1)) {
            retval += ln2<iw, fw, backingStorageType>();
            x.v >>= 1;
        }
        
        // ln(y) = ln(1-x) = \sum_{i=1}^{\infty} -x^i/i
        // y = 1-x, x = 1-y
        x = Type(1) - x;
        Type el = 1;
        Type eld = 1;
        for (int i = 1; eld != Type(0); i++) {
            el *= x;
            eld = el / i;
            retval -= eld;
        }
        return retval;
    }

    template<int iw, int fw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, fw, backingStorageType>
    pow(FixedPoint<iw, fw, backingStorageType> a, int b) {
        if (b < 0) {
            return FixedPoint<iw, fw, backingStorageType>(1) / pow(a, -b);
        }
        
        FixedPoint<iw, fw, backingStorageType> retval = 1;
        while(b) {
            if (b & 1) {
                retval *= a;
            }
            b >>= 1;
            a *= a;
        }
        return retval;
    }
    
    template<int iw, int fw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, fw, backingStorageType>
    pow(FixedPoint<iw, fw, backingStorageType> const &a, FixedPoint<iw, fw, backingStorageType> b) {
        using Type = FixedPoint<iw, fw, backingStorageType>;
        Type bloga = log(a) * b;
        Type val = bloga;
        Type series = val + Type(1);
        for(int i = 2; val != Type(0); i++) {
            val = val * bloga / i;
            series += val;
        }
        return series * pow(a, int(b));
    }
}
#endif // FIXEDPOINTSTD_HPP
