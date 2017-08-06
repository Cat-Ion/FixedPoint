/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file at the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef FIXEDPOINT_HPP
#define FIXEDPOINT_HPP
#define _USE_MATH_DEFINES
#include <assert.h>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <type_traits>

template<typename T> class make_bigger { };
template<> class make_bigger<uint8_t > { public: typedef uint16_t type; };
template<> class make_bigger<uint16_t> { public: typedef uint32_t type; };
template<> class make_bigger<uint32_t> { public: typedef uint64_t type; };

template<typename T> unsigned nlz(T x);
template<> unsigned nlz(uint32_t x) {
    unsigned r = 0;
    if (x <= 0x0000FFFF) { r += 16; x <<= 16; }
    if (x <= 0x00FFFFFF) { r +=  8; x <<=  8; }
    if (x <= 0x0FFFFFFF) { r +=  4; x <<=  4; }
    if (x <= 0x3FFFFFFF) { r +=  2; x <<=  2; }
    if (x <= 0x7FFFFFFF) { r +=  1; }
    return r;
}
template<> unsigned nlz(uint16_t x) {
    unsigned r = 0;
    if (x <= 0x00FF) { r +=  8; x <<=  8; }
    if (x <= 0x0FFF) { r +=  4; x <<=  4; }
    if (x <= 0x3FFF) { r +=  2; x <<=  2; }
    if (x <= 0x7FFF) { r +=  1; }
    return r;
}
template<> unsigned nlz(uint8_t x) {
    unsigned r = 0;
    if (x <= 0x0F) { r +=  4; x <<=  4; }
    if (x <= 0x3F) { r +=  2; x <<=  2; }
    if (x <= 0x7F) { r +=  1; }
    return r;
}

template<unsigned size, typename _storageType = uint32_t>
class MultiwordInteger
{
public:
    typedef _storageType storageType;
    static const constexpr size_t numWords = size;
    static const constexpr size_t storageSize = sizeof(storageType) * 8;

    template<unsigned otherSize, typename otherStorageType>
    friend class MultiwordInteger;

protected:
    typedef typename std::make_signed<_storageType>::type signedType;
    typedef typename make_bigger<_storageType>::type bigType;
    _storageType s[size];

public:
    constexpr MultiwordInteger() : s{0} {}
    template<unsigned otherSize> constexpr MultiwordInteger(MultiwordInteger<otherSize, storageType> const &o);
    constexpr MultiwordInteger(storageType const &v);
    constexpr MultiwordInteger(int64_t v);
    constexpr MultiwordInteger(double v);
    template<unsigned otherSize, typename otherStorageType> constexpr MultiwordInteger(MultiwordInteger<otherSize, otherStorageType> const &o);

    static constexpr
    MultiwordInteger<size, storageType>
    maxVal() {
        MultiwordInteger<size, storageType> r;
        for (unsigned i = 0; i < size-1; i++) {
            r.s[i] = (bigType(1)<<(storageSize)) - 1;;
        }
        r.s[size-1] = (1U<<(storageSize - 1)) - 1;
        return r;
    }

    static constexpr
    MultiwordInteger<size, storageType>
    minVal() {
        MultiwordInteger<size, storageType> r;
        for (unsigned i = 0; i < size-1; i++) {
            r.s[i] = 0;
        }
        r.s[size-1] = 1U<<(storageSize-1);
        return r;
    }

    constexpr MultiwordInteger<size, storageType>& operator+=(MultiwordInteger<size, storageType> const &o);
    constexpr MultiwordInteger<size, storageType>& operator-=(MultiwordInteger<size, storageType> const &o);
    template<unsigned otherSize> constexpr MultiwordInteger<size, storageType>& operator*=(MultiwordInteger<otherSize, storageType> const &o);
    template<unsigned otherSize> constexpr MultiwordInteger<size, storageType>& operator/=(MultiwordInteger<otherSize, storageType> const &o);
    constexpr MultiwordInteger<size, storageType>& operator%=(MultiwordInteger<size, storageType> const &o);

    constexpr MultiwordInteger<size, storageType>& operator++();
    constexpr MultiwordInteger<size, storageType> operator++(int);
    constexpr MultiwordInteger<size, storageType>& operator--();
    constexpr MultiwordInteger<size, storageType>& operator--(int);

    constexpr MultiwordInteger<size, storageType> operator-() const;

    constexpr MultiwordInteger<size, storageType>& operator &=(MultiwordInteger<size, storageType> const &o);
    constexpr MultiwordInteger<size, storageType>& operator |=(MultiwordInteger<size, storageType> const &o);
    constexpr MultiwordInteger<size, storageType>& operator ^=(MultiwordInteger<size, storageType> const &o);
    constexpr MultiwordInteger<size, storageType>  operator ~() const;
    constexpr MultiwordInteger<size, storageType>& operator <<=(unsigned n);
    constexpr MultiwordInteger<size, storageType>& operator >>=(unsigned n);

    constexpr bool operator<(MultiwordInteger<size, storageType> const &o) const;
    constexpr bool operator==(MultiwordInteger<size, storageType> const &o) const;

    constexpr void negate();
    constexpr void fill_leading_bits(unsigned num);

    constexpr unsigned leading_zeros() const;
    constexpr bool     is_negative() const;
    constexpr bool     is_positive() const;
    // Return true if bit at position is set, starting from the LSB.
    constexpr bool bit(size_t position) const;

    explicit constexpr operator bool() const;
    explicit constexpr operator int8_t() const;
    explicit constexpr operator int16_t() const;
    explicit constexpr operator int32_t() const;
    explicit constexpr operator int64_t() const;
    explicit constexpr operator double() const;

    template<unsigned otherSize, unsigned outSize> constexpr void
    mul(MultiwordInteger<otherSize, storageType> const &o,
        MultiwordInteger<outSize,   storageType>       *out) const;

    template<unsigned otherSize> constexpr void
    quotrem(MultiwordInteger<otherSize, storageType> divisor,
            MultiwordInteger<size,      storageType> &qotient,
            MultiwordInteger<otherSize, storageType> *remainder) const;

protected:
    /* Division by a divisor of length 1
     * The remainder is put into the dividend.
     */
    template<unsigned divisorSize> constexpr void
    unsignedQuotrem(MultiwordInteger<size+1,      storageType> &dividend,
                    MultiwordInteger<divisorSize, storageType> const &divisor,
                    MultiwordInteger<size,        storageType> &quotient,
                    unsigned dividendSize) const;

    /* Division by a longer divisor
     * The remainder is put into the dividend.
     */
    template<unsigned divisorSize> constexpr void
    unsignedQuotrem(MultiwordInteger<size+1,      storageType> &dividend,
                    MultiwordInteger<divisorSize, storageType>  divisor,
                    MultiwordInteger<size,        storageType> &quotient,
                    unsigned dividend_length, unsigned divisor_length,
                    unsigned divisor_nlz) const;
};

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType = uint32_t>
class FixedPoint
{
    static_assert((_integerWidth + _fractionalWidth + 1) % (sizeof(backingStorageType) * 8) == 0, "Integer and fractional width do not match size of storage type.");
public:
    typedef MultiwordInteger<(_integerWidth+_fractionalWidth+1)/(sizeof(backingStorageType)*8), backingStorageType> StorageType;
protected:
    typedef FixedPoint<_integerWidth,_fractionalWidth,backingStorageType> FP;
public:
    static const unsigned fractionalWidth = _fractionalWidth;
    static const int integerWidth = _integerWidth;
    StorageType v;

    constexpr
    FixedPoint()
    {
    }

    constexpr
    FixedPoint(double v)
    {
        if (v > maxVal<double>()) {
            v = maxVal<double>();
        } else if (v < minVal<double>()) {
            v = minVal<double>();
        }
        this->v = v * pow(2., (double)_fractionalWidth);
    }

    constexpr
    FixedPoint(int v)
    {
        if (v > 1 << _integerWidth) {
            this->v = StorageType::maxVal();
        } else if(v < -(1<<_integerWidth)) {
            this->v = StorageType::minVal();
        } else {
            this->v = backingStorageType(v);
            this->v <<= _fractionalWidth;
        }
    }

    constexpr
    FixedPoint(FixedPoint const &o) : v(o.v) {}

    constexpr
    FixedPoint(StorageType const &s) : v(s) {}

    template<int oiw, unsigned ofw, typename otherStorageType>
    constexpr
    FixedPoint(FixedPoint<oiw,ofw,otherStorageType> const &o) {
        if (ofw > _fractionalWidth) {
            typename FixedPoint<oiw, ofw, otherStorageType>::StorageType ov = o.v >> (ofw-_fractionalWidth);
            this->v = ov;
        } else {
            this->v = o.v;
            this->v <<= (_fractionalWidth - ofw);
        }
    }

    template<typename T>
    constexpr static
    T
    maxVal()
    {
        return T(FP(StorageType::maxVal()));
    }

    template<typename T>
    constexpr static
    T
    minVal()
    {
        return T(FP(StorageType::minVal()));
    }

    template<typename T>
    constexpr static
    T
    smallestVal()
    {
        return T(FP(StorageType(backingStorageType(1))));
    }

    explicit constexpr
    operator double() const
    {
        double dv = double(v);
        double p2 = pow(2., -(double) _fractionalWidth);
        double r = p2 * dv;
        return r;
    }

    constexpr static
    FP
    from_int(int v)
    {
        FP r;
        r.v = int64_t(v);
        return r;
    }

    constexpr
    FP
    operator+=(FP const &o) {
        StorageType z(backingStorageType(0));
        if (o.v.is_negative() && v.is_negative()) {
            StorageType d(StorageType::minVal());
            d -= v;
            assert (d <= o.v);
        } else if(o.v.is_positive() && v.is_positive()) {
            StorageType d(StorageType::maxVal());
            d -= v;
            assert (d >= o.v);
        }
        this->v += o.v;
        return *this;
    }

    constexpr
    FP
    operator-=(FP const &o) {
        *this += -o;
        return *this;
    }

    constexpr
    FP
    operator*=(FP const &o) {
        assert((double)*this * (double)o <= maxVal<double>());
        assert((double)*this * (double)o >= minVal<double>());
        MultiwordInteger<StorageType::numWords*2, backingStorageType> nv;
        nv = this->v * o.v;
        this->v = nv >> _fractionalWidth;
        if (_fractionalWidth > 0 && nv.bit(_fractionalWidth - 1)) {
            ++this->v;
        }
        return *this;
    }

    constexpr
    FP
    operator*=(int64_t const &o) {
        assert((double)*this * (double)o <= maxVal<double>());
        assert((double)*this * (double)o >= minVal<double>());
        this->v *= o;
        return *this;
    }

    constexpr
    FP&
    operator/=(FP const &o) {
        assert((double)*this / (double)o <= maxVal<double>());
        assert((double)*this / (double)o >= minVal<double>());
        MultiwordInteger<StorageType::numWords+(_fractionalWidth+StorageType::storageSize-1)/StorageType::storageSize, backingStorageType> nv(this->v);
        nv <<= _fractionalWidth;
        nv /= o.v;
        v = nv;
        return *this;
    }

    constexpr
    FP&
    operator/=(int64_t const &o) {
        assert((double)*this / (double)o <= maxVal<double>());
        assert((double)*this / (double)o >= minVal<double>());
        v /= o;
        return *this;
    }

    constexpr
    FP&
    operator%=(FP const &o) {
        v %= o.v;
        return *this;
    }

    constexpr
    FP&
    operator=(FP const &o)
    {
      v = o.v;
      return *this;
    }

    friend constexpr
    FP
    operator-(FP right) {
        right.v = -right.v;
        return right;
    }

    friend constexpr
    bool
    operator==(FP const &left, FP const &right) {
        return left.v == right.v;
    }

    friend constexpr
    bool
    operator<(FP const &left, FP const &right) {
        return left.v < right.v;
    }

    constexpr
    bool
    is_negative() const {
      return v.is_negative();
    }

    constexpr
    bool
    is_positive() const {
      return v.is_positive();
    }

    constexpr
    FP
    nabs() const {
        return v.is_positive() ? -*this : *this;
    }

    explicit constexpr operator int8_t()  const { return int8_t (v>>_fractionalWidth); }
    explicit constexpr operator int16_t() const { return int16_t(v>>_fractionalWidth); }
    explicit constexpr operator int32_t() const { return int32_t(v>>_fractionalWidth); }
    explicit constexpr operator int64_t() const { return int64_t(v>>_fractionalWidth); }
};

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator>>(MultiwordInteger<size, storageType> left, unsigned right) { return left >>= right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator<<(MultiwordInteger<size, storageType> left, unsigned right) { return left <<= right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator&(MultiwordInteger<size, storageType> left, MultiwordInteger<size, storageType> const &right) { return left &= right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator|(MultiwordInteger<size, storageType> left, MultiwordInteger<size, storageType> const &right) { return left |= right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator^(MultiwordInteger<size, storageType> left, MultiwordInteger<size, storageType> const &right) { return left ^= right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator+ (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return left += right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator- (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return left -= right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator* (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return left *= right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator/ (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return left /= right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator% (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return left %= right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr bool operator> (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &left,
                          FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return right < left; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr bool operator<=(FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &left,
                          FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return !(left > right); }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr bool operator>=(FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &left,
                          FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return !(left < right); }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr bool operator!=(FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &left,
                          FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return !(left == right); }


namespace std {
    template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, mfw, backingStorageType>
    abs(FixedPoint<iw,mfw,backingStorageType> const &y) {
        if (y.is_negative()) {
            if (y == FixedPoint<iw, mfw,backingStorageType>::template minVal<FixedPoint<iw, mfw,backingStorageType>>()) {
                return FixedPoint<iw, mfw, backingStorageType>::template maxVal<FixedPoint<iw, mfw, backingStorageType>>();
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
        FixedPoint<iw,mfw,backingStorageType>
                mx = mxmn.first,
                mn = mxmn.second,
                a = mn / mx,
                s = a*a,
                r(-0.0464964749);
        FixedPoint<iw, mfw, backingStorageType> const
                c1(0.15931422),
                c2(-0.327622764),
                pi2(M_PI_2),
                pi(M_PI);

        r = r * s + c1;
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

        FixedPoint<iw, mfw, backingStorageType> c1(.405284735), c2(1.27323954), c3(0.225);
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

        FixedPoint<iw, mfw, backingStorageType> c1(.405284735), c2(1.27323954), c3(0.225);
        FixedPoint<iw, mfw, backingStorageType> ret;
        ret = c1 * x.nabs();
        ret += c2;
        ret *= x;

        ret += c3 * (ret * std::abs(ret) - ret);

        return ret;
    }
}

#include "MultiwordIntegerConstructors.hpp"
#include "MultiwordIntegerArithmetics.hpp"
#include "MultiwordIntegerBitwise.hpp"
#include "MultiwordIntegerCasts.hpp"
#include "MultiwordIntegerComparisons.hpp"
#include "MultiwordIntegerUtility.hpp"
#endif // FIXEDPOINT_HPP
