/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file at the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef FIXEDPOINT_HPP
#define FIXEDPOINT_HPP
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
    MultiwordInteger() {}
    template<unsigned otherSize>
    constexpr
    MultiwordInteger(MultiwordInteger<otherSize, storageType> const &o) {
        unsigned num = size > otherSize ? otherSize : size;
        for (unsigned i = 0; i < num; i++) {
            s[i] = o.s[i];
        }
        if (static_cast<signedType>(o.s[otherSize-1]) < 0) {
            for (unsigned i = num; i < size; i++) {
                s[i] = ~static_cast<storageType>(0);
            }
        } else {
            for (unsigned i = num; i < size; i++) {
                s[i] = 0;
            }
        }
    }

    constexpr
    MultiwordInteger(storageType const &v) {
        s[0] = v;
        if (static_cast<signedType>(v) < 0) {
            for (unsigned i = 1; i < size; i++) {
                s[i] = static_cast<signedType>(-1);
            }
        } else {
            for (unsigned i = 1; i < size; i++) {
                s[i] = 0;
            }
        }
    }

    constexpr
    MultiwordInteger(int64_t v) {
        unsigned i = 0;
        typename std::make_unsigned<int64_t>::type uv = v;
        while (uv && i < size) {
            s[i] = uv & ((1UL<<storageSize) - 1);
            uv >>= storageSize;
            i++;
        }
        while (i < size) {
            s[i] = 0;
            i++;
        }
        if (v < 0) {
            fill_leading_bits(leading_zeros());
        }
    }

    constexpr
    MultiwordInteger(double v) {
        double sig = significand(v);
        int lg = ilogb(v);
        if (lg < 63) {
            int64_t i = sig*(1UL<<lg);
            *this = i;
        } else {
            *this = int64_t(sig * (1UL<<62));
            *this <<= lg - 62;
        }
    }

    template<unsigned otherSize, typename otherStorageType>
    constexpr
    MultiwordInteger(MultiwordInteger<otherSize, otherStorageType> const &o) {
        static_assert((sizeof(storageType) % sizeof(otherStorageType)) == 0
                      || (sizeof(otherStorageType) % sizeof(storageType)) == 0,
                      "Types must fit into each other without remainder.");
        if (sizeof(otherStorageType) < sizeof(storageType)) {
            unsigned shiftNum = sizeof (storageType) / sizeof (otherStorageType);
            unsigned shiftWidth = sizeof (otherStorageType) * 8;
            unsigned i = 0;
            for (i = 0; i < size && i*shiftNum < otherSize; i++) {
                s[i] = 0;
                unsigned start = shiftNum;
                if (i*shiftNum + start > otherSize) {
                    start = otherSize - i*shiftNum;
                    if (o.is_negative()) {
                        s[i] |= (~static_cast<storageType>(0)) << (start * shiftWidth);
                    }
                }
                for (unsigned j = start; j-- > 0; ) {
                    s[i] |= o.s[i*shiftNum + j] << (shiftWidth * j);
                }
            }
            if (o.is_negative()) {
                while (i < size) {
                    s[i] = ~static_cast<storageType>(0);
                    i++;
                }
            } else {
                while (i < size) {
                    s[i] = 0;
                    i++;
                }
            }
        } else {
            unsigned shiftNum = sizeof (otherStorageType) / sizeof (storageType);
            unsigned shiftWidth = sizeof (storageType) * 8;
            unsigned i = 0, j = 0;
            for (i = j = 0; i < size && j < otherSize; j++) {
                otherStorageType v = o.s[j];
                for (unsigned k = 0; k < shiftNum && i < size; k++, i++) {
                    s[i] = v;
                    v >>= shiftWidth;
                }
            }
            if (o.is_negative()) {
                while (i < size) {
                    s[i] = ~static_cast<storageType>(0);
                    i++;
                }
            } else {
                while (i < size) {
                    s[i] = 0;
                    i++;
                }
            }
        }
    }

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

    constexpr
    MultiwordInteger<size, storageType>&
    operator+=(MultiwordInteger<size, storageType> const &o) {
        storageType c = 0;
        for(unsigned i = 0; i < size; i++) {
            bigType t = this->s[i] + o.s[i] + c;
            this->s[i] = t;
            c = t >> storageSize;
        }
        return *this;
    }

    constexpr
    MultiwordInteger<size, storageType>&
    operator-=(MultiwordInteger<size, storageType> const &o) {
        bigType c = 0;
        for (unsigned i = 0; i < size; i++) {
            c += this->s[i];
            c -= o.s[i];
            this->s[i] = c;
            c = (signedType)(c>>storageSize);
        }
        return *this;
    }

    template<unsigned otherSize>
    constexpr
    MultiwordInteger<size, storageType>&
    operator*=(MultiwordInteger<otherSize, storageType> const &o) {
        MultiwordInteger<size, storageType> nv;
        mul<otherSize, size>(o, &nv);
        *this = nv;
        return *this;
    }

    template<unsigned otherSize>
    constexpr
    MultiwordInteger<size, storageType>&
    operator/=(MultiwordInteger<otherSize, storageType> const &o) {
        quotrem(o, *this, static_cast<MultiwordInteger<otherSize, storageType>*>(nullptr));
        return *this;
    }

    constexpr
    MultiwordInteger<size, storageType>&
    operator++() {
        bigType c = 1;
        for (unsigned i = 0; i < size && c; i++) {
            c += s[i];
            s[i] = c;
            c >>= storageSize;
        }
        return *this;
    }

    constexpr
    MultiwordInteger<size, storageType>
    operator++(int) {
        MultiwordInteger<size, storageType> r(*this);
        ++(*this);
        return r;
    }

    constexpr
    MultiwordInteger<size, storageType>&
    operator--() {
        bigType c = ~static_cast<bigType>(0);
        for (unsigned i = 0; i < size && c; i++) {
            c += s[i];
            s[i] = c;
            c >>= storageSize;
        }
        return *this;
    }

    constexpr
    MultiwordInteger<size, storageType>&
    operator--(int) {
        MultiwordInteger<size, storageType> r(*this);
        *this--;
        return r;
    }

    constexpr
    MultiwordInteger<size, storageType>&
    operator%=(MultiwordInteger<size, storageType> const &o) {
        *this = *this % o;
    }

    constexpr
    MultiwordInteger<size, storageType>&
    operator<<=(unsigned n) {
        size_t w = n / storageSize;
        n %= storageSize;

        if (w >= size) {
            *this = int64_t(0);
            return *this;
        }

        if (n == 0) {
            for (unsigned i = size; i-- > w; ) {
                s[i] = s[i-w];
            }
            for (unsigned i = 0; i < w; i++) {
                s[i] = 0;
            }
            return *this;
        }

        unsigned i = size;
        while (i-- > (1+w)) {
            this->s[i] = (this->s[i - w] << n)
                    + (this->s[i - w - 1] >> (storageSize - n));
        }
        this->s[i] = (this->s[i - w] << n) ;
        while (i-- > 0) {
            this->s[i] = 0;
        }

        return *this;
    }

    constexpr
    MultiwordInteger<size, storageType>&
    operator>>=(unsigned n) {
        size_t w = n / storageSize;
        n %= storageSize;

        bool adjust_leading_bits = this->is_negative();

        if (w >= size) {
            if (adjust_leading_bits) {
                for (unsigned i = 0; i < size; i++) {
                    s[i] = ~static_cast<storageType>(0);
                }
            } else {
                *this = storageType(0);
            }
            return *this;
        }

        if (n == 0) {
            for (unsigned i = 0; i < size-w; i++) {
                s[i] = s[i+w];
            }
            for (unsigned i = size-w; i < size; i++) {
                s[i] = 0;
            }
            if (adjust_leading_bits) {
                this->fill_leading_bits(w*storageSize);
            }
            return *this;
        }

        unsigned i = 0;
        while (i < (size - w - 1)) {
            this->s[i] = (this->s[i + w] >> n) + (this->s[i + w + 1] << (storageSize - n));
            i++;
        }
        this->s[i] = this->s[i + w] >> n;
        while (++i < size) {
            this->s[i] = 0;
        }

        if (adjust_leading_bits) {
            this->fill_leading_bits(w * storageSize + n);
        }
        return *this;
    }

    constexpr
    MultiwordInteger<size, storageType>
    operator-() const {
        MultiwordInteger<size, storageType> r(*this);
        r.negate();
        return r;
    }

    constexpr
    bool
    operator<(MultiwordInteger<size, storageType> const &o) const {
        if (static_cast<signedType>(s[size-1]) < static_cast<signedType>(o.s[size-1])) {
            return true;
        } else if(static_cast<signedType>(s[size-1]) > static_cast<signedType>(o.s[size-1])) {
            return false;
        }
        for (unsigned i = size - 1; i-- > 0; ) {
            if (s[i] < o.s[i]) {
                return true;
            } else if (s[i] > o.s[i]) {
                return false;
            }
        }
        return false;
    }

    constexpr
    bool
    operator==(MultiwordInteger<size, storageType> const &o) const {
        for (unsigned i = size; i-- > 0; ) {
            if (s[i] != o.s[i]) {
                return false;
            }
        }
        return true;
    }

    template<unsigned otherSize, unsigned outSize>
    constexpr
    void
    mul(MultiwordInteger<otherSize, storageType> const &o,
        MultiwordInteger<outSize,   storageType>       *out) const {
        *out = (int64_t) 0;

        unsigned limitThis = size > outSize ? outSize : size;
        for (unsigned i = 0; i < limitThis; i++) {
            unsigned limitOther = i + otherSize < outSize ? otherSize : outSize - i;
            storageType k = 0;
            for (unsigned j = 0; j < limitOther; j++) {
                bigType t = static_cast<bigType>(this->s[i]) * o.s[j] + out->s[i+j] + k;
                out->s[i+j] = t;
                k = t >> storageSize;
            }
            out->s[i + limitOther] = k;
        }

        // r has unsigned product, correct if this or o are less than zero
        if (static_cast<signedType>(o.s[otherSize-1]) < 0 && otherSize < outSize) {
            MultiwordInteger<outSize, storageType> tmp(*this);
            tmp <<= (storageSize * otherSize);
            *out -= tmp;
        }
        if(static_cast<signedType>(this->s[size-1]) < 0 && size < outSize) {
            MultiwordInteger<outSize, storageType> tmp(o);
            tmp <<= storageSize * size;
            *out -= tmp;
        }
    }

    template<unsigned otherSize>
    constexpr
    void
    quotrem(MultiwordInteger<otherSize, storageType> divisor,
            MultiwordInteger<size,      storageType> &qotient,
            MultiwordInteger<otherSize, storageType> *remainder) const {
        unsigned divisor_length = otherSize, dividend_length = size;
        size_t divisor_nlz = divisor.is_negative() ? ((-divisor).leading_zeros()) : divisor.leading_zeros();

        divisor_length -= divisor_nlz / storageSize;
        if (otherSize > 1) {
            divisor_nlz %= storageSize;
        }
        dividend_length -= this->leading_zeros()/storageSize;

        if (dividend_length == 0) {
            qotient = int64_t(0);
            if (remainder) {
                *remainder = *this;
            }
            return;
        }

        /* Dividing by zero. Set the dividend as close to infinity as we can,
         * and the remainder to zero.
         */
        if (divisor_length == 0) {
            if (this->is_negative()) {
                qotient = minVal();
            } else {
                qotient = maxVal();
            }
            if (remainder) {
                *remainder = int64_t(0);
            }
            return;
        }

        bool negate_result = false;
        bool negative_remainder = false;
        MultiwordInteger<size+1, storageType> unsigned_dividend(*this);

        if (this->is_negative()) {
            negate_result = true;
            unsigned_dividend.negate();
        }

        if (divisor.is_negative()) {
            negate_result = !negate_result;
            negative_remainder = true;
            divisor.negate();
        }

        if (divisor_length < 2) {
            unsignedQuotrem(unsigned_dividend, divisor, qotient, dividend_length);
        } else {
            unsignedQuotrem(unsigned_dividend, divisor, qotient, dividend_length, divisor_length, divisor_nlz);
        }

        if (negate_result) {
            qotient.negate();
            if (unsigned_dividend) {
                --qotient;
            }
        }
        if (!remainder) {
            return;
        }

        if (divisor_length >= 2) {
            unsigned_dividend >>= divisor_nlz;
        }

        *remainder = unsigned_dividend;
        if (negate_result && unsigned_dividend) {
            *remainder -= MultiwordInteger<size+1, storageType>(divisor);
            if (!negative_remainder) {
                // Remainder's negative, we want it to be positive
                remainder->negate();
            }
        } else if (negative_remainder) {
            remainder->negate();
        }
    }

protected:
    /* Division by a divisor of length 1
     * The remainder is put into the dividend.
     */
    template<unsigned divisorSize>
    constexpr
    void
    unsignedQuotrem(MultiwordInteger<size+1,      storageType> &dividend,
                    MultiwordInteger<divisorSize, storageType> const &divisor,
                    MultiwordInteger<size,        storageType> &quotient,
                    unsigned dividendSize) const {
        storageType k = 0;
        bigType b = ((bigType)1)<<storageSize;
        for (unsigned i = size; i-- > dividendSize; ) {
            quotient.s[i] = 0;
        }
        while (dividendSize-- > 0) {
            bigType t = b*k;
            t += dividend.s[dividendSize];
            quotient.s[dividendSize] = t / divisor.s[0];
            k = t - quotient.s[dividendSize] * divisor.s[0];
        }

        dividend = k;
    }

    /* Division by a longer divisor
     * The remainder is put into the dividend.
     */
    template<unsigned divisorSize>
    constexpr
    void
    unsignedQuotrem(MultiwordInteger<size+1,      storageType> &dividend,
                    MultiwordInteger<divisorSize, storageType>  divisor,
                    MultiwordInteger<size,        storageType> &quotient,
                    unsigned dividend_length, unsigned divisor_length,
                    unsigned divisor_nlz) const {
        divisor <<= divisor_nlz;
        dividend <<= divisor_nlz;

        quotient = int64_t(0);

        for (int j = dividend_length - divisor_length; j >= 0; j--) {
            bigType b = 1UL<<storageSize;
            bigType p = 0;
            bigType qhat = (dividend.s[j+divisor_length] * b + dividend.s[j+divisor_length-1]) / divisor.s[divisor_length-1];
            bigType rhat = (dividend.s[j+divisor_length] * b + dividend.s[j+divisor_length-1]) - qhat * divisor.s[divisor_length-1];

            bool retry = false;
            do {
                retry = false;
                if (qhat >= b || qhat * divisor.s[divisor_length-2] > b*rhat + dividend.s[j+divisor_length-2]) {
                    qhat--;
                    rhat += divisor.s[divisor_length-1];
                    if (rhat < b) {
                        retry = true;
                    }
                }
            } while (retry);

            typename std::make_signed<bigType>::type k = 0;
            typename std::make_signed<bigType>::type t = 0;
            for (unsigned i = 0; i < divisor_length; i++) {
                p = qhat * divisor.s[i];
                t = dividend.s[i+j] - k - (p & ((1UL<<storageSize)-1));
                dividend.s[i+j] = t;
                k = (p >> storageSize) - (t >> storageSize);
            }
            t = dividend.s[j+divisor_length] - k;
            dividend.s[j+divisor_length] = t;

            quotient.s[j] = qhat;
            if (t < 0) {
                quotient.s[j]--;
                k = 0;
                for (unsigned i = 0; i < divisor_length; i++) {
                    t = dividend.s[i+j] + divisor.s[i] + k;
                    dividend.s[i+j] = t;
                    k = t >> storageSize;
                }
                dividend.s[j+divisor_length] += k;
            }
        }
    }

public:
    constexpr
    void
    negate() {
        bigType c = 1;
        for (unsigned i = 0; i < size; i++) {
            c += static_cast<storageType>(~s[i]);
            s[i] = c;
            c >>= storageSize;
        }
    }

    constexpr
    operator bool() const {
        for (unsigned i = 0; i < size; i++) {
            if (s[i]) {
                return true;
            }
        }
        return false;
    }

    explicit constexpr
    operator double() const {
        double r = 0;
        double m = pow(2., storageSize);
        double n = 1;

        if (this->is_positive()) {
            return -double(-*this);
        }

        bigType c = 1;

        for (unsigned i = 0; i < size; i++) {
            c += static_cast<storageType>(~s[i]);
            storageType u = c;
            r += n * u;
            n *= m;
            c >>= storageSize;
        }

        return -r;
    }

    constexpr
    unsigned
    leading_zeros() const {
        unsigned r = 0;
        for (unsigned i = size; i-- > 0; ) {
            if (s[i] == 0) {
                r += storageSize;
            } else {
                return r + nlz(s[i]);
            }
        }
        return r;
    }

    constexpr
    bool
    is_negative() const {
        return static_cast<signedType>(s[size-1]) < 0;
    }

    constexpr
    bool
    is_positive() const {
        if (static_cast<signedType>(s[size-1]) > 0) {
            return true;
        } else if (static_cast<signedType>(s[size-1]) < 0) {
            return false;
        }
        for (unsigned i = size - 1; i-- > 0; ) {
            if (s[i]) {
                return true;
            }
        }
        return false;
    }

    constexpr
    void
    fill_leading_bits(unsigned num) {
        unsigned i = size - 1;
        while (num >= storageSize) {
            s[i] = ~static_cast<storageType>(0);
            num -= storageSize;
            i--;
        }
        if (num) {
            s[i] |= (~static_cast<storageType>(0)) << (storageSize - num);
        }
    }

    // Return true if bit at position is set, starting from the LSB.
    bool bit(size_t position) const {
        size_t word = position % storageSize;
        position %= storageSize;
        return s[word] & (1<<position);
    }
};

template<int integerWidth, unsigned _fractionalWidth, typename backingStorageType = uint32_t>
class FixedPoint
{
    static_assert((integerWidth + _fractionalWidth + 1) % (sizeof(backingStorageType) * 8) == 0, "Integer and fractional width do not match size of storage type.");
public:
    typedef MultiwordInteger<(integerWidth+_fractionalWidth+1)/(sizeof(backingStorageType)*8), backingStorageType> StorageType;
protected:
    typedef FixedPoint<integerWidth,_fractionalWidth,backingStorageType> FP;
public:
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
        if (v > 1 << integerWidth) {
            this->v = StorageType::maxVal();
        } else if(v < -(1<<integerWidth)) {
            this->v = StorageType::minVal();
        } else {
            this->v = (int64_t)v;
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

    constexpr static
    unsigned
    fractionalWidth() {
        return _fractionalWidth;
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
        return T(FP(StorageType(int64_t(1))));
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
        r.v = (int64_t)v;
        return r;
    }

    constexpr
    FP
    operator+=(FP const &o) {
        StorageType z(int64_t(0));
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
};

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator+ (MultiwordInteger<size, storageType> left,
           MultiwordInteger<size, storageType> const &right) { return left += right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator- (MultiwordInteger<size, storageType> left,
           MultiwordInteger<size, storageType> const &right) { return left -= right; }

template<unsigned leftSize, unsigned rightSize, typename storageType = uint32_t>
constexpr MultiwordInteger<leftSize+rightSize, storageType>
operator* (MultiwordInteger<leftSize, storageType> left,
           MultiwordInteger<rightSize, storageType> const &right) { MultiwordInteger<leftSize+rightSize, storageType> out; left.mul(right, &out); return out; }

template<unsigned leftSize, unsigned rightSize, typename storageType = uint32_t>
constexpr MultiwordInteger<leftSize, storageType>
operator/ (MultiwordInteger<leftSize, storageType> left,
           MultiwordInteger<rightSize, storageType> const &right) { return left /= right; }

template<unsigned leftSize, unsigned rightSize, typename storageType = uint32_t>
constexpr MultiwordInteger<rightSize, storageType>
operator% (MultiwordInteger<leftSize, storageType> left,
           MultiwordInteger<rightSize, storageType> const &right) {
    MultiwordInteger<leftSize, storageType> q;
    MultiwordInteger<rightSize, storageType> r;
    left.quotrem(right, q, &r);
    return r; }

template<unsigned size, typename storageType = uint32_t>
constexpr bool
operator>(MultiwordInteger<size, storageType> const &left,
          MultiwordInteger<size, storageType> const &right) { return right < left; }

template<unsigned size, typename storageType = uint32_t>
constexpr bool
operator<=(MultiwordInteger<size, storageType> const &left,
           MultiwordInteger<size, storageType> const &right) { return !(left > right); }

template<unsigned size, typename storageType = uint32_t>
constexpr bool
operator>=(MultiwordInteger<size, storageType> const &left,
           MultiwordInteger<size, storageType> const &right) { return !(left < right); }

template<unsigned size, typename storageType = uint32_t>
constexpr bool
operator!=(MultiwordInteger<size, storageType> const &left,
           MultiwordInteger<size, storageType> const &right) { return !(left == right); }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator>>(MultiwordInteger<size, storageType> left, unsigned right) { return left >>= right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator<<(MultiwordInteger<size, storageType> left, unsigned right) { return left <<= right; }

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
}

#endif // FIXEDPOINT_HPP
