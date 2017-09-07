/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef MULTIWORDINTEGERARITHMETICS_HPP
#define MULTIWORDINTEGERARITHMETICS_HPP
#include "MultiwordInteger.hpp"

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator+=(MultiwordInteger<size, storageType> const &o)
{
    storageType c = 0;
    for(unsigned i = 0; i < size; i++) {
        bigType t = this->s[i] + o.s[i] + c;
        this->s[i] = t;
        c = t >> storageSize;
    }
    return *this;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator-=(MultiwordInteger<size, storageType> const &o)
{
    bigType c = 0;
    for (unsigned i = 0; i < size; i++) {
        c += this->s[i];
        c -= o.s[i];
        this->s[i] = c;
        c = signedType(c>>storageSize);
    }
    return *this;
}
template<unsigned size, typename storageType> template<unsigned otherSize> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator*=(MultiwordInteger<otherSize, storageType> const &o) {
    MultiwordInteger<size, storageType> nv;
    this->mul<otherSize, size>(o, &nv);
    *this = nv;
    return *this;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator*=(int64_t const &o) {
    MultiwordInteger<size, storageType> nv(o);
    nv *= *this;
    *this = nv;
    return *this;
}
template<unsigned size, typename storageType> template<unsigned otherSize> constexpr MultiwordInteger<size, storageType>&
MultiwordInteger<size, storageType>::operator/=(MultiwordInteger<otherSize, storageType> const &o) {
    quotrem(o, *this, static_cast<MultiwordInteger<otherSize, storageType>*>(nullptr));
    return *this;
}

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>&
MultiwordInteger<size, storageType>::operator%=(MultiwordInteger<size, storageType> const &o) { *this = *this % o; return *this; }


template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>
operator+ (MultiwordInteger<size, storageType> left, MultiwordInteger<size, storageType> const &right) { return left += right; }

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>
operator- (MultiwordInteger<size, storageType> left, MultiwordInteger<size, storageType> const &right) { return left -= right; }

template<unsigned size, typename storageType, unsigned otherSize> constexpr MultiwordInteger<size+otherSize, storageType>
operator* (MultiwordInteger<size, storageType> left, MultiwordInteger<otherSize, storageType> const &right) {
    MultiwordInteger<size+otherSize, storageType> out;
    left.mul(right, &out);
    return out; }

template<unsigned size, typename storageType> constexpr MultiwordInteger<size+8/sizeof(storageType), storageType>
operator* (MultiwordInteger<size, storageType> left, int64_t const &right) {
    MultiwordInteger<size+8/sizeof(storageType), storageType> out(right);
    out *= *left;
    return out; }

template<unsigned size, typename storageType, unsigned otherSize> constexpr MultiwordInteger<size, storageType>
operator/ (MultiwordInteger<size, storageType> left, MultiwordInteger<otherSize, storageType> const &right) { return left /= right; }

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>
operator/ (MultiwordInteger<size, storageType> left, int64_t const &right) { return left /= right; }

template<unsigned size, typename storageType, unsigned otherSize> constexpr MultiwordInteger<otherSize, storageType>
operator% (MultiwordInteger<size, storageType> left, MultiwordInteger<otherSize, storageType> const &right) {
    MultiwordInteger<size, storageType> q;
    MultiwordInteger<otherSize, storageType> r;
    left.quotrem(right, q, &r);
    return r; }


template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>&
MultiwordInteger<size, storageType>::operator++() {
    bigType c = 1;
    for (unsigned i = 0; i < size && c; i++) {
        c += s[i];
        s[i] = c;
        c >>= storageSize;
    }
    return *this;
}

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>
MultiwordInteger<size, storageType>::operator++(int) {
    MultiwordInteger<size, storageType> r(*this);
    ++(*this);
    return r;
}

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>&
MultiwordInteger<size, storageType>::operator--() {
    bigType c = ~static_cast<bigType>(0);
    for (unsigned i = 0; i < size && c; i++) {
        c += s[i];
        s[i] = c;
        c >>= storageSize;
    }
    return *this;
}

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>
MultiwordInteger<size, storageType>::operator--(int) {
    MultiwordInteger<size, storageType> r(*this);
    --(*this);
    return r;
}


template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>
MultiwordInteger<size, storageType>::operator- () const {
    MultiwordInteger<size, storageType> r(*this);
    r.negate();
    return r;
}


template<unsigned size, typename storageType>
template<unsigned otherSize, unsigned outSize>
constexpr void MultiwordInteger<size, storageType>::mul (
        MultiwordInteger<otherSize, storageType> const &o,
        MultiwordInteger<outSize,   storageType>       *out) const{
    *out = storageType(0);

    unsigned limitThis = size > outSize ? outSize : size;
    for (unsigned i = 0; i < limitThis; i++) {
        unsigned limitOther = i + otherSize < outSize ? otherSize : outSize - i;
        storageType k = 0;
        for (unsigned j = 0; j < limitOther; j++) {
            bigType t = static_cast<bigType>(this->s[i]) * o.s[j] + out->s[i+j] + k;
            out->s[i+j] = t;
            k = t >> storageSize;
        }
        if (i + limitOther < outSize) {
            out->s[i + limitOther] = k;
        }
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

template<unsigned size, typename storageType> template<unsigned otherSize> constexpr void
MultiwordInteger<size, storageType>::quotrem(
        MultiwordInteger<otherSize, storageType> divisor,
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
        qotient = storageType(0);
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
            qotient = _minVal();
        } else {
            qotient = _maxVal();
        }
        if (remainder) {
            *remainder = storageType(0);
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

template<unsigned size, typename storageType> template<unsigned divisorSize> constexpr void
MultiwordInteger<size, storageType>::unsignedQuotrem(
        MultiwordInteger<size+1,      storageType> &dividend,
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

template<unsigned size, typename storageType> template<unsigned divisorSize> constexpr void
MultiwordInteger<size, storageType>::unsignedQuotrem(
        MultiwordInteger<size+1,      storageType> &dividend,
        MultiwordInteger<divisorSize, storageType>  divisor,
        MultiwordInteger<size,        storageType> &quotient,
        unsigned dividend_length, unsigned divisor_length,
        unsigned divisor_nlz) const {
    divisor <<= divisor_nlz;
    dividend <<= divisor_nlz;

    quotient = storageType(0);

    for (int j = dividend_length - divisor_length; j >= 0; j--) {
        bigType b = bigType(1)<<storageSize;
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
            t = dividend.s[i+j] - k - (p & ((bigType(1)<<storageSize)-1));
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

#endif // MULTIWORDINTEGERARITHMETICS_HPP
