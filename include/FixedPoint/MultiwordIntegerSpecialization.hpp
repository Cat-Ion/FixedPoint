#ifndef MULTIWORDINTEGERSPECIALIZATION_HPP
#define MULTIWORDINTEGERSPECIALIZATION_HPP
#include "MultiwordInteger.hpp"
#include "FixedPointHelpers.hpp"

template<typename _storageType>
class MultiwordInteger<1, _storageType> {
    template<unsigned otherSize, typename otherStorageType>
    friend class MultiwordInteger;
    static const constexpr unsigned size = 1;

public:
    typedef _storageType storageType;
    static const constexpr size_t numWords = 1;
    static const constexpr size_t storageSize = sizeof(storageType) * 8;

protected:
    typedef typename std::make_signed<_storageType>::type signedType;
    typedef typename FixedPointHelpers::make_bigger<_storageType>::type bigType;
    _storageType s[1];

public:
    static const MultiwordInteger<1, storageType> maxVal;
    static const MultiwordInteger<1, storageType> minVal;
    static constexpr MultiwordInteger<1, storageType> _maxVal() { return MultiwordInteger<size, storageType>(storageType((storageType(1)<<(storageSize-1)) - 1)); }
    static constexpr MultiwordInteger<1, storageType> _minVal() { return MultiwordInteger<size, storageType>(storageType((storageType(1)<<(storageSize-1)))); }

    constexpr MultiwordInteger() : s{0} {}
    constexpr MultiwordInteger(storageType const &v) : s{v} {}
    constexpr MultiwordInteger(int64_t v) : s{storageType(signedType(v))} {}
    constexpr MultiwordInteger(double v) : s{storageType(signedType(v))} {}

    template<unsigned otherSize, typename otherStorageType> constexpr
    MultiwordInteger(MultiwordInteger<otherSize, otherStorageType> const &o) : s{0}{
        for (unsigned i = 0; i < otherSize && i * 8 * sizeof(otherStorageType) < storageSize; i++) {
            s[0] |= o.s[i] << (i * sizeof(otherStorageType) * 8);
        }
        if (o.is_negative() && o.storageSize < storageSize) {
            fill_leading_bits(unsigned(storageSize - o.storageSize));
        }
    }

    constexpr MultiwordInteger<size, storageType>& operator+=(MultiwordInteger<size, storageType> const &o) { s[0] += o.s[0]; return *this;}
    constexpr MultiwordInteger<size, storageType>& operator-=(MultiwordInteger<size, storageType> const &o) { s[0] -= o.s[0]; return *this;}
    template<unsigned otherSize> constexpr MultiwordInteger<size, storageType>& operator*=(MultiwordInteger<otherSize, storageType> const &o) { s[0] *= o.s[0]; return *this; }
    constexpr MultiwordInteger<size, storageType>& operator*=(int64_t const &o) { s[0] *= storageType(o); return *this; }
    template<unsigned otherSize> constexpr MultiwordInteger<size, storageType>& operator/=(MultiwordInteger<otherSize, storageType> const &o) { quotrem(o, *this, static_cast<MultiwordInteger<otherSize, storageType>*>(0)); return *this; }
    constexpr MultiwordInteger<size, storageType>& operator/=(int64_t const &o) { MultiwordInteger<size, storageType> div(o); quotrem(div, *this, static_cast<MultiwordInteger<size, storageType>*>(0)); return *this; }
    constexpr MultiwordInteger<size, storageType>& operator%=(MultiwordInteger<size, storageType> const &o) {
        *this = *this % o;
        return *this; }
    constexpr MultiwordInteger<size, storageType> operator%(MultiwordInteger<size, storageType> const &o) {
        MultiwordInteger<size, storageType> q, r;
        quotrem(o, q, &r);
        return r;
    }

    constexpr MultiwordInteger<size, storageType>& operator++() { s[0]++; return *this; }
    constexpr MultiwordInteger<size, storageType>  operator++(int) { MultiwordInteger<size, storageType> r(*this); ++*this; return r; }
    constexpr MultiwordInteger<size, storageType>& operator--() { s[0]--; return *this; }
    constexpr MultiwordInteger<size, storageType>  operator--(int) { MultiwordInteger<size, storageType> r(*this); --*this; return r; }

    constexpr MultiwordInteger<size, storageType> operator-() const { return MultiwordInteger<size, storageType>(storageType(-s[0])); }

    constexpr MultiwordInteger<size, storageType>& operator &=(MultiwordInteger<size, storageType> const &o) { s[0] &= o.s[0]; return *this; }
    constexpr MultiwordInteger<size, storageType>& operator |=(MultiwordInteger<size, storageType> const &o) { s[0] |= o.s[0]; return *this; }
    constexpr MultiwordInteger<size, storageType>& operator ^=(MultiwordInteger<size, storageType> const &o) { s[0] ^= o.s[0]; return *this; }
    constexpr MultiwordInteger<size, storageType>  operator ~() const { return MultiwordInteger<size, storageType>(storageType(~s[0])); }
    constexpr MultiwordInteger<size, storageType>& operator <<=(unsigned n) {
        if (n >= storageSize) {
            s[0] = 0;
        } else {
            s[0]<<=n;
        }
        return *this;
    }
    constexpr MultiwordInteger<size, storageType>& operator >>=(unsigned n) { s[0] = storageType(signedType(s[0])>>n); return *this; }

    constexpr bool operator<(MultiwordInteger<size, storageType> const &o) const { return signedType(s[0]) < signedType(o.s[0]); }
    constexpr bool operator==(MultiwordInteger<size, storageType> const &o) const { return s[0] == o.s[0]; }

    constexpr void negate() { s[0] = ~s[0] + 1; }
    constexpr void fill_leading_bits(unsigned num) { s[0] |= ~storageType(0) << (storageSize - num); }

    constexpr unsigned leading_zeros() const { return FixedPointHelpers::nlz(s[0]); }
    constexpr bool     is_negative() const { return signedType(s[0]) < 0; }
    constexpr bool     is_positive() const { return signedType(s[0]) > 0; }
    // Return true if bit at position is set, starting from the LSB.
    constexpr bool bit(size_t position) const { return (s[0] >> position) & 1; }

    explicit constexpr operator bool() const    { return signedType(s[0]); }
    explicit constexpr operator int8_t() const  { return signedType(s[0]); }
    explicit constexpr operator int16_t() const { return signedType(s[0]); }
    explicit constexpr operator int32_t() const { return signedType(s[0]); }
    explicit constexpr operator int64_t() const { return signedType(s[0]); }
    explicit constexpr operator double() const  { return signedType(s[0]); }

    template<unsigned otherSize, unsigned outSize> constexpr void
    mul(MultiwordInteger<otherSize, storageType> const &o,
        MultiwordInteger<outSize,   storageType>       *out) const {
        if (otherSize == 1) {
            typename std::make_signed<bigType>::type c = signedType(s[0]);
            c *= signedType(o.s[0]);
            *out = MultiwordInteger<outSize, storageType>(storageType(c));
            if (outSize > 1) {
                out->s[1] = c >> storageSize;
            }
        } else {
            o.mul(*this, out);
        }
    }

    template<unsigned otherSize> constexpr void
    quotrem(MultiwordInteger<otherSize, storageType> divisor,
            MultiwordInteger<size,      storageType> &quotient,
            MultiwordInteger<otherSize, storageType> *remainder) const {

        if (divisor == MultiwordInteger<otherSize, storageType>(storageType(0))) {
            if (signedType(s[0]) < 0) {
                quotient = _minVal();
            } else {
                quotient = _maxVal();
            }
            if (remainder) {
                *remainder = storageType(0);
            }
            return;
        }

        if (otherSize > 1 && divisor.leading_zeros() <= (divisor.numWords * divisor.storageSize - storageSize)) {
            quotient = storageType(0);
            if (remainder) {
                *remainder = s[0];
            }
        } else {
            signedType s0 = s[0];
            quotient = MultiwordInteger<size, storageType>(storageType(signedType(s[0])/signedType(divisor.s[0])));
            if (signedType(s0 ^ divisor.s[0]) < 0 && (s0 % signedType(divisor.s[0])) != 0) {
                --quotient;
                if (remainder) {
                    *remainder = MultiwordInteger<otherSize, storageType>(storageType(s0 - quotient.s[0]*divisor.s[0]));
                }
            } else {
                if (remainder) {
                    *remainder = MultiwordInteger<otherSize, storageType>(storageType(s0 % signedType(divisor.s[0])));
                }
            }
        }
    }

    MultiwordInteger<size, storageType> nabs() const {
        if (is_positive()) {
            return -*this;
        } else {
            return *this;
        }
    }
};

template<typename storageType> const constexpr
MultiwordInteger<1, storageType> MultiwordInteger<1, storageType>::maxVal = MultiwordInteger<1, storageType>::_maxVal();

template<typename storageType> const constexpr
MultiwordInteger<1, storageType> MultiwordInteger<1, storageType>::minVal = MultiwordInteger<1, storageType>::_minVal();

#endif // MULTIWORDINTEGERSPECIALIZATION_HPP
