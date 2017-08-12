#ifndef MULTIWORDINTEGER_HPP
#define MULTIWORDINTEGER_HPP
#include <algorithm>
#include <stdint.h>
#include "FixedPointHelpers.hpp"

template<unsigned size, typename _storageType = uint32_t>
class MultiwordInteger
{
    template<unsigned otherSize, typename otherStorageType>
    friend class MultiwordInteger;

public:
    typedef _storageType storageType;
    static const constexpr size_t numWords = size;
    static const constexpr size_t storageSize = sizeof(storageType) * 8;

protected:
    typedef typename std::make_signed<_storageType>::type signedType;
    typedef typename FixedPointHelpers::make_bigger<_storageType>::type bigType;
    _storageType s[size];

public:
    static constexpr MultiwordInteger<size, storageType> _maxVal();
    static constexpr MultiwordInteger<size, storageType> _minVal();
    static const MultiwordInteger<size, storageType> maxVal;
    static const MultiwordInteger<size, storageType> minVal;

    constexpr MultiwordInteger() : s{0} {}
    template<unsigned otherSize> constexpr MultiwordInteger(MultiwordInteger<otherSize, storageType> const &o);
    constexpr MultiwordInteger(storageType const &v);
    constexpr MultiwordInteger(int64_t v);
    constexpr MultiwordInteger(double v);
    template<unsigned otherSize, typename otherStorageType> constexpr MultiwordInteger(MultiwordInteger<otherSize, otherStorageType> const &o);

    constexpr MultiwordInteger<size, storageType>& operator+=(MultiwordInteger<size, storageType> const &o);
    constexpr MultiwordInteger<size, storageType>& operator-=(MultiwordInteger<size, storageType> const &o);
    template<unsigned otherSize> constexpr MultiwordInteger<size, storageType>& operator*=(MultiwordInteger<otherSize, storageType> const &o);
    constexpr MultiwordInteger<size, storageType>& operator*=(int64_t const &o);
    template<unsigned otherSize> constexpr MultiwordInteger<size, storageType>& operator/=(MultiwordInteger<otherSize, storageType> const &o);
    constexpr MultiwordInteger<size, storageType>& operator/=(int64_t const &o);
    constexpr MultiwordInteger<size, storageType>& operator%=(MultiwordInteger<size, storageType> const &o);

    constexpr MultiwordInteger<size, storageType>& operator++();
    constexpr MultiwordInteger<size, storageType>  operator++(int);
    constexpr MultiwordInteger<size, storageType>& operator--();
    constexpr MultiwordInteger<size, storageType>  operator--(int);

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

template<unsigned size, typename storageType> const constexpr
MultiwordInteger<size, storageType> MultiwordInteger<size, storageType>::maxVal = MultiwordInteger<size, storageType>::_maxVal();

template<unsigned size, typename storageType> const constexpr
MultiwordInteger<size, storageType> MultiwordInteger<size, storageType>::minVal = MultiwordInteger<size, storageType>::_minVal();

#include "MultiwordIntegerConstructors.hpp"
#include "MultiwordIntegerArithmetics.hpp"
#include "MultiwordIntegerBitwise.hpp"
#include "MultiwordIntegerCasts.hpp"
#include "MultiwordIntegerComparisons.hpp"
#include "MultiwordIntegerUtility.hpp"
#include "MultiwordIntegerSpecialization.hpp"
#endif // MULTIWORDINTEGER_HPP
