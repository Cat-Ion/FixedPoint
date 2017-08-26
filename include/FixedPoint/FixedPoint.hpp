/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file at the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef FIXEDPOINT_HPP
#define FIXEDPOINT_HPP
#include "MultiwordInteger.hpp"

template<int _integerWidth, int _fractionalWidth, typename backingStorageType = uint32_t>
class FixedPoint
{
    static_assert((_integerWidth + _fractionalWidth + 1) % (sizeof(backingStorageType) * 8) == 0, "Integer and fractional width do not match size of storage type.");
    static_assert((_integerWidth + _fractionalWidth + 1) >= (sizeof(backingStorageType) * 8), "Integer and fractional width must sum up to at least the size of the backing storage type.");
public:
    typedef MultiwordInteger<(_integerWidth+_fractionalWidth+1)/(sizeof(backingStorageType)*8), backingStorageType> StorageType;

protected:
    typedef FixedPoint<_integerWidth,_fractionalWidth,backingStorageType> FP;

public:
    static const int fractionalWidth = _fractionalWidth;
    static const int integerWidth = _integerWidth;
    StorageType v;

    constexpr FixedPoint() : v(int64_t(0)) {}
    constexpr FixedPoint(double v);
    constexpr FixedPoint(int v);
    constexpr FixedPoint(FixedPoint const &o);
    constexpr FixedPoint(StorageType const &s);
    template<int oiw, int ofw, typename otherStorageType> constexpr FixedPoint(FixedPoint<oiw,ofw,otherStorageType> const &o);

    static constexpr FP _maxVal();
    static constexpr FP _minVal();
    static constexpr FP _smallestVal();

    static const FP maxVal;
    static const FP minVal;
    static const FP smallestVal;

    constexpr FP& operator+=(FP const &o);
    constexpr FP& operator-=(FP const &o);
    constexpr FP& operator*=(FP const &o);
    constexpr FP& operator*=(int64_t const &o);
    constexpr FP& operator/=(FP const &o);
    constexpr FP& operator/=(int64_t const &o);
    constexpr FP& operator%=(FP const &o);

    constexpr FP operator-() const;

    constexpr bool operator==(FP const &right) const;
    constexpr bool operator<(FP const &right) const;

    constexpr bool is_negative() const;
    constexpr bool is_positive() const;
    constexpr FP nabs() const;

    explicit constexpr operator double()  const;
    explicit constexpr operator int8_t()  const;
    explicit constexpr operator int16_t() const;
    explicit constexpr operator int32_t() const;
    explicit constexpr operator int64_t() const;

    constexpr void get_raw(uint8_t *out) const;
};
template<int integerWidth, int fractionalWidth, typename storageType>
constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>
FixedPoint<integerWidth, fractionalWidth, storageType>::maxVal = FixedPoint<integerWidth, fractionalWidth, storageType>::StorageType::_maxVal();

template<int integerWidth, int fractionalWidth, typename storageType>
constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>
FixedPoint<integerWidth, fractionalWidth, storageType>::minVal = FixedPoint<integerWidth, fractionalWidth, storageType>::StorageType::_minVal();

template<int integerWidth, int fractionalWidth, typename storageType>
constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>
FixedPoint<integerWidth, fractionalWidth, storageType>::smallestVal = typename FixedPoint<integerWidth, fractionalWidth, storageType>::StorageType(int64_t(1));


#include "FixedPointConstructors.hpp"
#include "FixedPointCasts.hpp"
#include "FixedPointArithmetics.hpp"
#include "FixedPointComparisons.hpp"
#include "FixedPointUtility.hpp"
#include "FixedPointStd.hpp"
#endif // FIXEDPOINT_HPP
