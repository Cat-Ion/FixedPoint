/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file at the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef FIXEDPOINT_HPP
#define FIXEDPOINT_HPP
#include "MultiwordInteger.hpp"

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

    constexpr FixedPoint() : v(int64_t(0)) {}
    constexpr FixedPoint(double v);
    constexpr FixedPoint(int v);
    constexpr FixedPoint(FixedPoint const &o);
    constexpr FixedPoint(StorageType const &s);
    template<int oiw, unsigned ofw, typename otherStorageType> constexpr FixedPoint(FixedPoint<oiw,ofw,otherStorageType> const &o);

    constexpr static FP maxVal = FP(StorageType::maxVal);
    constexpr static FP minVal = FP(StorageType::minVal);
    constexpr static FP smallestVal = FP(StorageType(int64_t(1)));

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
};
template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType> FixedPoint<integerWidth, fractionalWidth, storageType>::maxVal;
template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType> FixedPoint<integerWidth, fractionalWidth, storageType>::minVal;
template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType> FixedPoint<integerWidth, fractionalWidth, storageType>::smallestVal;


#include "FixedPointConstructors.hpp"
#include "FixedPointCasts.hpp"
#include "FixedPointArithmetics.hpp"
#include "FixedPointComparisons.hpp"
#include "FixedPointUtility.hpp"
#include "FixedPointStd.hpp"
#endif // FIXEDPOINT_HPP
