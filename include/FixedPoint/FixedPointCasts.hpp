/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef FIXEDPOINTCASTS_HPP
#define FIXEDPOINTCASTS_HPP
#include "FixedPoint.hpp"
#include "FixedPointHelpers.hpp"

template<int integerWidth, int fractionalWidth, typename storageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::operator double() const
{
    double dv = double(v);
    double p2 = FixedPointHelpers::dipow(2., -fractionalWidth);
    double r = p2 * dv;
    return r;
}

template<int integerWidth, int fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::operator int8_t()  const { return int8_t (v>>fractionalWidth); }
template<int integerWidth, int fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::operator int16_t() const { return int16_t(v>>fractionalWidth); }
template<int integerWidth, int fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::operator int32_t() const { return int32_t(v>>fractionalWidth); }
template<int integerWidth, int fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::operator int64_t() const { return int64_t(v>>fractionalWidth); }

#endif // FIXEDPOINTCASTS_HPP
