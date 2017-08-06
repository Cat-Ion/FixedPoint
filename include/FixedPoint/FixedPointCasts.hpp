#ifndef FIXEDPOINTCASTS_HPP
#define FIXEDPOINTCASTS_HPP
#include "FixedPoint.hpp"

template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::operator double() const
{
    double dv = double(v);
    double p2 = pow(2., -double(fractionalWidth));
    double r = p2 * dv;
    return r;
}

template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::operator int8_t()  const { return int8_t (v>>fractionalWidth); }
template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::operator int16_t() const { return int16_t(v>>fractionalWidth); }
template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::operator int32_t() const { return int32_t(v>>fractionalWidth); }
template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::operator int64_t() const { return int64_t(v>>fractionalWidth); }

#endif // FIXEDPOINTCASTS_HPP
