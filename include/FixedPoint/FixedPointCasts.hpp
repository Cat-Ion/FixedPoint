/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef FIXEDPOINTCASTS_HPP
#define FIXEDPOINTCASTS_HPP
#include "FixedPoint.hpp"
#include "FixedPointHelpers.hpp"

#define FLOAT_CAST(T) \
template<int integerWidth, int fractionalWidth, typename storageType> constexpr \
FixedPoint<integerWidth, fractionalWidth, storageType>::operator T() const \
{ \
    T dv = (T)(v); \
    T p2 = FixedPointHelpers::dipow<T>(2., -fractionalWidth); \
    T r = p2 * dv; \
    return r; \
}

FLOAT_CAST(float)
FLOAT_CAST(double)
FLOAT_CAST(long double)
#undef FLOAT_CAST

#define INT_CAST(T) \
template<int integerWidth, int fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::operator T()  const { return T (v>>fractionalWidth); }

INT_CAST(int8_t)
INT_CAST(int16_t)
INT_CAST(int32_t)
INT_CAST(int64_t)

#endif // FIXEDPOINTCASTS_HPP
