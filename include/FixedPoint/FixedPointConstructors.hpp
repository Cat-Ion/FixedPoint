/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef FIXEDPOINTCONSTRUCTORS_HPP
#define FIXEDPOINTCONSTRUCTORS_HPP
#include "FixedPoint.hpp"
#include "FixedPointHelpers.hpp"

template<int integerWidth, int fractionalWidth, typename storageType>
template<typename T>
constexpr void FixedPoint<integerWidth, fractionalWidth, storageType>::construct_from_float(T v)
{
    if (v > T(_maxVal())) {
        this->v = StorageType::_maxVal();
    } else if (v < T(_minVal())) {
        this->v = StorageType::_minVal();
    } else {
        v *= FixedPointHelpers::dipow<T>(T(2.), fractionalWidth);
        this->v = v;
    }
}

template<int integerWidth, int fractionalWidth, typename storageType>
template<typename T>
constexpr void FixedPoint<integerWidth, fractionalWidth, storageType>::construct_from_int(T v)
{
    constexpr int shiftAmount = (integerWidth >= sizeof(T)*8 - 1) ? 0 : integerWidth;
    if (shiftAmount && v > (T(1) << shiftAmount)) {
        this->v = StorageType::_maxVal();
    } else if(shiftAmount && v < -(T(1)<<shiftAmount)) {
        this->v = StorageType::_minVal();
    } else {
        this->v = storageType(v);
        this->v <<= fractionalWidth;
    }
}

template<int integerWidth, int fractionalWidth, typename storageType>
template<typename T>
constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(T const &v) {
    constexpr bool is_floating_type = FixedPointHelpers::is_one_of<T, float, double, long double>::value;
    constexpr bool is_signed_int_type = FixedPointHelpers::is_one_of<T, char, short, int, long, long long>::value;
    static_assert(is_floating_type || is_signed_int_type, "T must be floating point or signed integer type");
    if constexpr(is_floating_type) {
        construct_from_float<T>(v);
    } else if constexpr(is_signed_int_type) {
        construct_from_int<T>(v);
    }
}

#define CONSTRUCTOR template<int integerWidth, int fractionalWidth, typename storageType> constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint

CONSTRUCTOR(StorageType const &s) : v(s) {}

CONSTRUCTOR(unsigned long long int v)
{
    *this = FixedPoint<integerWidth, fractionalWidth, storageType>((long long)v);
}

CONSTRUCTOR(FixedPoint const &o) : v(o.v) {}

#undef CONSTRUCTOR

template<int integerWidth, int fractionalWidth, typename storageType>
template<int oiw, int ofw, typename ost>
constexpr FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(FixedPoint<oiw,ofw,ost> const &o) {
    if (ofw > fractionalWidth) {
        typename FixedPoint<oiw, ofw, ost>::StorageType ov = o.v >> (ofw-fractionalWidth);
        this->v = ov;
        if (o.v.bit(ofw-fractionalWidth-1)) {
            ++this->v;
        }
    } else {
        this->v = o.v;
        if (ofw < fractionalWidth) {
            this->v <<= (fractionalWidth - ofw);
        }
    }
}

#endif // FIXEDPOINTCONSTRUCTORS_HPP
