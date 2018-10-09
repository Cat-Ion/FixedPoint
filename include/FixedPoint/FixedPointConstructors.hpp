/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef FIXEDPOINTCONSTRUCTORS_HPP
#define FIXEDPOINTCONSTRUCTORS_HPP
#include "FixedPoint.hpp"
#include "FixedPointHelpers.hpp"

template<int integerWidth, int fractionalWidth, typename storageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(double v)
{
    if (v > double(_maxVal())) {
        this->v = StorageType::_maxVal();
    } else if (v < double(_minVal())) {
        this->v = StorageType::_minVal();
    } else {
        v *= FixedPointHelpers::dipow(2., fractionalWidth);
        this->v = v;
    }
}

template<int integerWidth, int fractionalWidth, typename storageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(int v)
{
    constexpr int shiftAmount = (integerWidth >= sizeof(int)*8 - 1) ? 0 : integerWidth;
    if (shiftAmount && v > (1 << shiftAmount)) {
        this->v = StorageType::_maxVal();
    } else if(shiftAmount && v < -(1<<shiftAmount)) {
        this->v = StorageType::_minVal();
    } else {
        this->v = storageType(v);
        this->v <<= fractionalWidth;
    }
}

template<int integerWidth, int fractionalWidth, typename storageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(FixedPoint const &o) : v(o.v) {}

template<int integerWidth, int fractionalWidth, typename storageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(StorageType const &s) : v(s) {}

template<int integerWidth, int fractionalWidth, typename storageType>
template<int oiw, int ofw, typename otherStorageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(FixedPoint<oiw,ofw,otherStorageType> const &o) {
    if (ofw > fractionalWidth) {
        typename FixedPoint<oiw, ofw, otherStorageType>::StorageType ov = o.v >> (ofw-fractionalWidth);
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
