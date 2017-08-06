#ifndef FIXEDPOINTCONSTRUCTORS_HPP
#define FIXEDPOINTCONSTRUCTORS_HPP
#include "FixedPoint.hpp"

template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(double v)
{
    if (v > double(maxVal)) {
        this->v = StorageType::_maxVal();
    } else if (v < double(minVal)) {
        this->v = StorageType::_minVal();
    } else {
        this->v = v * pow(2., double(fractionalWidth));
    }
}

template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(int v)
{
    if (v > 1 << integerWidth) {
        this->v = StorageType::_maxVal();
    } else if(v < -(1<<integerWidth)) {
        this->v = StorageType::_minVal();
    } else {
        this->v = storageType(v);
        this->v <<= fractionalWidth;
    }
}

template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(FixedPoint const &o) : v(o.v) {}

template<int integerWidth, unsigned fractionalWidth, typename storageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(StorageType const &s) : v(s) {}

template<int integerWidth, unsigned fractionalWidth, typename storageType>
template<int oiw, unsigned ofw, typename otherStorageType> constexpr
FixedPoint<integerWidth, fractionalWidth, storageType>::FixedPoint(FixedPoint<oiw,ofw,otherStorageType> const &o) {
    if (ofw > fractionalWidth) {
        typename FixedPoint<oiw, ofw, otherStorageType>::StorageType ov = o.v >> (ofw-fractionalWidth);
        this->v = ov;
    } else {
        this->v = o.v;
        this->v <<= (fractionalWidth - ofw);
    }
}
#endif // FIXEDPOINTCONSTRUCTORS_HPP
