#ifndef FIXEDPOINTUTILITY_HPP
#define FIXEDPOINTUTILITY_HPP
#include "FixedPoint.hpp"

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType> constexpr
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::_maxVal() {
    return StorageType::_maxVal();
}

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType> constexpr
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::_minVal() {
    return StorageType::_minVal();
}

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType> constexpr
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::_smallestVal() {
    return StorageType::_smallestVal();
}

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType> constexpr
bool
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::is_negative() const {
  return v.is_negative();
}

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType> constexpr
bool
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::is_positive() const {
  return v.is_positive();
}

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType> constexpr
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::nabs() const {
    return v.is_positive() ? -*this : *this;
}

#endif // FIXEDPOINTUTILITY_HPP
