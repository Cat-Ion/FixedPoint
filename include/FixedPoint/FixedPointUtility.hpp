/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef FIXEDPOINTUTILITY_HPP
#define FIXEDPOINTUTILITY_HPP
#include "FixedPoint.hpp"

template<int _integerWidth, int _fractionalWidth, typename backingStorageType> constexpr
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::_maxVal() {
    return StorageType::_maxVal();
}

template<int _integerWidth, int _fractionalWidth, typename backingStorageType> constexpr
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::_minVal() {
    return StorageType::_minVal();
}

template<int _integerWidth, int _fractionalWidth, typename backingStorageType> constexpr
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::_smallestVal() {
    return StorageType::_smallestVal();
}

template<int _integerWidth, int _fractionalWidth, typename backingStorageType> constexpr
bool
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::is_negative() const {
  return v.is_negative();
}

template<int _integerWidth, int _fractionalWidth, typename backingStorageType> constexpr
bool
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::is_positive() const {
  return v.is_positive();
}

template<int _integerWidth, int _fractionalWidth, typename backingStorageType> constexpr
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::nabs() const {
    return v.is_positive() ? -*this : *this;
}

template<int _integerWidth, int _fractionalWidth, typename backingStorageType> constexpr
void
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::get_raw(uint8_t *out) const {
    v.get_raw(out);
}

#endif // FIXEDPOINTUTILITY_HPP
