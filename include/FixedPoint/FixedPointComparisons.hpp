#ifndef FIXEDPOINTCOMPARISONS_HPP
#define FIXEDPOINTCOMPARISONS_HPP
#include "FixedPoint.hpp"

template<int _integerWidth, int _fractionalWidth, typename backingStorageType>
constexpr bool FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::operator==(
        FixedPoint<_integerWidth, _fractionalWidth, backingStorageType> const &right) const {
    return this->v == right.v;
}

template<int _integerWidth, int _fractionalWidth, typename backingStorageType>
constexpr bool FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::operator<(
        FixedPoint<_integerWidth, _fractionalWidth, backingStorageType> const &right) const {
    return this->v < right.v;
}

template<int integerWidth, int minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr bool operator> (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &left,
                          FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return right < left; }

template<int integerWidth, int minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr bool operator<=(FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &left,
                          FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return !(left > right); }

template<int integerWidth, int minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr bool operator>=(FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &left,
                          FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return !(left < right); }

template<int integerWidth, int minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr bool operator!=(FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &left,
                          FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return !(left == right); }

#endif // FIXEDPOINTCOMPARISONS_HPP
