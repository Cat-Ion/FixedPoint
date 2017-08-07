#ifndef FIXEDPOINTARITHMETICS_HPP
#define FIXEDPOINTARITHMETICS_HPP
#include "FixedPoint.hpp"

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType>
constexpr FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>&
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::operator+=(FixedPoint<_integerWidth, _fractionalWidth, backingStorageType> const &o) {
    this->v += o.v;
    return *this;
}

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType>
constexpr FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>&
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::operator-=(FixedPoint<_integerWidth, _fractionalWidth, backingStorageType> const &o) {
    *this += -o;
    return *this;
}

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType>
constexpr FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>&
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::operator*=(FixedPoint<_integerWidth, _fractionalWidth, backingStorageType> const &o) {
    MultiwordInteger<StorageType::numWords*2, backingStorageType> nv;
    nv = this->v * o.v;
    this->v = nv >> fractionalWidth;
    if (fractionalWidth > 0 && nv.bit(fractionalWidth - 1)) {
        ++this->v;
    }
    return *this;
}

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType>
constexpr FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>&
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::operator*=(int64_t const &o) {
    this->v *= o;
    return *this;
}

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType>
constexpr FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>&
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::operator/=(FixedPoint<_integerWidth, _fractionalWidth, backingStorageType> const &o) {
    MultiwordInteger<StorageType::numWords+(fractionalWidth+StorageType::storageSize-1)/StorageType::storageSize, backingStorageType> nv(this->v);
    nv <<= fractionalWidth;
    nv /= o.v;
    v = nv;
    return *this;
}

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType>
constexpr FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>&
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::operator/=(int64_t const &o) {
    v /= o;
    return *this;
}

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType>
constexpr FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>&
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::operator%=(FixedPoint<_integerWidth, _fractionalWidth, backingStorageType> const &o) {
    v %= o.v;
    return *this;
}

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator+ (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return left += right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator- (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return left -= right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator* (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return left *= right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator* (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           int64_t const &right) { return left *= right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator/ (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return left /= right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator/ (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           int64_t const &right) { return left /= right; }

template<int integerWidth, unsigned minimumFractionalWidth, typename backingStorageType = uint32_t>
constexpr FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType>
operator% (FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> left,
           FixedPoint<integerWidth, minimumFractionalWidth, backingStorageType> const &right) { return left %= right; }

template<int _integerWidth, unsigned _fractionalWidth, typename backingStorageType>
constexpr FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>
FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>::operator-() const {
    return FixedPoint<_integerWidth, _fractionalWidth, backingStorageType>(-v);
}

#endif // FIXEDPOINTARITHMETICS_HPP
