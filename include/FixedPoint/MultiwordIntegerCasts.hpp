/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef MULTIWORDINTEGERCASTS_HPP
#define MULTIWORDINTEGERCASTS_HPP
#include "MultiwordInteger.hpp"
#include "FixedPointHelpers.hpp"

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::operator bool() const {
    for (unsigned i = 0; i < size; i++) {
        if (s[i]) {
            return true;
        }
    }
    return false;
}

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::operator int8_t() const {
    return s[0] & 0xFF;
}

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::operator int16_t() const {
    int16_t r = 0;
    for (unsigned i = 0, j = 0; i < 16; i += storageSize) {
        int16_t v = s[j++];
        r |= v << i;
    }
    return r;
}

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::operator int32_t() const {
    int32_t r = 0;
    for (unsigned i = 0, j = 0; i < 32; i += storageSize) {
        int32_t v = s[j++];
        r |= v << i;
    }
    return r;
}

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::operator int64_t() const {
    int64_t r = 0;
    for (unsigned i = 0, j = 0; i < 64; i += storageSize) {
        int64_t v = s[j++];
        r |= v << i;
    }
    return r;
}

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::operator double() const {
    double r = 0;
    double m = FixedPointHelpers::dipow(2., storageSize);
    double n = 1;

    if (this->is_positive()) {
        return -double(-*this);
    }

    bigType c = 1;

    for (unsigned i = 0; i < size; i++) {
        c += static_cast<storageType>(~s[i]);
        storageType u = c;
        r += n * u;
        n *= m;
        c >>= storageSize;
    }

    return -r;
}

#endif // MULTIWORDINTEGERCASTS_HPP
