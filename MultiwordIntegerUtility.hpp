#ifndef MULTIWORDINTEGERUTILITY_HPP
#define MULTIWORDINTEGERUTILITY_HPP
#include "FixedPoint.hpp"

template<unsigned size, typename storageType> constexpr
void MultiwordInteger<size, storageType>::negate() {
    bigType c = 1;
    for (unsigned i = 0; i < size; i++) {
        c += static_cast<storageType>(~s[i]);
        s[i] = c;
        c >>= storageSize;
    }
}

template<unsigned size, typename storageType> constexpr
void MultiwordInteger<size, storageType>::fill_leading_bits(unsigned num) {
    unsigned i = size - 1;
    while (num >= storageSize) {
        s[i] = ~static_cast<storageType>(0);
        num -= storageSize;
        i--;
    }
    if (num) {
        s[i] |= (~static_cast<storageType>(0)) << (storageSize - num);
    }
}

template<unsigned size, typename storageType> constexpr
bool MultiwordInteger<size, storageType>::is_negative() const {
    return static_cast<signedType>(s[size-1]) < 0;
}

template<unsigned size, typename storageType> constexpr
bool MultiwordInteger<size, storageType>::is_positive() const {
    if (static_cast<signedType>(s[size-1]) > 0) {
        return true;
    } else if (static_cast<signedType>(s[size-1]) < 0) {
        return false;
    }
    for (unsigned i = size - 1; i-- > 0; ) {
        if (s[i]) {
            return true;
        }
    }
    return false;
}

template<unsigned size, typename storageType> constexpr
bool MultiwordInteger<size, storageType>::bit(size_t position) const {
    size_t word = position / storageSize;
    position %= storageSize;
    return s[word] & (1<<position);
}
#endif // MULTIWORDINTEGERUTILITY_HPP
