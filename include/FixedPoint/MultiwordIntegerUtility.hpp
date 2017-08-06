#ifndef MULTIWORDINTEGERUTILITY_HPP
#define MULTIWORDINTEGERUTILITY_HPP
#include "MultiwordInteger.hpp"

template<unsigned size, typename storageType> constexpr
MultiwordInteger<size, storageType> MultiwordInteger<size, storageType>::_maxVal() {
    MultiwordInteger<size, storageType> r;
    for (unsigned i = 0; i < size-1; i++) {
        r.s[i] = (bigType(1)<<(storageSize)) - 1;;
    }
    r.s[size-1] = (1U<<(storageSize - 1)) - 1;
    return r;
}

template<unsigned size, typename storageType> constexpr
MultiwordInteger<size, storageType> MultiwordInteger<size, storageType>::_minVal() {
    MultiwordInteger<size, storageType> r;
    for (unsigned i = 0; i < size-1; i++) {
        r.s[i] = 0;
    }
    r.s[size-1] = 1U<<(storageSize-1);
    return r;
}

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
unsigned MultiwordInteger<size, storageType>::leading_zeros() const {
    unsigned r = 0;
    for (unsigned i = size; i-- > 0; ) {
        if (s[i] == 0) {
            r += storageSize;
        } else {
            return r + nlz(s[i]);
        }
    }
    return r;
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
