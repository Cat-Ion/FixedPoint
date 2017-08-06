#ifndef MULTIWORDINTEGERBITWISE_HPP
#define MULTIWORDINTEGERBITWISE_HPP
#include "FixedPoint.hpp"

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator &=(MultiwordInteger<size, storageType> const &o) {
    for (unsigned i = 0; i < size; i++) {
        s[i] &= o.s[i];
    }
    return *this;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator |=(MultiwordInteger<size, storageType> const &o) {
    for (unsigned i = 0; i < size; i++) {
        s[i] |= o.s[i];
    }
    return *this;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator ^=(MultiwordInteger<size, storageType> const &o) {
    for (unsigned i = 0; i < size; i++) {
        s[i] ^= o.s[i];
    }
    return *this;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>  MultiwordInteger<size, storageType>::operator ~() const {
    MultiwordInteger<size, storageType> r;
    for (unsigned i = 0; i < size; i++) {
        r.s[i] = ~s[i];
    }
    return r;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator <<=(unsigned n) {
    size_t w = n / storageSize;
    n %= storageSize;

    if (w >= size) {
        *this = storageType(0);
        return *this;
    }

    if (n == 0) {
        for (unsigned i = size; i-- > w; ) {
            s[i] = s[i-w];
        }
        for (unsigned i = 0; i < w; i++) {
            s[i] = 0;
        }
        return *this;
    }

    unsigned i = size;
    while (i-- > (1+w)) {
        this->s[i] = (this->s[i - w] << n)
                + (this->s[i - w - 1] >> (storageSize - n));
    }
    this->s[i] = (this->s[i - w] << n) ;
    while (i-- > 0) {
        this->s[i] = 0;
    }

    return *this;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator >>=(unsigned n) {
    size_t w = n / storageSize;
    n %= storageSize;

    bool adjust_leading_bits = this->is_negative();

    if (w >= size) {
        if (adjust_leading_bits) {
            for (unsigned i = 0; i < size; i++) {
                s[i] = ~static_cast<storageType>(0);
            }
        } else {
            *this = storageType(0);
        }
        return *this;
    }

    if (n == 0) {
        for (unsigned i = 0; i < size-w; i++) {
            s[i] = s[i+w];
        }
        for (unsigned i = size-w; i < size; i++) {
            s[i] = 0;
        }
        if (adjust_leading_bits) {
            this->fill_leading_bits(w*storageSize);
        }
        return *this;
    }

    unsigned i = 0;
    while (i < (size - w - 1)) {
        this->s[i] = (this->s[i + w] >> n) + (this->s[i + w + 1] << (storageSize - n));
        i++;
    }
    this->s[i] = this->s[i + w] >> n;
    while (++i < size) {
        this->s[i] = 0;
    }

    if (adjust_leading_bits) {
        this->fill_leading_bits(w * storageSize + n);
    }
    return *this;
}

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator>>(MultiwordInteger<size, storageType> left, unsigned right) { return left >>= right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator<<(MultiwordInteger<size, storageType> left, unsigned right) { return left <<= right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator&(MultiwordInteger<size, storageType> left, MultiwordInteger<size, storageType> const &right) { return left &= right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator|(MultiwordInteger<size, storageType> left, MultiwordInteger<size, storageType> const &right) { return left |= right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator^(MultiwordInteger<size, storageType> left, MultiwordInteger<size, storageType> const &right) { return left ^= right; }

#endif // MULTIWORDINTEGERBITWISE_HPP
