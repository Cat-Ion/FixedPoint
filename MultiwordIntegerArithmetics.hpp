#ifndef MULTIWORDINTEGERARITHMETICS_HPP
#define MULTIWORDINTEGERARITHMETICS_HPP
#include "FixedPoint.hpp"

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator+=(MultiwordInteger<size, storageType> const &o)
{
    storageType c = 0;
    for(unsigned i = 0; i < size; i++) {
        bigType t = this->s[i] + o.s[i] + c;
        this->s[i] = t;
        c = t >> storageSize;
    }
    return *this;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator-=(MultiwordInteger<size, storageType> const &o)
{
    bigType c = 0;
    for (unsigned i = 0; i < size; i++) {
        c += this->s[i];
        c -= o.s[i];
        this->s[i] = c;
        c = signedType(c>>storageSize);
    }
    return *this;
}
template<unsigned size, typename storageType> template<unsigned otherSize> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator*=(MultiwordInteger<otherSize, storageType> const &o) {
    MultiwordInteger<size, storageType> nv;
    this->mul<otherSize, size>(o, &nv);
    *this = nv;
    return *this;
}
template<unsigned size, typename storageType> template<unsigned otherSize> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator/=(MultiwordInteger<otherSize, storageType> const &o) {
    quotrem(o, *this, static_cast<MultiwordInteger<otherSize, storageType>*>(nullptr));
    return *this;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator%=(MultiwordInteger<size, storageType> const &o) { *this = *this % o; return *this; }

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType> operator+ (MultiwordInteger<size, storageType> left, MultiwordInteger<size, storageType> const &right) { return left += right; }
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType> operator- (MultiwordInteger<size, storageType> left, MultiwordInteger<size, storageType> const &right) { return left -= right; }
template<unsigned size, typename storageType, unsigned otherSize> constexpr MultiwordInteger<size+otherSize, storageType> operator* (MultiwordInteger<size, storageType> left, MultiwordInteger<otherSize, storageType> const &right) { MultiwordInteger<size+otherSize, storageType> out; left.mul(right, &out); return out; }
template<unsigned size, typename storageType, unsigned otherSize> constexpr MultiwordInteger<size, storageType> operator/ (MultiwordInteger<size, storageType> left, MultiwordInteger<otherSize, storageType> const &right) { return left /= right; }
template<unsigned size, typename storageType, unsigned otherSize> constexpr MultiwordInteger<otherSize, storageType> operator% (MultiwordInteger<size, storageType> left, MultiwordInteger<otherSize, storageType> const &right) {
    MultiwordInteger<size, storageType> q;
    MultiwordInteger<otherSize, storageType> r;
    left.quotrem(right, q, &r);
    return r; }

template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator++() {
    bigType c = 1;
    for (unsigned i = 0; i < size && c; i++) {
        c += s[i];
        s[i] = c;
        c >>= storageSize;
    }
    return *this;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType> MultiwordInteger<size, storageType>::operator++(int) {
    MultiwordInteger<size, storageType> r(*this);
    ++(*this);
    return r;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator--() {
    bigType c = ~static_cast<bigType>(0);
    for (unsigned i = 0; i < size && c; i++) {
        c += s[i];
        s[i] = c;
        c >>= storageSize;
    }
    return *this;
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>& MultiwordInteger<size, storageType>::operator--(int) {
    MultiwordInteger<size, storageType> r(*this);
    *this--;
    return r;
}

#endif // MULTIWORDINTEGERARITHMETICS_HPP
