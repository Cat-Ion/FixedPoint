/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef MULTIWORDINTEGERCOMPARISONS_HPP
#define MULTIWORDINTEGERCOMPARISONS_HPP
#include "MultiwordInteger.hpp"

template<unsigned size, typename storageType> constexpr bool MultiwordInteger<size, storageType>::operator<(MultiwordInteger<size, storageType> const &o) const {
    if (static_cast<signedType>(s[size-1]) < static_cast<signedType>(o.s[size-1])) {
        return true;
    } else if(static_cast<signedType>(s[size-1]) > static_cast<signedType>(o.s[size-1])) {
        return false;
    }
    for (unsigned i = size - 1; i-- > 0; ) {
        if (s[i] < o.s[i]) {
            return true;
        } else if (s[i] > o.s[i]) {
            return false;
        }
    }
    return false;
}
template<unsigned size, typename storageType> constexpr bool MultiwordInteger<size, storageType>::operator==(MultiwordInteger<size, storageType> const &o) const {
    for (unsigned i = size; i-- > 0; ) {
        if (s[i] != o.s[i]) {
            return false;
        }
    }
    return true;
}

template<unsigned size, typename storageType>
constexpr bool
operator>(MultiwordInteger<size, storageType> const &left,
          MultiwordInteger<size, storageType> const &right) { return right < left; }

template<unsigned size, typename storageType>
constexpr bool
operator<=(MultiwordInteger<size, storageType> const &left,
           MultiwordInteger<size, storageType> const &right) { return !(left > right); }

template<unsigned size, typename storageType>
constexpr bool
operator>=(MultiwordInteger<size, storageType> const &left,
           MultiwordInteger<size, storageType> const &right) { return !(left < right); }

template<unsigned size, typename storageType>
constexpr bool
operator!=(MultiwordInteger<size, storageType> const &left,
           MultiwordInteger<size, storageType> const &right) { return !(left == right); }

#endif // MULTIWORDINTEGERCOMPARISONS_HPP
