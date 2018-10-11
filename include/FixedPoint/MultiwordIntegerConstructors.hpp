/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef MULTIWORDINTEGERCONSTRUCTORS_HPP
#define MULTIWORDINTEGERCONSTRUCTORS_HPP
#include <cmath>
#include "MultiwordInteger.hpp"
#include "FixedPointHelpers.hpp"

template<unsigned size, typename storageType>
template<unsigned otherSize>
constexpr MultiwordInteger<size, storageType>::MultiwordInteger(MultiwordInteger<otherSize, storageType> const &o) : s{0} {
    unsigned num = size > otherSize ? otherSize : size;
    for (unsigned i = 0; i < num; i++) {
        s[i] = o.s[i];
    }
    if (static_cast<signedType>(o.s[otherSize-1]) < 0) {
        for (unsigned i = num; i < size; i++) {
            s[i] = ~static_cast<storageType>(0);
        }
    } else {
        for (unsigned i = num; i < size; i++) {
            s[i] = 0;
        }
    }
}

template<unsigned size, typename storageType>
constexpr MultiwordInteger<size, storageType>::MultiwordInteger(storageType const &v) : s{v} {
    if (static_cast<signedType>(v) < 0) {
        for (unsigned i = 1; i < size; i++) {
            s[i] = static_cast<signedType>(-1);
        }
    }
}

template<unsigned size, typename storageType>
template<unsigned otherSize, typename otherStorageType>
constexpr MultiwordInteger<size, storageType>::MultiwordInteger(MultiwordInteger<otherSize, otherStorageType> const &o) : s{0} {
    static_assert((sizeof(storageType) % sizeof(otherStorageType)) == 0
                  || (sizeof(otherStorageType) % sizeof(storageType)) == 0,
                  "Types must fit into each other without remainder.");
    if (sizeof(otherStorageType) < sizeof(storageType)) {
        unsigned shiftNum = sizeof (storageType) / sizeof (otherStorageType);
        unsigned shiftWidth = sizeof (otherStorageType) * 8;
        unsigned i = 0;
        for (i = 0; i < size && i*shiftNum < otherSize; i++) {
            s[i] = 0;
            unsigned start = shiftNum;
            if (i*shiftNum + start > otherSize) {
                start = otherSize - i*shiftNum;
                if (o.is_negative()) {
                    s[i] |= (~static_cast<storageType>(0)) << (start * shiftWidth);
                }
            }
            for (unsigned j = start; j-- > 0; ) {
                s[i] |= o.s[i*shiftNum + j] << (shiftWidth * j);
            }
        }
        if (o.is_negative()) {
            while (i < size) {
                s[i] = ~static_cast<storageType>(0);
                i++;
            }
        } else {
            while (i < size) {
                s[i] = 0;
                i++;
            }
        }
    } else {
        unsigned shiftNum = sizeof (otherStorageType) / sizeof (storageType);
        unsigned shiftWidth = sizeof (storageType) * 8;
        unsigned i = 0, j = 0;
        for (i = j = 0; i < size && j < otherSize; j++) {
            otherStorageType v = o.s[j];
            for (unsigned k = 0; k < shiftNum && i < size; k++, i++) {
                s[i] = v;
                v >>= shiftWidth;
            }
        }
        if (o.is_negative()) {
            while (i < size) {
                s[i] = ~static_cast<storageType>(0);
                i++;
            }
        } else {
            while (i < size) {
                s[i] = 0;
                i++;
            }
        }
    }
}

template<unsigned size, typename storageType>
template<typename T>
constexpr MultiwordInteger<size, storageType>::MultiwordInteger(T v) : s{0} {
    static_assert(FixedPointHelpers::is_one_of<T, int8_t, int16_t, int32_t, int64_t, float, double, long double, unsigned long long>::value, "T must be a signed int or floating point type, or unsigned long long");
    if constexpr (FixedPointHelpers::is_one_of<T, int8_t, int16_t, int32_t>::value) {
        *this = MultiwordInteger<size, typename std::make_unsigned<T>::type>(v);
    }
    else if constexpr(FixedPointHelpers::is_one_of<T, int64_t>::value) {
        uint64_t uv = v;
        unsigned i = 0;
        for (; i < size && uv; i++) {
            s[i] = uv;
            uv >>= storageSize;
        }
        for (; i < size; i++) {
            s[i] = 0;
        }
        if(v < 0) {
            fill_leading_bits(leading_zeros());
        }
    }
    else if constexpr(FixedPointHelpers::is_one_of<T, float, double>::value) {
        if (v != 0.) {
            using IntType = typename std::conditional<std::is_same<T, float>::value, int32_t, int64_t>::type;
            if (v > std::numeric_limits<IntType>::min() && v < std::numeric_limits<IntType>::max()) {
                IntType i = IntType(v);
                *this = i;
            } else {
                int lg = FixedPointHelpers::ilogb(v);
                v *= FixedPointHelpers::dipow<T>(2, sizeof(IntType)*8-2-lg);
                *this = MultiwordInteger(int64_t(v));
                *this <<= lg + 2 - sizeof(IntType)*8;
            }
        }
    }
    else if constexpr(FixedPointHelpers::is_one_of<T, long double>::value) {
        if (v != 0.) {
            if (v > INT64_MIN && v < INT64_MAX) {
                int64_t i = int64_t(v);
                *this = i;
            } else {
                while (std::abs(v) >= 1) {
                    double tv = v;
                    MultiwordInteger<size, storageType> tmwi = tv;
                    *this += tmwi;
                    v -= tv;
                }
            }
        }
    }
    else if constexpr(FixedPointHelpers::is_one_of<T, unsigned long long>::value) {
        *this = MultiwordInteger<size, storageType>((long long)v);
    }
}

#endif // MULTIWORDINTEGERCONSTRUCTORS_HPP
