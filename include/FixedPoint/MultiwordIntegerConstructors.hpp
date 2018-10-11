/* Copyright 2017 Valentin Ochs <a at 0au dot de>.
 * See the LICENSE file in the top-level directory and at
 * https://github.com/Cat-Ion/FixedPoint/blob/master/LICENSE. */
#ifndef MULTIWORDINTEGERCONSTRUCTORS_HPP
#define MULTIWORDINTEGERCONSTRUCTORS_HPP
#include <cmath>
#include "MultiwordInteger.hpp"
#include "FixedPointHelpers.hpp"

template<unsigned size, typename storageType> template<unsigned otherSize> constexpr MultiwordInteger<size, storageType>::MultiwordInteger(MultiwordInteger<otherSize, storageType> const &o) : s{0} {
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
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::MultiwordInteger(storageType const &v) : s{v} {
    if (static_cast<signedType>(v) < 0) {
        for (unsigned i = 1; i < size; i++) {
            s[i] = static_cast<signedType>(-1);
        }
    }
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::MultiwordInteger(int8_t v) : s{0} {
    *this = MultiwordInteger<size, std::make_unsigned<decltype(v)>::type>(v);
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::MultiwordInteger(int16_t v) : s{0} {
    *this = MultiwordInteger<size, std::make_unsigned<decltype(v)>::type>(v);
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::MultiwordInteger(int32_t v) : s{0} {
    *this = MultiwordInteger<size, std::make_unsigned<decltype(v)>::type>(v);
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::MultiwordInteger(int64_t v) : s{0} {
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
template<unsigned size, typename storageType> 
template<typename U, typename EN>
constexpr MultiwordInteger<size, storageType>::MultiwordInteger(unsigned long long int v, typename EN::type *) {
    *this = MultiwordInteger<size, storageType>((long long)v);
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::MultiwordInteger(float v) : s{0} {
    if (v != 0.) {
        if (v > INT32_MIN && v < INT32_MAX) {
            int32_t i = int32_t(v);
            *this = i;
        } else {
            int lg = FixedPointHelpers::ilogb(v);
            v *= FixedPointHelpers::dipow<float>(2, 30-lg);
            *this = int64_t(v);
            *this <<= lg - 30;
        }
    }
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::MultiwordInteger(double v) : s{0} {
    if (v != 0.) {
        if (v > double(INT64_MIN) && v < double(INT64_MAX)) {
            int64_t i = int64_t(v);
            *this = i;
        } else {
            int lg = FixedPointHelpers::ilogb(v);
            v *= FixedPointHelpers::dipow<double>(2, 62-lg);
            *this = MultiwordInteger(int64_t(v));
            *this <<= lg - 62;
        }
    }
}
template<unsigned size, typename storageType> constexpr MultiwordInteger<size, storageType>::MultiwordInteger(long double v) : s{0} {
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
template<unsigned size, typename storageType> template<unsigned otherSize, typename otherStorageType> constexpr MultiwordInteger<size, storageType>::MultiwordInteger(MultiwordInteger<otherSize, otherStorageType> const &o) : s{0} {
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

#endif // MULTIWORDINTEGERCONSTRUCTORS_HPP
