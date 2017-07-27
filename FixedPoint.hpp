#ifndef FIXEDPOINT_HPP
#define FIXEDPOINT_HPP
#define NDEBUG
#include <assert.h>
#undef NDEBUG
#include <cmath>
#include <cstring>
#include <string>
#include <type_traits>

template<typename T> class make_bigger { };
template<> class make_bigger<uint8_t > { public: typedef uint16_t type; };
template<> class make_bigger<uint16_t> { public: typedef uint32_t type; };
template<> class make_bigger<uint32_t> { public: typedef uint64_t type; };

template<unsigned size, typename _storageType = uint32_t>
class MultiwordInteger
{
public:
    typedef _storageType storageType;
    static const constexpr size_t numWords = size;
    static const constexpr size_t storageSize = sizeof(storageType) * 8;
    _storageType s[size];

protected:
    typedef typename std::make_signed<_storageType>::type signedType;
    typedef typename make_bigger<_storageType>::type bigType;

public:

    MultiwordInteger() {}
    template<unsigned otherSize>
    constexpr
    MultiwordInteger(MultiwordInteger<otherSize, storageType> const &o) {
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

    constexpr
    MultiwordInteger(storageType v) {
        s[0] = v;
        for (unsigned i = 1; i < size; i++) {
            s[i] = 0;
        }
    }

    constexpr
    MultiwordInteger(int64_t v) {
        unsigned i = 0;
        while (v && i < size) {
            s[i] = v & ((1UL<<storageSize) - 1);
            v >>= storageSize;
            i++;
        }
        while (i < size) {
            s[i] = 0;
            i++;
        }
    }

    constexpr
    MultiwordInteger(double v) {
        double sig = significand(v);
        int lg = ilogb(v);
        if (lg < 52) {
            int64_t i = sig*(1UL<<lg);
            *this = i;
        } else {
            *this = int64_t(sig * (1UL<<52));
            *this <<= lg - 52;
        }
    }

    static constexpr
    MultiwordInteger<size, storageType>
    maxVal() {
        MultiwordInteger<size, storageType> r;
        for (unsigned i = 0; i < size-1; i++) {
            r.s[0] = (bigType(1)<<(storageSize)) - 1;;
        }
        r.s[size-1] = (1U<<(storageSize - 1)) - 1;
        return r;
    }

    static constexpr
    MultiwordInteger<size, storageType>
    minVal() {
        MultiwordInteger<size, storageType> r;
        for (unsigned i = 0; i < size-1; i++) {
            r.s[0] = 0;
        }
        r.s[size-1] = 1U<<(storageSize-1);
        return r;
    }

    constexpr
    MultiwordInteger<size, storageType>
    operator+=(MultiwordInteger<size, storageType> const &o) {
        storageType c = 0;
        unsigned i = size;
        while (i-- > 0) {
            bigType t = this->s[i] + o.s[i] + c;
            this->s[i] = t;
            c = t >> storageSize;
        }
        return *this;
    }

    constexpr
    MultiwordInteger<size, storageType>
    operator-=(MultiwordInteger<size, storageType> const &o) {
        bigType c = 0;
        for (unsigned i = 0; i < size; i++) {
            c += this->s[i];
            c -= o.s[i];
            this->s[i] = c;
            c = (signedType)(c>>storageSize);
        }
        return *this;
    }

    template<unsigned otherSize>
    constexpr
    MultiwordInteger<size, storageType>
    operator*=(MultiwordInteger<otherSize, storageType> const &o) {
        MultiwordInteger<size> nv;
        mul<otherSize, size>(o, &nv);
        *this = nv;
    }

    template<unsigned otherSize>
    constexpr
    MultiwordInteger<size, storageType>
    operator/=(MultiwordInteger<otherSize, storageType> const &o) {
        *this = MultiwordInteger<size, storageType>(*this / o);
        return *this;
    }

    constexpr
    MultiwordInteger<size, storageType>
    operator%=(MultiwordInteger<size, storageType> const &o) {
        *this = *this % o;
    }

    constexpr
    MultiwordInteger<size, storageType>
    operator<<=(size_t n) {
        size_t w = n / storageSize;
        n %= storageSize;

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

    constexpr
    MultiwordInteger<size, storageType>
    operator>>=(size_t n) {
        size_t w = n / storageSize;
        n %= storageSize;

        if (n == 0) {
            for (unsigned i = 0; i < size-w; i++) {
                s[i] = s[i+w];
            }
            for (unsigned i = size-w; i < size; i++) {
                s[i] = 0;
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

        return *this;
    }

    constexpr
    MultiwordInteger<size, storageType>
    operator-() const {
        MultiwordInteger<size, storageType> r(int64_t(0));
        r -= *this;
        return r;
    }

    constexpr
    bool
    operator<(MultiwordInteger<size, storageType> const &o) const {
        if (static_cast<signedType>(s[size-1]) < static_cast<signedType>(o.s[size-1])) {
            return true;
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

    constexpr
    bool
    operator==(MultiwordInteger<size, storageType> const &o) const {
        for (unsigned i = size; i-- > 0; ) {
            if (s[i] != o.s[i]) {
                return false;
            }
        }
        return true;
    }

    template<unsigned otherSize, unsigned outSize>
    constexpr
    void
    mul(MultiwordInteger<otherSize, storageType> const &o,
        MultiwordInteger<outSize,   storageType>       *out) const {
        *out = (int64_t) 0;

        unsigned limitThis = size > outSize ? outSize : size;
        for (unsigned i = 0; i < limitThis; i++) {
            unsigned limitOther = i + otherSize < outSize ? otherSize : outSize - i;
            storageType k = 0;
            for (unsigned j = 0; j < limitOther; j++) {
                bigType t = static_cast<bigType>(this->s[i]) * o.s[j] + out->s[i+j] + k;
                out->s[i+j] = t;
                k = t >> storageSize;
            }
            out->s[i + limitOther] = k;
        }

        // r has unsigned product, correct if this or o are less than zero
        if (static_cast<signedType>(o.s[otherSize-1]) < 0 && otherSize < outSize) {
            MultiwordInteger<outSize, storageType> tmp(*this);
            tmp <<= (storageSize * otherSize);
            *out -= tmp;
        }
        if(static_cast<signedType>(this->s[size-1]) < 0 && size < outSize) {
            MultiwordInteger<outSize, storageType> tmp(o);
            tmp <<= storageSize * size;
            *out -= tmp;
        }
    }

    template<unsigned otherSize>
    constexpr
    void
    quotrem(MultiwordInteger<otherSize, storageType> const &o_,
            MultiwordInteger<size,      storageType>       *q,
            MultiwordInteger<otherSize, storageType>       *r) const {
        unsigned n = otherSize, m = size;
        size_t s = o_.is_negative() ? ((-o_).nlz()) : o_.nlz();
        n -= s / storageSize;
        s %= storageSize;
        m -= this->nlz()/storageSize;

        for(unsigned i = size; i-- > m; ) {
            q->s[i] = 0;
        }

        if (n < 2) {
            if (n == 0) {
                if (this->is_negative()) {
                    *q = minVal();
                } else {
                    *q = maxVal();
                }
                return;
            }

            bool negate = false;
            if (m > 0 && static_cast<signedType>(o_.s[n-1] ^ this->s[m-1]) < 0) {
                negate = true;
            }

            storageType o = o_.is_negative() ? -static_cast<signedType>(o_.s[0]) : o_.s[0];
            MultiwordInteger<size, storageType> p;

            p = this->is_negative() ? -*this : *this;

            int k = 0;
            unsigned j = m;
            bigType b = ((bigType)1)<<storageSize;
            while (j-- > 0) {
                bigType t = b*k;
                t += p.s[j];
                t /= o;
                q->s[j] = (k * b + p.s[j]) / o;
                k = (k * b + p.s[j]) - q->s[j] * o;
            }
            if (r) {
                r->s[0] = k;
            }
            if (negate) {
                q->negate();
                if (r) {
                    r->s[0] = -k;
                }
            }
        } else {

            MultiwordInteger<size+1, storageType> nu;
            MultiwordInteger<otherSize, storageType> nv;

            if (this->is_negative()) {
                nu = -*this;
            } else {
                nu = *this;
            }
            if (o_.is_negative()) {
                nv = -o_;
            } else {
                nv = o_;
            }

            nu <<= s;
            nv <<= s;

            for (int j = m - n; j >= 0; j--) {
                bigType b = 1UL<<storageSize;
                bigType p = 0;
                bigType qhat = (nu.s[j+n] * b + nu.s[j+n-1]) / nv.s[n-1];
                bigType rhat = (nu.s[j+n] * b + nu.s[j+n-1]) - qhat * nv.s[n-1];

                bool retry = false;
                do {
                    retry = false;
                    if (qhat >= b || qhat * nv.s[n-2] > b*rhat + nu.s[j+n-2]) {
                        qhat--;
                        rhat += nv.s[n-1];
                        if (rhat < b) {
                            retry = true;
                        }
                    }
                } while (retry);

                int k = 0;
                bigType t = 0;
                for (unsigned i = 0; i < n; i++) {
                    p = qhat * nv.s[i];
                    t = nu.s[i+j] - k - (p & ((1UL<<storageSize)-1));
                    nu.s[i+j] = t;
                    k = (p >> storageSize) - (t >> storageSize);
                }
                t = nu.s[j+n] - k;
                nu.s[j+n] = t;

                q->s[j] = qhat;
                if (t < 0) {
                    q->s[j]--;
                    k = 0;
                    for (unsigned i = 0; i < n; i++) {
                        t = nu.s[i+j] + nv.s[i] + k;
                        nu.s[i+j] = t;
                        k = t >> storageSize;
                    }
                    nu.s[j+n] += k;
                }
            }

            if (r) {
                *r = nu;
                *r >>= s;
            }
        }
    }

    constexpr
    void
    negate() {
        bigType c = 1;
        for (unsigned i = 0; i < size; i++) {
            c += ~s[i];
            s[i] = c;
            c >>= storageSize;
        }
    }

    explicit constexpr
    operator double() const {
        double r = 0;
        double m = pow(2., storageSize);
        double n = 1;

        if (this->is_positive()) {
            return -double(-*this);
        }

        bigType c = 1;

        for (unsigned i = 0; i < size; i++) {
            bigType t = ~s[i] + c;
            storageType u = t;
            r += n * u;
            n *= m;
            c = t >> storageSize;
        }

        return -r;
    }

    constexpr
    size_t
    nlz() const {
        size_t r = 0;
        for (unsigned i = size; i-- > 0; ) {
            if (s[i] == 0) {
                r += storageSize;
            } else {
                unsigned x = s[i];
                if (x <= 0x0000FFFF) { r += 16; x <<= 16; }
                if (x <= 0x00FFFFFF) { r +=  8; x <<=  8; }
                if (x <= 0x0FFFFFFF) { r +=  4; x <<=  4; }
                if (x <= 0x3FFFFFFF) { r +=  2; x <<=  2; }
                if (x <= 0x7FFFFFFF) { r +=  1; x <<=  1; }
                r -= 32 - storageSize;
                return r;
            }
        }
        return r;
    }

    constexpr
    bool
    is_negative() const {
        return static_cast<signedType>(s[size-1]) < 0;
    }

    constexpr
    bool
    is_positive() const {
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
};

template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
class FixedPoint
{
protected:
    typedef FixedPoint<iw,mfw,backingStorageType> FP;
    typedef MultiwordInteger<(iw+mfw+sizeof(backingStorageType)*8-1)/(sizeof(backingStorageType)*8), backingStorageType> StorageType;
    static constexpr size_t fw = StorageType::storageSize * StorageType::numWords - iw - 1;
public:
    StorageType v;

    template<typename T>
    constexpr static
    T
    maxVal()
    {
        return T(FP(StorageType::maxVal()));
    }

    template<typename T>
    constexpr static
    T
    minVal()
    {
        return T(FP(StorageType::minVal()));
    }

    template<typename T>
    constexpr static
    T
    smallestVal()
    {
        return T(FP(StorageType(int64_t(1))));
    }

    constexpr
    FixedPoint()
    {
    }

    constexpr
    FixedPoint(double v)
    {
        if (v > maxVal<double>()) {
            v = maxVal<double>();
        } else if (v < minVal<double>()) {
            v = minVal<double>();
        }
        this->v = v * pow(2., (double)fw);
    }

    constexpr
    FixedPoint(int v)
    {
        if (v > 1 << iw) {
            this->v = StorageType::maxVal();
        } else if(v < -(1<<iw)) {
            this->v = StorageType::minVal();
        } else {
            this->v = (int64_t)v;
            this->v <<= fw;
        }
    }

    constexpr
    FixedPoint(FixedPoint const &o) {
        this->v = o.v;
    }

    constexpr
    FixedPoint(StorageType const &s) {
        this->v = s;
    }

    explicit constexpr
    operator double() const
    {
        double dv = v;
        double p2 = pow(2., -(double) fw);
        double r = p2 * dv;
        return r;
    }

    constexpr static
    FP
    from_int(int v)
    {
        FP r;
        r.v = (int64_t)v;
        return r;
    }

    constexpr
    FP
    operator+=(FP const &o) {
        StorageType z(int64_t(0));
        if (o.v.is_negative() && v.is_negative()) {
            StorageType d(StorageType::minVal());
            d -= v;
            assert (d <= o.v);
        } else if(o.v.is_positive() && v.is_positive()) {
            StorageType d(StorageType::maxVal());
            d -= v;
            assert (d >= o.v);
        }
        this->v += o.v;
        return *this;
    }

    constexpr
    FP
    operator-=(FP const &o) {
        *this += -o;
        return *this;
    }

    constexpr
    FP
    operator*=(FP const &o) {
        assert((double)*this * (double)o <= maxVal<double>());
        assert((double)*this * (double)o >= minVal<double>());
        MultiwordInteger<StorageType::numWords*2, backingStorageType> nv;
        nv = this->v * o.v;
        this->v = nv >> fw;
        return *this;
    }

    constexpr
    FP
    operator*=(int64_t const &o) {
        assert((double)*this * (double)o <= maxVal<double>());
        assert((double)*this * (double)o >= minVal<double>());
        this->v *= o;
        return *this;
    }

    constexpr
    FP&
    operator/=(FP const &o) {
        assert((double)*this / (double)o <= maxVal<double>());
        assert((double)*this / (double)o >= minVal<double>());
        MultiwordInteger<StorageType::numWords*2-iw/StorageType::storageSize, backingStorageType> nv(this->v);
        nv <<= fw;
        nv /= o.v;
        v = nv;
        return *this;
    }

    constexpr
    FP&
    operator/=(int64_t const &o) {
        assert((double)*this / (double)o <= maxVal<double>());
        assert((double)*this / (double)o >= minVal<double>());
        v /= o;
        return *this;
    }

    constexpr
    FP&
    operator=(FP const &o)
    {
      v = o.v;
      return *this;
    }

    friend constexpr
    FP
    operator-(FP right) {
        right.v = -right.v;
        return right;
    }

    friend constexpr
    bool
    operator==(FP const &left, FP const &right) {
        return left.v == right.v;
    }

    friend constexpr
    bool
    operator<(FP const &left, FP const &right) {
        return left.v < right.v;
    }

    constexpr
    bool
    is_negative() const {
      return v.is_negative();
    }

    constexpr
    bool
    is_positive() const {
      return v.is_positive();
    }

    constexpr
    FP
    nabs() const {
        return v.is_positive() ? -*this : *this;
    }
};

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator+ (MultiwordInteger<size, storageType> left,
           MultiwordInteger<size, storageType> const &right) { return left += right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator- (MultiwordInteger<size, storageType> left,
           MultiwordInteger<size, storageType> const &right) { return left -= right; }

template<unsigned leftSize, unsigned rightSize, typename storageType = uint32_t>
constexpr MultiwordInteger<leftSize+rightSize, storageType>
operator* (MultiwordInteger<leftSize, storageType> left,
           MultiwordInteger<rightSize, storageType> const &right) { MultiwordInteger<leftSize+rightSize, storageType> out; left.mul(right, &out); return out; }

template<unsigned leftSize, unsigned rightSize, typename storageType = uint32_t>
constexpr MultiwordInteger<leftSize, storageType>
operator/ (MultiwordInteger<leftSize, storageType> left,
           MultiwordInteger<rightSize, storageType> const &right) { return left /= right; }

template<unsigned leftSize, unsigned rightSize, typename storageType = uint32_t>
constexpr MultiwordInteger<rightSize, storageType>
operator% (MultiwordInteger<leftSize, storageType> left,
           MultiwordInteger<leftSize, storageType> const &right) {
    MultiwordInteger<leftSize, storageType> q;
    MultiwordInteger<rightSize, storageType> r;
    left.quotrem(right, &q, &r);
    return r; }

template<unsigned size, typename storageType = uint32_t>
constexpr bool
operator>(MultiwordInteger<size, storageType> const &left,
          MultiwordInteger<size, storageType> const &right) { return right < left; }

template<unsigned size, typename storageType = uint32_t>
constexpr bool
operator<=(MultiwordInteger<size, storageType> const &left,
           MultiwordInteger<size, storageType> const &right) { return !(left > right); }

template<unsigned size, typename storageType = uint32_t>
constexpr bool
operator>=(MultiwordInteger<size, storageType> const &left,
           MultiwordInteger<size, storageType> const &right) { return !(left < right); }

template<unsigned size, typename storageType = uint32_t>
constexpr bool
operator!=(MultiwordInteger<size, storageType> const &left,
           MultiwordInteger<size, storageType> const &right) { return !(left == right); }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator>>(MultiwordInteger<size, storageType> left, size_t right) { return left >>= right; }

template<unsigned size, typename storageType = uint32_t>
constexpr MultiwordInteger<size, storageType>
operator<<(MultiwordInteger<size, storageType> left, size_t right) { return left <<= right; }

template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
constexpr FixedPoint<iw, mfw, backingStorageType>
operator+ (FixedPoint<iw, mfw, backingStorageType> left,
           FixedPoint<iw, mfw, backingStorageType> const &right) { return left += right; }

template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
constexpr FixedPoint<iw, mfw, backingStorageType>
operator- (FixedPoint<iw, mfw, backingStorageType> left,
           FixedPoint<iw, mfw, backingStorageType> const &right) { return left -= right; }

template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
constexpr FixedPoint<iw, mfw, backingStorageType>
operator* (FixedPoint<iw, mfw, backingStorageType> left,
           FixedPoint<iw, mfw, backingStorageType> const &right) { return left *= right; }

template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
constexpr FixedPoint<iw, mfw, backingStorageType>
operator/ (FixedPoint<iw, mfw, backingStorageType> left,
           FixedPoint<iw, mfw, backingStorageType> const &right) { return left /= right; }

template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
constexpr bool operator> (FixedPoint<iw, mfw, backingStorageType> const &left,
                          FixedPoint<iw, mfw, backingStorageType> const &right) { return right < left; }

template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
constexpr bool operator<=(FixedPoint<iw, mfw, backingStorageType> const &left,
                          FixedPoint<iw, mfw, backingStorageType> const &right) { return !(left > right); }

template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
constexpr bool operator>=(FixedPoint<iw, mfw, backingStorageType> const &left,
                          FixedPoint<iw, mfw, backingStorageType> const &right) { return !(left < right); }

template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
constexpr bool operator!=(FixedPoint<iw, mfw, backingStorageType> const &left,
                          FixedPoint<iw, mfw, backingStorageType> const &right) { return !(left == right); }


namespace std {
    template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw, mfw, backingStorageType>
    abs(FixedPoint<iw,mfw,backingStorageType> const &y) {
        if (y.is_negative()) {
            if (y == FixedPoint<iw, mfw,backingStorageType>::template minVal<FixedPoint<iw, mfw,backingStorageType>>()) {
                return FixedPoint<iw, mfw, backingStorageType>::template maxVal<FixedPoint<iw, mfw, backingStorageType>>();
            } else {
                return -y;
            }
        } else {
            return y;
        }
    }

    template<int iw, unsigned mfw, typename backingStorageType = uint32_t>
    constexpr FixedPoint<iw,mfw, backingStorageType>
    atan2(FixedPoint<iw,mfw,backingStorageType> const &y, FixedPoint<iw,mfw,backingStorageType> const &x)
    {
        FixedPoint<iw,mfw,backingStorageType>
                mx = std::min(y.nabs(), x.nabs()),
                mn = std::max(y.nabs(), x.nabs()),
                a = mn / mx,
                s = a*a,
                r(-0.0464964749);

        r = r * s + FixedPoint<iw,mfw,backingStorageType>(0.15931422);
        r = r * s + FixedPoint<iw,mfw,backingStorageType>(-0.327622764);
        r = r * s * a + a;

        if (y.nabs() < x.nabs()) {
            r = FixedPoint<iw,mfw,backingStorageType>(M_PI_2) - r;
        }
        if (x.is_negative()) {
            r = FixedPoint<iw,mfw,backingStorageType>(M_PI) - r;
        }
        if (y.is_negative()) {
            r = -r;
        }

        return r;
    }
}

#endif // FIXEDPOINT_HPP
