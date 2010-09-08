#ifndef _CMPX_H_
#define _CMPX_H_

#include <stdio.h>
#include "mydefs.hpp"

class Cmpx {
// data members
private:
    fptype re;   // z = re + j*im
    fptype im;

public:
// constructors
    Cmpx();
    Cmpx(const fptype x, const fptype y);

// member functions
// get functions
    inline fptype get_re()  const  { return re;  }
    inline fptype get_im()  const  { return im;  }
    inline fptype get_mag() const  { return sqrt(re*re + im*im); }
    inline fptype get_ang() const  { return atan2(im, re); }

// set functions
    // Cmpx& set_re(const fptype x);
    // Cmpx& set_im(const fptype x);

// complex number operations
    inline Cmpx conjugate() const { return Cmpx(re, -im); }

// overloaded operators
    void operator+=(const Cmpx &z);
    void operator-=(const Cmpx &z);
    void operator*=(const Cmpx &z);
    void operator/=(const Cmpx &z);
    void operator*=(const fptype k);

    // Cmpx& operator+=(const Cmpx &z);
    // Cmpx& operator-=(const Cmpx &z);
    // Cmpx& operator*=(const Cmpx &z);
    // Cmpx& operator/=(const Cmpx &z);
    // Cmpx& operator*=(const fptype k);
    // inline Cmpx operator+(const Cmpx &z) const  { Cmpx z1 = *this; return z1 += z; }
    // inline Cmpx operator+(const Cmpx &z) const  { return Cmpx(*this) += z; }
    // inline Cmpx operator-(const Cmpx &z) const  { return Cmpx(*this) -= z; }
    // inline Cmpx operator*(const Cmpx &z) const  { return Cmpx(*this) *= z; }
    // inline Cmpx operator/(const Cmpx &z) const  { return Cmpx(*this) /= z; }
    // inline Cmpx operator*(const fptype k) const  { return Cmpx(*this) *= k; }

// print methods
    char* polar() const;
    char* cartesian() const;
};

// inline Cmpx operator*(const fptype k, const Cmpx z) {
    // return z * k;
// }

#endif // #ifndef  _CMPX_H_
