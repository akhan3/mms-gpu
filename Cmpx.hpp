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
    HOSTDEVICE Cmpx();
    HOSTDEVICE Cmpx(const fptype x, const fptype y);

// member functions
// get functions
    HOSTDEVICE inline fptype get_re()  const  { return re;  }
    HOSTDEVICE inline fptype get_im()  const  { return im;  }
    HOSTDEVICE inline fptype get_mag() const  { return sqrt(re*re + im*im); }
    HOSTDEVICE inline fptype get_ang() const  { return atan2(im, re); }

// set functions
    // Cmpx& set_re(const fptype x);
    // Cmpx& set_im(const fptype x);

// complex number operations
    HOSTDEVICE inline Cmpx conjugate() const { return Cmpx(re, -im); }

// overloaded operators
    HOSTDEVICE void operator+=(const Cmpx &z);
    HOSTDEVICE void operator-=(const Cmpx &z);
    HOSTDEVICE void operator*=(const Cmpx &z);
    HOSTDEVICE void operator/=(const Cmpx &z);
    HOSTDEVICE void operator*=(const fptype k);

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
    HOST char* polar() const;
    HOST char* cartesian() const;
};

// inline Cmpx operator*(const fptype k, const Cmpx z) {
    // return z * k;
// }

#ifdef __CUDACC__
#include "Cmpx.cu"
#endif

#endif // #ifndef  _CMPX_H_
