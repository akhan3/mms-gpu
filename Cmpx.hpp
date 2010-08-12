#ifndef _CMPX_H_
#define _CMPX_H_

#include <stdio.h>
#include "mydefs.hpp"

class Cmpx {
// data members
private:
    fptype re;   // z = re + j*im
    fptype im;
    fptype mag;  // z = mag * exp(j*ang)
    fptype ang;

public:
// constructors
    Cmpx();
    Cmpx(const fptype x, const fptype y, const int mode); // mode=0 for cartesian

// member functions
// get functions
    inline fptype get_re()  const  { return re;  }
    inline fptype get_im()  const  { return im;  }
    inline fptype get_mag() const  { return mag; }
    inline fptype get_ang() const  { return ang; }

// set functions
    Cmpx& init(const fptype x, const fptype y, const int mode); // mode=0 for cartesian
    Cmpx& set_re(const fptype x);
    Cmpx& set_im(const fptype x);
    Cmpx& set_mag(const fptype x);
    Cmpx& set_ang(const fptype x);

// complex number operations
    Cmpx  conjugate() const;

// overloaded operators
    Cmpx& operator+=(const Cmpx &z);
    Cmpx& operator-=(const Cmpx &z);
    Cmpx& operator*=(const Cmpx &z);
    Cmpx& operator/=(const Cmpx &z);
    Cmpx& operator*=(const fptype k);
    // inline Cmpx operator+(const Cmpx &z) const  { Cmpx z1 = *this; return z1 += z; }
    inline Cmpx operator+(const Cmpx &z) const  { return Cmpx(*this) += z; }
    inline Cmpx operator-(const Cmpx &z) const  { return Cmpx(*this) -= z; }
    inline Cmpx operator*(const Cmpx &z) const  { return Cmpx(*this) *= z; }
    inline Cmpx operator/(const Cmpx &z) const  { return Cmpx(*this) /= z; }
    inline Cmpx operator*(const fptype k) const  { return Cmpx(*this) *= k; }

// print methods
    inline char* polar() const  { char *ch = new char[100]; sprintf(ch, "%g<%gr", mag, ang); return ch; }
    inline char* cartesian() const  { char *ch = new char[100]; sprintf(ch, "%g+%gj", re, im); return ch; }

private:
    inline fptype magnitude() const  { return sqrt(re*re + im*im); }
    inline fptype angle() const      { return atan2(im, re); }
    inline fptype real() const       { return mag * cos(ang); }
    inline fptype imaginary() const  { return mag * sin(ang); }
};

inline Cmpx operator*(const fptype k, const Cmpx z) {
    return z * k;
}

#endif // #ifndef  _CMPX_H_
