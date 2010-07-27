#ifndef _CMPX_H_
#define _CMPX_H_

#include <stdio.h>

class Cmpx {
// data members
public:
    float re;   // z = re + j*im
    float im;
    float mag;  // z = mag * exp(j*ang)
    float ang;

// constructor
    Cmpx();
    Cmpx(const float x, const float y, const int mode); // mode=0 for cartesian

// member functions
private:
    inline float magnitude() const  { return sqrt(re*re + im*im); }
    inline float angle() const      { return atan2(im, re); }
    inline float real() const       { return mag * cos(ang); }
    inline float imaginary() const  { return mag * sin(ang); }

public:
    Cmpx conjugate() const;
    Cmpx operator+=(const Cmpx &z);
    Cmpx operator-=(const Cmpx &z);
    Cmpx operator*=(const Cmpx &z);
    Cmpx operator/=(const Cmpx &z);
    Cmpx operator*=(const float k);
    // inline Cmpx operator+(const Cmpx &z) const  { Cmpx z1 = *this; return z1 += z; }
    inline Cmpx operator+(const Cmpx &z) const  { return Cmpx(*this) += z; }
    inline Cmpx operator-(const Cmpx &z) const  { return Cmpx(*this) -= z; }
    inline Cmpx operator*(const Cmpx &z) const  { return Cmpx(*this) *= z; }
    inline Cmpx operator/(const Cmpx &z) const  { return Cmpx(*this) /= z; }
    inline Cmpx operator*(const float k) const  { return Cmpx(*this) *= k; }

    inline char* polar() const  { char *ch = new char[100]; sprintf(ch, "%10.6f<%4.0fd", mag, ang*180/M_PI); return ch; }
    inline char* cartesian() const  { char *ch = new char[100]; sprintf(ch, "%g+%gj", re, im); return ch; }
};

inline Cmpx operator*(const float k, const Cmpx z) {
    return z * k;
}

#endif // #ifndef  _CMPX_H_
