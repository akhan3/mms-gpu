#ifndef _CMPX_H_
#define _CMPX_H_

#include <stdio.h>

class Cmpx {
// data members
public:
    float x;   // z = x + j*y
    float y;

public:
// constructors
    Cmpx();
    Cmpx(const float x, const float y);

// member functions
    Cmpx& init(const float x, const float y);
                                                                                        // // get functions
                                                                                            // inline float get_re()  const  { return re;  }
                                                                                            // inline float get_im()  const  { return im;  }
                                                                                            // inline float get_mag() const  { return mag; }
                                                                                            // inline float get_ang() const  { return ang; }

                                                                                        // // set functions
                                                                                            // Cmpx& set_re(const float x);
                                                                                            // Cmpx& set_im(const float x);
                                                                                            // Cmpx& set_mag(const float x);
                                                                                            // Cmpx& set_ang(const float x);

// complex number operations
    Cmpx  conjugate() const;
    inline float magnitude() const  { return sqrt(x*x + y*y); }
    inline float angle() const      { return atan2(y, x); }

// overloaded operators
    Cmpx& operator+=(const Cmpx &z);
    Cmpx& operator-=(const Cmpx &z);
    Cmpx& operator*=(const Cmpx &z);
    // Cmpx& operator/=(const Cmpx &z);
    Cmpx& operator*=(const float k);
    // inline Cmpx operator+(const Cmpx &z) const  { Cmpx z1 = *this; return z1 += z; }
    inline Cmpx operator+(const Cmpx &z) const  { return Cmpx(*this) += z; }
    inline Cmpx operator-(const Cmpx &z) const  { return Cmpx(*this) -= z; }
    inline Cmpx operator*(const Cmpx &z) const  { return Cmpx(*this) *= z; }
    // inline Cmpx operator/(const Cmpx &z) const  { return Cmpx(*this) /= z; }
    inline Cmpx operator*(const float k) const  { return Cmpx(*this) *= k; }

// print methods
    inline char* cartesian() const  { char *ch = new char[100]; sprintf(ch, "%g+%gj", x, y); return ch; }
    inline char* polar() const  { char *ch = new char[100]; sprintf(ch, "%g<%gr", magnitude(), angle()); return ch; }

};

inline Cmpx operator*(const float k, const Cmpx z) {
    return z * k;
}

#endif // #ifndef  _CMPX_H_
