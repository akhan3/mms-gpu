#include <cmath>
#include <assert.h>
#include "Cmpx.hpp"

// constructor
Cmpx::Cmpx() {
    re = 0;
    im = 0;
    mag = 0;
    ang = 0;
}

Cmpx::Cmpx(const float x, const float y, const int mode) {
    init(x, y, mode);
}

Cmpx& Cmpx::init(const float x, const float y, const int mode) {
    if (mode == 0) {
        re = x;
        im = y;
        mag = magnitude();
        ang = angle();
    }
    else {
        // assert(x >= 0);
        if(x < 0) {
            mag = -x;
            ang = fmod(y + M_PI, 2*M_PI);
        }
        else {
            mag = x;
            ang = fmod(y, 2*M_PI);
        }
        re = real();
        im = imaginary();
    }
    return *this;
}

Cmpx& Cmpx::set_re(const float x) {
    re = x;
    mag = magnitude();
    ang = angle();
    return *this;
}

Cmpx& Cmpx::set_im(const float x) {
    im = x;
    mag = magnitude();
    ang = angle();
    return *this;
}

Cmpx& Cmpx::set_mag(const float x) {
    assert(x >= 0);
    mag = x;
    re = real();
    im = imaginary();
    return *this;
}

Cmpx& Cmpx::set_ang(const float x) {
    ang = fmod(x, 2*M_PI);
    re = real();
    im = imaginary();
    return *this;
}

Cmpx Cmpx::conjugate() const {
    return Cmpx(re, -im, 0);
}

Cmpx& Cmpx::operator+=(const Cmpx &z) {
    re += z.re;
    im += z.im;
    mag = magnitude();
    ang = angle();
    return *this;
}

Cmpx& Cmpx::operator-=(const Cmpx &z) {
    re -= z.re;
    im -= z.im;
    mag = magnitude();
    ang = angle();
    return *this;
}

Cmpx& Cmpx::operator*=(const Cmpx &z) {
    mag *= z.mag;
    ang = fmod(ang + z.ang, 2*M_PI);
    re = real();
    im = imaginary();
    return *this;
}

Cmpx& Cmpx::operator/=(const Cmpx &z) {
    mag /= z.mag;
    ang = fmod(ang - z.ang, 2*M_PI);
    re = real();
    im = imaginary();
    return *this;
}

Cmpx& Cmpx::operator*=(const float k) {
    if(k < 0) {
        mag *= -k;
        ang = fmod(ang + M_PI, 2*M_PI);
    }
    else {
        mag *= k;
    }
    re = real();
    im = imaginary();
    return *this;
}
