#include <cmath>
#include <assert.h>
#include "Cmpx.hpp"

// constructor
Cmpx::Cmpx() {
    re = 0;
    im = 0;
}

Cmpx::Cmpx(const fptype x, const fptype y) {
    re = x;
    im = y;
}

// Cmpx& Cmpx::set_re(const fptype x) {
    // re = x;
    // return *this;
// }

// Cmpx& Cmpx::set_im(const fptype x) {
    // im = x;
    // return *this;
// }

// Cmpx Cmpx::conjugate() const {
    // return Cmpx(re, -im);
// }

// Cmpx& Cmpx::operator+=(const Cmpx &z) {
void Cmpx::operator+=(const Cmpx &z) {
    re += z.re;
    im += z.im;
    // return *this;
}

// Cmpx& Cmpx::operator-=(const Cmpx &z) {
    // re -= z.re;
    // im -= z.im;
    // return *this;
// }

// Cmpx& Cmpx::operator*=(const Cmpx &z) {
void Cmpx::operator*=(const Cmpx &z) {
    fptype x = re*z.re - im*z.im;
    fptype y = re*z.im + im*z.re;
    re = x;
    im = y;
    // return *this;
}

// Cmpx& Cmpx::operator*=(const fptype k) {
void Cmpx::operator*=(const fptype k) {
    re *= k;
    im *= k;
    // return *this;
}

char* Cmpx::polar() const {
    char *ch = new char[100];
    sprintf(ch, "%g<%gr", get_mag(), get_ang()); return ch;
}

char* Cmpx::cartesian() const {
    char *ch = new char[100];
    sprintf(ch, "%g+%gj", re, im); return ch;
}
