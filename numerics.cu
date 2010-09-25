#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cmath>
#include "Vector3.hpp"
#include "numerics.hpp"

HOSTDEVICE
int factorial(int x)
{
// Loop implementation
    int fac = 1;
    for(int i = 2; i <= x; i++)
        fac *= i;
    return fac;

// Recursive implementation
    // return (x == 0) ? 1 : x*factorial(x-1);

// LUT implementation
    // assert(x >= 0 && x <= 8);
    // switch (x) {
        // case 0: return     1;
        // case 1: return     1;
        // case 2: return     2;
        // case 3: return     6;
        // case 4: return    24;
        // case 5: return   120;
        // case 6: return   720;
        // case 7: return  5040;
        // case 8: return 40320;
        // default:
            // return (x == 0) ? 1 : x*factorial(x-1);
    // }
}

// Legendre function (recursive implementation)
HOSTDEVICE
fptype legendre(int k, fptype x)
{
    // assert(k >= 0);
    // assert(x >= -1 && x <= 1);
    switch (k) {
        case 0:
            return 1;
        case 1:
            return x;
        case 2:
            return 1/2.0 * (3 * x*x - 1);
        case 3:
            return 1/2.0 * (5 * x*x*x - 3 * x);
        default:
            return ((2*k-1) * x * legendre(k-1, x) - (k-1) * legendre(k-2, x)) / k;
    }
}

// Associated Legendre function
HOSTDEVICE
fptype associated_legendre(int l, int m, fptype x)
{
    // assert(l >= 0 && l <= 4);
    // assert(abs(m) <= l);
    // assert(m == abs(m));
    switch (l) {
        case 0:
            switch (m) {
                case 0:  return 1;
                default: return 0;
            }
        case 1:
            switch (m) {
                case 0:  return x;
                case 1:  return -sqrt(1 - x*x);
                default: return 0;
            }
        case 2:
            switch (m) {
                case 0:  return 1/2.0 * (3*x*x - 1);
                case 1:  return -3 * x * sqrt(1 - x*x);
                case 2:  return 3 * (1 - x*x);
                default: return 0;
            }
        case 3:
            switch (m) {
                case 0:  return 1/2.0 * x * (5*x*x - 3);
                case 1:  return -3/2.0 * (5*x*x - 1) * sqrt(1 - x*x);
                case 2:  return 15 * x * (1 - x*x);
                case 3:  return -15 * powf(1 - x*x, 1.5);
                default: return 0;
            }
        case 4:
            switch (m) {
                case 0:  return 1/8.0 * (x*x * (35*x*x - 30) + 3);
                case 1:  return -5/2.0 * x * (7*x*x - 3) * sqrt(1 - x*x);
                case 2:  return 15/2.0 * (7*x*x - 1) * (1 - x*x);
                case 3:  return -105 * x * powf(1 - x*x, 1.5);
                case 4:  return 105 * (1 - x*x) * (1 - x*x);
                default: return 0;
            }
        default:
            // printf("FATAL ERROR؛ associated_legendre(l=%d, m=%d, x=%f) is not implemented\n", l,m,x);
            return 0.0/0.0; // NaN
    }
}

// Spherical harmonics
HOSTDEVICE
Cmpx spherical_harmonic(int l, int m, fptype theta, fptype phi)
{
    // return Cmpx(sqrt((1.0*factorial(l-abs(m))) / factorial(l+abs(m))) * associated_legendre(l,abs(m),cos(theta)), m*phi, 1);
    fptype mag = associated_legendre(l, abs(m), cos(theta));
    fptype ang = m*phi;
    return Cmpx(mag*cos(ang), mag*sin(ang));
}


// random  number functions
HOST
int rand_atob(const int a, const int b) {
    double r = rand() / (double)RAND_MAX;
    r = a + (b-a+1) * r;
    return (int)r;
}

HOST
fptype frand_atob(const fptype a, const fptype b) {
    double r = rand() / (double)RAND_MAX;
    r = a + (b-a) * r;
    return (fptype)r;
}
