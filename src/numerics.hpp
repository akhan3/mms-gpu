#ifndef _NUMERICS_H_
#define _NUMERICS_H_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cmath>
// #include "Vector3.hpp"
#include "Cmpx.hpp"

HOSTDEVICE  int factorial(int x);
HOSTDEVICE  fptype legendre(int k, fptype x);
HOSTDEVICE  fptype associated_legendre(int l, int m, fptype x);
HOSTDEVICE  Cmpx spherical_harmonic(int l, int m, fptype theta, fptype phi);
HOST        int rand_atob(const int a, const int b);
HOST        fptype frand_atob(const fptype a, const fptype b);

#ifdef __CUDACC__
#include "numerics.cu"
#endif

#endif // #ifndef  _NUMERICS_H_
