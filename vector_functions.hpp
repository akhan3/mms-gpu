#ifndef _VECTOR_FUNCTIONS_H_
#define _VECTOR_FUNCTIONS_H_

#include "Cmpx.hpp"
#include "Vector3.hpp"


void divergence_3d( const Vector3 *V,
                    const int xdim, const int ydim, const int zdim,
                    float *S );

void gradient_3d(   const float *S,
                    const int xdim, const int ydim, const int zdim,
                    Vector3 *V );

void add_exchange_field(
                    const Vector3 *M,
                    const float Ms, const float Aexch, const float mu,
                    const int xdim, const int ydim, const int zdim, const float meshwidth,
                    Vector3 *H );

#endif // #ifndef  _VECTOR_FUNCTIONS_H_
