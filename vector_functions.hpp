#ifndef _VECTOR_FUNCTIONS_H_
#define _VECTOR_FUNCTIONS_H_

#include "mydefs.hpp"
#include "Cmpx.hpp"
#include "Vector3.hpp"


void divergence_3d( const Vector3 *V,
                    const int xdim, const int ydim, const int zdim,
                    const fptype dx, const fptype dy, const fptype dz,
                    fptype *S );

void gradient_3d(   const fptype *S,
                    const int xdim, const int ydim, const int zdim,
                    const fptype dx, const fptype dy, const fptype dz,
                    const fptype constant_multiple,
                    Vector3 *V );

void add_exchange_field(
                    const Vector3 *M,
                    const fptype Ms, const fptype Aexch, const fptype mu,
                    const int xdim, const int ydim, const int zdim,
                    const fptype dx, const fptype dy, const fptype dz,
                    Vector3 *H );

#endif // #ifndef  _VECTOR_FUNCTIONS_H_
