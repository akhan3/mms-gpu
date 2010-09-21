#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cmath>
#include "vector_functions.hpp"


void divergence_3d( const Vector3 *V,
                    const int xdim, const int ydim, const int zdim,
                    fptype *S )
{
    // printf("calculating divergence... \n");
    fptype meshwidth = 1.0;
    for(int z = 0; z < zdim; z++)
        #ifdef _OPENMP
        // #pragma omp parallel for
        #endif
        for(int y = 0; y < ydim; y++)
            for(int x = 0; x < xdim; x++) {
                fptype x2 = (x != xdim-1) ? V[(z  )*ydim*xdim + (y  )*xdim + (x+1)].x : 0;
                fptype x0 = (x != 0     ) ? V[(z  )*ydim*xdim + (y  )*xdim + (x-1)].x : 0;
                fptype y2 = (y != ydim-1) ? V[(z  )*ydim*xdim + (y+1)*xdim + (x  )].y : 0;
                fptype y0 = (y != 0     ) ? V[(z  )*ydim*xdim + (y-1)*xdim + (x  )].y : 0;
                fptype z2 = (z != zdim-1) ? V[(z+1)*ydim*xdim + (y  )*xdim + (x  )].z : 0;
                fptype z0 = (z != 0     ) ? V[(z-1)*ydim*xdim + (y  )*xdim + (x  )].z : 0;
                S[z*ydim*xdim + y*xdim + x] = (x2 - x0 + y2 - y0 + z2 - z0) / (2*meshwidth);
            }
}


void gradient_3d(   const fptype *S,
                    const int xdim, const int ydim, const int zdim,
                    const fptype constant_multiple,
                    Vector3 *V )
{
    // printf("calculating gradient... \n");
    fptype meshwidth = 1.0;
    for(int z = 0; z < zdim; z++)
        #ifdef _OPENMP
        // #pragma omp parallel for
        #endif
        for(int y = 0; y < ydim; y++)
            for(int x = 0; x < xdim; x++) {
                fptype x2 = (x != xdim-1) ? S[(z  )*ydim*xdim + (y  )*xdim + (x+1)] : 0;
                fptype x0 = (x != 0     ) ? S[(z  )*ydim*xdim + (y  )*xdim + (x-1)] : 0;
                fptype y2 = (y != ydim-1) ? S[(z  )*ydim*xdim + (y+1)*xdim + (x  )] : 0;
                fptype y0 = (y != 0     ) ? S[(z  )*ydim*xdim + (y-1)*xdim + (x  )] : 0;
                fptype z2 = (z != zdim-1) ? S[(z+1)*ydim*xdim + (y  )*xdim + (x  )] : 0;
                fptype z0 = (z != 0     ) ? S[(z-1)*ydim*xdim + (y  )*xdim + (x  )] : 0;
                V[z*ydim*xdim + y*xdim + x].x = (x2 - x0) / (2*meshwidth) * constant_multiple;
                V[z*ydim*xdim + y*xdim + x].y = (y2 - y0) / (2*meshwidth) * constant_multiple;
                V[z*ydim*xdim + y*xdim + x].z = (z2 - z0) / (2*meshwidth) * constant_multiple;
            }
}


void add_exchange_field(const Vector3 *M,
                        const fptype Ms, const fptype Aexch, const fptype mu,
                        const int xdim, const int ydim, const int zdim, fptype meshwidth,
                        Vector3 *H )
{
    // printf("calculating exchange field... \n");
    const fptype constant_multiple = 2 * Aexch / (mu * Ms*Ms) / (meshwidth*meshwidth);
    for(int z = 0; z < zdim; z++)
        #ifdef _OPENMP
        // #pragma omp parallel for
        #endif
        for(int y = 0; y < ydim; y++)
            for(int x = 0; x < xdim; x++) {
                Vector3 x2 = (x != xdim-1) ? M[(z  )*ydim*xdim + (y  )*xdim + (x+1)] : Vector3(0,0,0);
                Vector3 x0 = (x != 0     ) ? M[(z  )*ydim*xdim + (y  )*xdim + (x-1)] : Vector3(0,0,0);
                Vector3 y2 = (y != ydim-1) ? M[(z  )*ydim*xdim + (y+1)*xdim + (x  )] : Vector3(0,0,0);
                Vector3 y0 = (y != 0     ) ? M[(z  )*ydim*xdim + (y-1)*xdim + (x  )] : Vector3(0,0,0);
                Vector3 z2 = (z != zdim-1) ? M[(z+1)*ydim*xdim + (y  )*xdim + (x  )] : Vector3(0,0,0);
                Vector3 z0 = (z != 0     ) ? M[(z-1)*ydim*xdim + (y  )*xdim + (x  )] : Vector3(0,0,0);
                H[z*ydim*xdim + y*xdim + x] += constant_multiple * (x2 + x0 + y2 + y0 + z2 + z0);
                // H[z*ydim*xdim + y*xdim + x] = H[z*ydim*xdim + y*xdim + x] + constant_multiple * (x2 + x0 + y2 + y0 + z2 + z0);
                // printf("Hexch_normalized = "); (1/(constant_multiple*Ms) * H[z*ydim*xdim + y*xdim + x]).print(stdout);
            }
}
