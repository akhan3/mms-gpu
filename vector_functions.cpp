#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cmath>
#include "vector_functions.hpp"


void divergence_3d( const Vector3 *V,
                    const int xdim, const int ydim, const int zdim,
                    const fptype dx, const fptype dy, const fptype dz,
                    fptype *S )
{
    // printf("calculating divergence... \n");
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
                S[z*ydim*xdim + y*xdim + x] = ((x2-x0)/dx + (y2-y0)/dy + (z2-z0)/dz) / 2.0;
            }
}


void gradient_3d(   const fptype *S,
                    const int xdim, const int ydim, const int zdim,
                    const fptype dx, const fptype dy, const fptype dz,
                    const fptype constant_multiple,
                    Vector3 *V )
{
    // printf("calculating gradient... \n");
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
                V[z*ydim*xdim + y*xdim + x].x = (x2 - x0) / (2*dx) * constant_multiple;
                V[z*ydim*xdim + y*xdim + x].y = (y2 - y0) / (2*dy) * constant_multiple;
                V[z*ydim*xdim + y*xdim + x].z = (z2 - z0) / (2*dz) * constant_multiple;
            }
}


void add_exchange_field(const Vector3 *M,
                        const fptype Ms, const fptype Aexch, const fptype mu,
                        const int xdim, const int ydim, const int zdim,
                        const fptype dx, const fptype dy, const fptype dz,
                        Vector3 *H )
{
    // printf("calculating exchange field... \n");
    const fptype constant_multiple = 2 * Aexch / (mu * Ms*Ms);
    for(int z = 0; z < zdim; z++)
        #ifdef _OPENMP
        // #pragma omp parallel for
        #endif
        for(int y = 0; y < ydim; y++)
            for(int x = 0; x < xdim; x++) {
                fptype L_Mx = ( ((x != xdim-1) ? M[(z  )*ydim*xdim + (y  )*xdim + (x+1)].x : 0)
                               +((x != 0     ) ? M[(z  )*ydim*xdim + (y  )*xdim + (x-1)].x : 0)
                               -(            2 * M[(z  )*ydim*xdim + (y  )*xdim + (x  )].x    ) ) / (dx*dx)
                             +( ((y != ydim-1) ? M[(z  )*ydim*xdim + (y+1)*xdim + (x  )].x : 0)
                               +((y != 0     ) ? M[(z  )*ydim*xdim + (y-1)*xdim + (x  )].x : 0)
                               -(            2 * M[(z  )*ydim*xdim + (y  )*xdim + (x  )].x    ) ) / (dy*dy)
                             +( ((z != zdim-1) ? M[(z+1)*ydim*xdim + (y  )*xdim + (x  )].x : 0)
                               +((z != 0     ) ? M[(z-1)*ydim*xdim + (y  )*xdim + (x  )].x : 0)
                               -(            2 * M[(z  )*ydim*xdim + (y  )*xdim + (x  )].x    ) ) / (dz*dz);

                fptype L_My = ( ((x != xdim-1) ? M[(z  )*ydim*xdim + (y  )*xdim + (x+1)].y : 0)
                               +((x != 0     ) ? M[(z  )*ydim*xdim + (y  )*xdim + (x-1)].y : 0)
                               -(            2 * M[(z  )*ydim*xdim + (y  )*xdim + (x  )].y    ) ) / (dx*dx)
                             +( ((y != ydim-1) ? M[(z  )*ydim*xdim + (y+1)*xdim + (x  )].y : 0)
                               +((y != 0     ) ? M[(z  )*ydim*xdim + (y-1)*xdim + (x  )].y : 0)
                               -(            2 * M[(z  )*ydim*xdim + (y  )*xdim + (x  )].y    ) ) / (dy*dy)
                             +( ((z != zdim-1) ? M[(z+1)*ydim*xdim + (y  )*xdim + (x  )].y : 0)
                               +((z != 0     ) ? M[(z-1)*ydim*xdim + (y  )*xdim + (x  )].y : 0)
                               -(            2 * M[(z  )*ydim*xdim + (y  )*xdim + (x  )].y    ) ) / (dz*dz);

                fptype L_Mz = ( ((x != xdim-1) ? M[(z  )*ydim*xdim + (y  )*xdim + (x+1)].z : 0)
                               +((x != 0     ) ? M[(z  )*ydim*xdim + (y  )*xdim + (x-1)].z : 0)
                               -(            2 * M[(z  )*ydim*xdim + (y  )*xdim + (x  )].z    ) ) / (dx*dx)
                             +( ((y != ydim-1) ? M[(z  )*ydim*xdim + (y+1)*xdim + (x  )].z : 0)
                               +((y != 0     ) ? M[(z  )*ydim*xdim + (y-1)*xdim + (x  )].z : 0)
                               -(            2 * M[(z  )*ydim*xdim + (y  )*xdim + (x  )].z    ) ) / (dy*dy)
                             +( ((z != zdim-1) ? M[(z+1)*ydim*xdim + (y  )*xdim + (x  )].z : 0)
                               +((z != 0     ) ? M[(z-1)*ydim*xdim + (y  )*xdim + (x  )].z : 0)
                               -(            2 * M[(z  )*ydim*xdim + (y  )*xdim + (x  )].z    ) ) / (dz*dz);

                H[z*ydim*xdim + y*xdim + x] +=  constant_multiple * Vector3(L_Mx, L_My, L_Mz);
            }
}
