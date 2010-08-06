#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include "vector_functions.hpp"
#define NEWLINE printf("\n");


void divergence_3d( const Vector3 *V,
                    const int xdim, const int ydim, const int zdim,
                    float *S )
{
    // printf("calculating divergence... \n");
    float meshwidth = 1.0;
    for(int z = 0; z < zdim; z++)
        for(int y = 0; y < ydim; y++)
            for(int x = 0; x < xdim; x++) {
                float x2 = (x != xdim-1) ? V[(z  )*ydim*xdim + (y  )*xdim + (x+1)].x : 0;
                float x0 = (x != 0     ) ? V[(z  )*ydim*xdim + (y  )*xdim + (x-1)].x : 0;
                float y2 = (y != ydim-1) ? V[(z  )*ydim*xdim + (y+1)*xdim + (x  )].y : 0;
                float y0 = (y != 0     ) ? V[(z  )*ydim*xdim + (y-1)*xdim + (x  )].y : 0;
                float z2 = (z != zdim-1) ? V[(z+1)*ydim*xdim + (y  )*xdim + (x  )].z : 0;
                float z0 = (z != 0     ) ? V[(z-1)*ydim*xdim + (y  )*xdim + (x  )].z : 0;
                S[z*ydim*xdim + y*xdim + x] = (x2 - x0 + y2 - y0 + z2 - z0) / (2*meshwidth);
            }
}


void gradient_3d(   const float *S,
                    const int xdim, const int ydim, const int zdim,
                    Vector3 *V )
{
    // printf("calculating gradient... \n");
    float meshwidth = 1.0;
    for(int z = 0; z < zdim; z++)
        for(int y = 0; y < ydim; y++)
            for(int x = 0; x < xdim; x++) {
                float x2 = (x != xdim-1) ? S[(z  )*ydim*xdim + (y  )*xdim + (x+1)] : 0;
                float x0 = (x != 0     ) ? S[(z  )*ydim*xdim + (y  )*xdim + (x-1)] : 0;
                float y2 = (y != ydim-1) ? S[(z  )*ydim*xdim + (y+1)*xdim + (x  )] : 0;
                float y0 = (y != 0     ) ? S[(z  )*ydim*xdim + (y-1)*xdim + (x  )] : 0;
                float z2 = (z != zdim-1) ? S[(z+1)*ydim*xdim + (y  )*xdim + (x  )] : 0;
                float z0 = (z != 0     ) ? S[(z-1)*ydim*xdim + (y  )*xdim + (x  )] : 0;
                V[z*ydim*xdim + y*xdim + x].x = (x2 - x0) / (2*meshwidth);
                V[z*ydim*xdim + y*xdim + x].y = (y2 - y0) / (2*meshwidth);
                V[z*ydim*xdim + y*xdim + x].z = (z2 - z0) / (2*meshwidth);
            }
}


void add_exchange_field(const Vector3 *M,
                        const float Ms, const float Aexch, const float mu,
                        const int xdim, const int ydim, const int zdim, float meshwidth,
                        Vector3 *H )
{
    // printf("calculating exchange field... \n");
    const float constant_multiple = 2 * Aexch / (mu * Ms*Ms) / (meshwidth*meshwidth);
    float x2, x0, y2, y0, z2, z0;
    for(int z = 0; z < zdim; z++)
        for(int y = 0; y < ydim; y++)
            for(int x = 0; x < xdim; x++) {
                x2 = (x != xdim-1) ? M[(z  )*ydim*xdim + (y  )*xdim + (x+1)].x : 0;
                x0 = (x != 0     ) ? M[(z  )*ydim*xdim + (y  )*xdim + (x-1)].x : 0;
                y2 = (y != ydim-1) ? M[(z  )*ydim*xdim + (y+1)*xdim + (x  )].x : 0;
                y0 = (y != 0     ) ? M[(z  )*ydim*xdim + (y-1)*xdim + (x  )].x : 0;
                z2 = (z != zdim-1) ? M[(z+1)*ydim*xdim + (y  )*xdim + (x  )].x : 0;
                z0 = (z != 0     ) ? M[(z-1)*ydim*xdim + (y  )*xdim + (x  )].x : 0;
                H[z*ydim*xdim + y*xdim + x].x += constant_multiple * (x2 + x0 + y2 + y0 + z2 + z0);

                x2 = (x != xdim-1) ? M[(z  )*ydim*xdim + (y  )*xdim + (x+1)].y : 0;
                x0 = (x != 0     ) ? M[(z  )*ydim*xdim + (y  )*xdim + (x-1)].y : 0;
                y2 = (y != ydim-1) ? M[(z  )*ydim*xdim + (y+1)*xdim + (x  )].y : 0;
                y0 = (y != 0     ) ? M[(z  )*ydim*xdim + (y-1)*xdim + (x  )].y : 0;
                z2 = (z != zdim-1) ? M[(z+1)*ydim*xdim + (y  )*xdim + (x  )].y : 0;
                z0 = (z != 0     ) ? M[(z-1)*ydim*xdim + (y  )*xdim + (x  )].y : 0;
                H[z*ydim*xdim + y*xdim + x].y += constant_multiple * (x2 + x0 + y2 + y0 + z2 + z0);

                x2 = (x != xdim-1) ? M[(z  )*ydim*xdim + (y  )*xdim + (x+1)].x : 0;
                x0 = (x != 0     ) ? M[(z  )*ydim*xdim + (y  )*xdim + (x-1)].z : 0;
                y2 = (y != ydim-1) ? M[(z  )*ydim*xdim + (y+1)*xdim + (x  )].z : 0;
                y0 = (y != 0     ) ? M[(z  )*ydim*xdim + (y-1)*xdim + (x  )].z : 0;
                z2 = (z != zdim-1) ? M[(z+1)*ydim*xdim + (y  )*xdim + (x  )].z : 0;
                z0 = (z != 0     ) ? M[(z-1)*ydim*xdim + (y  )*xdim + (x  )].z : 0;
                H[z*ydim*xdim + y*xdim + x].z += constant_multiple * (x2 + x0 + y2 + y0 + z2 + z0);
            }
}
