#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <assert.h>
#include <iostream>
#include <cmath>
#include <FreeImage.h>
#include "Box.hpp"
#include "Cmpx.hpp"
#include "Vector3.hpp"
#include "helper_functions.hpp"
#include "vector_functions.hpp"
#include "ode_functions.hpp"
#define NEWLINE printf("\n");
using std::cout;
using std::endl;


//*************************************************************************//
//******************** Main function **************************************//
//*************************************************************************//
int main(int argc, char **argv)
{
    const int verbose_level = 2;
    int status = 0;
    char filename_arg[1000];
    unsigned int seed = time(NULL);
    unsigned int P = 3;
    assert(argc >= 2);
    if(argc >= 2) {
        sscanf(argv[1], "%s", filename_arg);
    }
    if(argc >= 3) {
        sscanf(argv[2], "%u", &P);
        assert(P <= 4);
    }
    if(argc >= 4) {
        sscanf(argv[3], "%u", &seed);
    }
    srand(seed);


// Material parameters
// ================================================
    const float mu_0 = 4 * M_PI * 1e-7; // permeability of vacuum
    const float Ms = 8.6e5;             // saturation magnetization (permalloy)
    const float Aexch = 1.3e-11;        // exchange constant (permalloy)
    const float alfa = 0.5;             // damping coefficient (permalloy)
    const float gamma = 2.21e5;         // gyromagnetic ratio (permalloy)

// Mask configuration for magnetic material
    BYTE *mask = NULL; // mask matrix
    unsigned int xdim = 0;
    unsigned int ydim = 0;

    char filename[1000];
    sprintf(filename, "%s", filename_arg);
// read the mask from file
    load_mask(filename, &mask, &xdim, &ydim);
    assert(xdim == ydim);
    unsigned int zdim = 1;
    printf("(xdim, ydim, zdim) = (%d, %d, %d)\n", xdim, ydim, zdim);

// generate the initial magnetization distribution
    Vector3 *M = new Vector3[zdim*ydim*xdim]();  // magnetization matrix
    if(M == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

    for(unsigned int z = 0; z < zdim; z++)
    // unsigned int z = 0;
    {
        for(unsigned int y = 0; y < ydim; y++) {
            for(unsigned int x = 0; x < xdim; x++) {
                if (!mask[y*xdim + x])
                {
                    // float theta = frand_atob(90-30, 90+30) * M_PI/180;
                    float theta = frand_atob(0, 180) * M_PI/180;
                    float phi   = frand_atob(0, 360) * M_PI/180;
                    M[z*ydim*xdim + y*xdim + x] = Ms * Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
                }
            }
        }
    }
    delete []mask;

// write M field to file
    status |= save_vector3d(M, zdim, ydim, xdim, "M", verbose_level);
    if(status) return EXIT_FAILURE;

// set up grid and simulation time
    const float meshwidth = 1e-9;
    const float dt = 1e-14;           // initial time step = 1ps (will be adapted on the fly)
    const float finaltime = 1e-10;

// magnetization dynamics
// ===================================================================
    printf("If you want to see the initial state of M, now is the time! "); fflush(NULL);
    getchar();

    status |= time_marching(    M,
                                dt, finaltime,
                                xdim, ydim, zdim, meshwidth,
                                mu_0, Ms, Aexch, alfa, gamma,
                                verbose_level );
    if(status) return EXIT_FAILURE;

// write M vectorfield to file
    status |= save_vector3d(M, zdim, ydim, xdim, "M", verbose_level);
    if(status) return EXIT_FAILURE;

// closing
    delete []M;

    printf("SEED = %d\n", seed);
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
