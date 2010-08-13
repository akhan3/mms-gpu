#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <assert.h>
#include <iostream>
#include <cmath>
#include "Box.hpp"
#include "Cmpx.hpp"
#include "Vector3.hpp"
#include "helper_functions.hpp"
#include "vector_functions.hpp"
#include "ode_functions.hpp"
using std::cout;
using std::endl;


//*************************************************************************//
//******************** Main function **************************************//
//*************************************************************************//
int main(int argc, char **argv)
{
    const int verbose_level = 2;
    int status = 0;
// intial command line arguments
    char filename_arg[1000] = "verysmall_16x16.png";
    unsigned int P = 3;
    fptype finaltime = 1e-9;
    fptype timestep = 1e-14;
    fptype meshwidth = 1e-9;
    int coupling = true;
    int exchange = true;
    int external = false;
    int use_fmm = false;
    // unsigned int seed = time(NULL);
    unsigned int seed = 1985;
// read command line arguments
    if(argc >= 2)
        sscanf(argv[1], "%s", filename_arg);
    if(argc >= 3) {
        sscanf(argv[2], "%u", &P);
        assert(P <= 4);
    }
    if(argc >= 4)
        sscanf(argv[3], "%f", &finaltime);
    if(argc >= 5)
        sscanf(argv[4], "%f", &timestep);
    if(argc >= 6)
        sscanf(argv[5], "%f", &meshwidth);
    if(argc >= 7)
        sscanf(argv[6], "%d", &coupling);
    if(argc >= 8)
        sscanf(argv[7], "%d", &exchange);
    if(argc >= 9)
        sscanf(argv[8], "%d", &external);
    if(argc >= 10)
        sscanf(argv[9], "%d", &use_fmm);
    if(argc >= 11)
        sscanf(argv[10], "%u", &seed);
    srand(seed);

// Material parameters
// ================================================
    const fptype mu_0 = 4 * M_PI * 1e-7; // permeability of vacuum
    const fptype Ms = 8.6e5;             // saturation magnetization (permalloy)
    const fptype Aexch = 1.3e-11;        // exchange constant (permalloy)
    const fptype alfa = 0.5;             // damping coefficient (permalloy)
    const fptype gamma = 2.21e5;         // gyromagnetic ratio (permalloy)

// Mask configuration for magnetic material
    unsigned int xdim = 0;
    unsigned int ydim = 0;
    unsigned int zdim = 0;
#ifdef USE_FREEIMAGE
    BYTE *mask = NULL; // mask matrix
    char filename[1000];
    sprintf(filename, "%s", filename_arg);
// read the mask from file
    load_mask(filename, &mask, &xdim, &ydim);
    zdim = 1;
#else
    xdim = 64;
    ydim = xdim;
    zdim = 32;
    byte *mask = new byte[ydim*xdim](); // mask matrix
    for(unsigned int y = 1; y < ydim-1; y++)
        for(unsigned int x = 1; x < xdim-1; x++)
            mask[y*xdim + x] = 1;
#endif
    assert(xdim == ydim);
    printf("(xdim, ydim, zdim) = (%d, %d, %d)\n", xdim, ydim, zdim);

// generate the initial magnetization distribution
    Vector3 *M = new Vector3[zdim*ydim*xdim]();  // magnetization matrix
    if(M == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

    unsigned int zstart;
    unsigned int zend;
    if(zdim == 1) {
        zstart = 0;
        zend = 0;
    }
    else {
        zstart = 1;
        zend = zdim-2;
    }
    for(unsigned int z = zstart; z <= zend; z++) {
        for(unsigned int y = 0; y < ydim; y++) {
            for(unsigned int x = 0; x < xdim; x++) {
                if (!mask[y*xdim + x])
                {
                    fptype theta = frand_atob(90-30, 90+30) * M_PI/180;
                    fptype phi   = frand_atob(0, 90) * M_PI/180;
                    // fptype theta = frand_atob(0, 180) * M_PI/180;
                    // fptype phi   = frand_atob(0, 360) * M_PI/180;
                    // fptype theta = M_PI/2;
                    // fptype phi   = 0;
                    M[z*ydim*xdim + y*xdim + x] = Ms * Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
                }
            }
        }
    }
    delete []mask;

// write M field to file
    status |= save_vector3d(M, zdim, ydim, xdim, "M.dat", verbose_level);
    if(status) return EXIT_FAILURE;

// set up grid and simulation time
    // const fptype meshwidth = 1e-9;
    // const fptype finaltime = 10e-9;

// magnetization dynamics
// ===================================================================
    status |= time_marching(    M,
                                finaltime, timestep,
                                xdim, ydim, zdim, meshwidth, P,
                                mu_0, Ms, Aexch, alfa, gamma,
                                coupling, exchange, external, use_fmm,
                                verbose_level );
    if(status) return EXIT_FAILURE;

// // write M vectorfield to file
    // status |= save_vector3d(M, zdim, ydim, xdim, "M.dat", verbose_level);
    // if(status) return EXIT_FAILURE;

// closing
    delete []M;

    printf("SEED = %d\n", seed);
    // printf("%s\n", status ? "failed to complete" : "successfuly completed");
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
