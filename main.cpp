#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <assert.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include "Box.hpp"
#include "Cmpx.hpp"
#include "Vector3.hpp"
#include "numerics.hpp"
#include "helper_functions.hpp"
#include "vector_functions.hpp"
#include "ode_functions.hpp"

//*************************************************************************//
//******************** Main function **************************************//
//*************************************************************************//
int main(int argc, char **argv)
{
    int status = 0;
    const int verbosity = 4;

    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);

// intial command line arguments
    char filename_arg[1000] = "verysmall_16x16.png";
    unsigned int P = 3;
    fptype finaltime = 1e-9;
    fptype timestep = 1e-14;
    fptype sample_width  = 10e-9;
    fptype sample_height = 10e-9;
    fptype sample_depth  = 1e-9;
    int xdim = 16;
    int ydim = 16;
    int zdim = 3;
    int demag = true;
    int exchange = true;
    int external = false;
    int use_fmm = false;
    int use_gpu = false;
    char sim_name[1000] = "sim_untitled";
    unsigned int seed = time(NULL);
    unsigned int IC = 0; // 0 = SD, 1 = Vortex, 2 = Random
// read command line arguments
    if(argc >= 2)
        sscanf(argv[1], "%s", filename_arg);
    if(argc >= 3) {
        sscanf(argv[2], "%u", &P);
        assert(P <= 4); // very important
    }
    if(argc >=  4)  sscanf(argv[ 3], "%f", &finaltime);
    if(argc >=  5)  sscanf(argv[ 4], "%f", &timestep);
    if(argc >=  6)  sscanf(argv[ 5], "%f", &sample_width);
    if(argc >=  7)  sscanf(argv[ 6], "%f", &sample_height);
    if(argc >=  8)  sscanf(argv[ 7], "%f", &sample_depth);
    if(argc >=  9)  sscanf(argv[ 8], "%d", &xdim);
    if(argc >= 10)  sscanf(argv[ 9], "%d", &ydim);
    if(argc >= 11)  sscanf(argv[10], "%d", &zdim);
    if(argc >= 12)  sscanf(argv[11], "%d", &demag);
    if(argc >= 13)  sscanf(argv[12], "%d", &exchange);
    if(argc >= 14)  sscanf(argv[13], "%d", &external);
    if(argc >= 15)  sscanf(argv[14], "%d", &use_fmm);
    if(argc >= 16)  sscanf(argv[15], "%d", &use_gpu);
    if(argc >= 17)  sscanf(argv[16], "%s", sim_name);
    if(argc >= 18)  sscanf(argv[17], "%u", &seed);
    if(argc >= 19)  sscanf(argv[18], "%u", &IC);
    srand(seed);
// print command line arguments
#ifdef _OPENMP
    printf("Compiled with OpenMP and running with %s threads.\n", getenv("OMP_NUM_THREADS"));
#endif
    if(verbosity >= 2) {
        printf("imagefile = %s \n", filename_arg);
        printf("P = %d \n", P);
        printf("finaltime = %g \n", finaltime);
        printf("timestep = %g \n", timestep);
        printf("sample_width  = %g \n", sample_width);
        printf("sample_height = %g \n", sample_height);
        printf("sample_depth  = %g \n", sample_depth);
        printf("xdim = %d \n", xdim);
        printf("ydim = %d \n", ydim);
        printf("zdim = %d \n", zdim);
        printf("demag = %d \n", demag);
        printf("exchange = %d \n", exchange);
        printf("external = %d \n", external);
        printf("use_fmm = %d \n", use_fmm);
        printf("use_gpu = %d \n", use_gpu);
        printf("sim_name = %s \n", sim_name);
        printf("SEED = %d \n", seed);
        printf("IC = %d \n", IC);
    }

// create directory to hold results
    // char dirname[1000];
    // sprintf(dirname, "%s/%s", "results_square", sim_name);
    // mkdir(dirname, 0755);
    mkdir(sim_name, 0755);
    char filename[1000];
    // if(status) return EXIT_FAILURE;

// Material parameters
// ================================================
    const fptype mu_0 = 4 * M_PI * 1e-7; // permeability of vacuum
    const fptype Ms = 8.6e5;             // saturation magnetization (permalloy)
    const fptype Aexch = 1.3e-11;        // exchange constant (permalloy)
    const fptype alfa = 0.5;             // damping coefficient (permalloy)
    const fptype gamma = 2.21e5;         // gyromagnetic ratio (permalloy)
    const fptype dx = sample_width  / xdim;
    const fptype dy = sample_height / ydim;
    const fptype dz = sample_depth  / zdim;
    assert(dx == dy);
    assert(dy == dz);
    assert(dx <= 2.5e-9);

// Mask configuration for magnetic material
#ifdef USE_FREEIMAGE
    BYTE *mask = NULL; // mask matrix
    sprintf(filename, "%s", filename_arg);
// read the mask from file
    load_mask(filename, &mask, &xdim, &ydim);
#else
    byte *mask = new byte[ydim*xdim](); // mask matrix
// specimen magnet 20x20x20
    for(int y = 0; y < ydim; y++)
        for(int x = 0; x < xdim; x++)
            mask[y*xdim + x] = 1;   // all white (no material)
    // for(unsigned int y = 21; y <= 41; y++)
        // for(unsigned int x = 21; x <= 41; x++)
    for(int y = 1; y < ydim-1; y++)
        for(int x = 1; x < xdim-1; x++)
            mask[y*xdim + x] = 0;   // selected black (material)
#endif
    // assert(xdim == ydim);
    printf("(xdim, ydim, zdim) = (%d, %d, %d)\n", xdim, ydim, zdim);
    printf("(sample_width, sample_height, sample_depth) = (%g, %g, %g)\n", sample_width, sample_height, sample_depth);
    printf("(dx, dy, dz) = (%g, %g, %g)\n", dx, dy, dz);

// determine initial condition
    int IC_singledomain = 0;
    int IC_vortex = 0;
    int IC_random = 0;
    switch(IC) {
        case 0: // single domain
            IC_singledomain = 1;
            printf("Initial condition: Single Domain\n");
            break;
        case 1: // vortex
            IC_vortex = 1;
            printf("Initial condition: Vortex\n");
            break;
        case 2: // random
            IC_random = 1;
            printf("Initial condition: Random\n");
            break;
        default:
            fprintf(stderr, "ERROR: Unknown Initial Condition!\n");
            return EXIT_FAILURE;
    }

// generate the initial magnetization distribution
    byte *material = new byte[zdim*ydim*xdim]();  // material matrix
    Vector3 *M = new Vector3[zdim*ydim*xdim]();  // magnetization matrix
    if(material == NULL || M == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

// if(zdim == 1)
    // for(unsigned int z = 0; z < zdim; z++) {
        // for(unsigned int y = 0; y < ydim; y++) {
            // for(unsigned int x = 0; x < xdim; x++) {
                // if(!mask[y*xdim + x])
                // {
                    // // fptype theta = frand_atob(0, 180) * M_PI/180;
                    // fptype phi   = frand_atob(0, 360) * M_PI/180;
                    // fptype theta = M_PI/2;
                    // // fptype phi   = 0;
                    // M[z*ydim*xdim + y*xdim + x] = Ms * Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
                    // material[z*ydim*xdim + y*xdim + x] = 1;
                // }
            // }
        // }
    // }
// else if(zdim >= 3)

    assert(zdim >= 3);
    for(int z = 1; z < zdim-1; z++) {
        for(int y = 0; y < ydim; y++) {
            for(int x = 0; x < xdim; x++) {
                if(!mask[y*xdim + x])
                {
                    fptype theta = M_PI/2;
                    fptype phi   = 0;
                    // fptype theta = frand_atob(0, 180) * M_PI/180;
                    // fptype phi   = frand_atob(0, 360) * M_PI/180;
                    // fptype theta = M_PI/2 + frand_atob(-10, 10) * M_PI/180;
                    // fptype phi   = 0 + frand_atob(-90, 90) * M_PI/180;
                    if(IC_singledomain)
                        phi = 0;
                    else if(IC_random)
                        phi   = frand_atob(0, 360) * M_PI/180;
                    else if(IC_vortex) {
                        if((x-xdim/2 > 0) && (y-ydim/2 > 0))
                            phi = -M_PI/4;
                        else if((x-xdim/2 <= 0) && (y-ydim/2 > 0))
                            phi = M_PI/4;
                        else if((x-xdim/2 <= 0) && (y-ydim/2 <= 0))
                            phi = 3*M_PI/4;
                        else if((x-xdim/2 > 0) && (y-ydim/2 <= 0))
                            phi = -3*M_PI/4;
                    }

                    M[z*ydim*xdim + y*xdim + x] = Ms * Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
                    material[z*ydim*xdim + y*xdim + x] = 1;
                }
            }
        }
        // fptype theta = 0;
        // fptype phi   = frand_atob(0, 360) * M_PI/180;
        // M[z*ydim*xdim + ydim/2*xdim + xdim/2] = Ms * Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
    }

    delete []mask;

// write M field to file
    sprintf(filename, "%s/%s", sim_name, "M.dat");
    status |= save_vector3d(M, zdim, ydim, xdim, filename, verbosity);
    if(status) return EXIT_FAILURE;

// write material field to file
    fptype *m = new fptype[xdim*ydim];
    int z = 1;
    for(int y = 0; y < ydim; y++)
        for(int x = 0; x < xdim; x++)
            m[y*xdim + x] = (fptype)material[z*ydim*xdim + y*xdim + x];
    sprintf(filename, "%s/%s", sim_name, "material.dat");
    status |= matrix2file(m, ydim, xdim, filename, verbosity);
    if(status) return EXIT_FAILURE;
    delete []m;

// magnetization dynamics
// ===================================================================
    status |= time_marching(    material, M,
                                finaltime, timestep,
                                xdim, ydim, zdim, dx, dy, dz,
                                P, mu_0, Ms, Aexch, alfa, gamma,
                                demag, exchange, external, use_fmm,
                                use_gpu, sim_name, verbosity );

// closing
    delete []M;
    delete []material;

    printf("SEED = %d\n", seed);
    // printf("%s\n", status ? "failed to complete" : "successfuly completed");

    status |= gettimeofday(&time2, NULL);
    double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    printf("Simulation completed in %f seconds.\n", deltatime);

    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
