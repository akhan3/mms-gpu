#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <assert.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <errno.h>
#include "globals.hpp"
#include "Box.hpp"
#include "Cmpx.hpp"
#include "Vector3.hpp"
#include "numerics.hpp"
#include "helper_functions.hpp"
#include "vector_functions.hpp"
#include "ode_functions.hpp"

//*************************************************************************//
//******************** intial command line arguments **********************//
//*************************************************************************//
    // DEFINE_string   (maskImg, "", "Text format mask file for the structure of magnetic material");
    // DEFINE_string   (maskTxt, "", "PNG format mask file for the structure of magnetic material");
    // DEFINE_double   (sample_width,  100e-9, "width (X) of magnetic sample");
    // DEFINE_double   (sample_height, 100e-9, "height (Y) of magnetic sample");
    // DEFINE_double   (sample_depth,  5e-9, "depth (Z) of magnetic sample");
    DEFINE_int32    (Nx, 13, "number of mesh points in X-direction");
    DEFINE_int32    (Ny, 13, "number of mesh points in Y-direction");
    DEFINE_int32    (Nlayers, 3, "number of layers of magnetic materials");
    DEFINE_double   (cellSize, 5e-9, "mesh cell size (dx,dy,dz)");
    // DEFINE_int32    (zdim, 3, "number of mesh points in Z-direction");
    DEFINE_string   (Minit_file, "", "Text file in matrix format for intial state of M");
    DEFINE_int32    (P, 3, "order of multipole expansion");
    DEFINE_double   (finaltime, 10e-9, "final time for simulation to terminate");
    DEFINE_double   (timestep, 1e-12, "simulation time step");
    DEFINE_bool     (demag, true, "account for demagnetization field");
    DEFINE_bool     (exchange, true, "account for exchange field");
    DEFINE_bool     (external, false, "account for external field");
    DEFINE_bool     (use_fmm, false, "use Fast-Multipole-Method for potential calculation");
    DEFINE_bool     (use_gpu, false, "use Graphics Processor for potential calculation");
    DEFINE_string   (sim_name, "sim_untitled", "name of this simulation");
    DEFINE_int32    (seed, time(NULL), "seed for random number generator");

//*************************************************************************//
//******************** Main function **************************************//
//*************************************************************************//
int main(int argc, char **argv)
{
    int status = 0;
    const int verbosity = 4;

    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);

// read command line arguments
    google::ParseCommandLineFlags(&argc, &argv, true);

    // char *maskImg = (char*)FLAGS_maskImg.c_str();
    // char *maskTxt = (char*)FLAGS_maskTxt.c_str();
    unsigned int P = FLAGS_P;
    fptype finaltime = FLAGS_finaltime;
    fptype timestep = FLAGS_timestep;
    fptype cellSize  = FLAGS_cellSize;
    int Nx = FLAGS_Nx;
    int Ny = FLAGS_Ny;
    int Nlayers = FLAGS_Nlayers;
    // fptype sample_width  = FLAGS_sample_width;
    // fptype sample_height = FLAGS_sample_height;
    // fptype sample_depth  = FLAGS_sample_depth;
    // int xdim = FLAGS_xdim;
    // int ydim = FLAGS_ydim;
    // int zdim = FLAGS_zdim;
    int demag = FLAGS_demag;
    int exchange = FLAGS_exchange;
    int external = FLAGS_external;
    int use_fmm = FLAGS_use_fmm;
    int use_gpu = FLAGS_use_gpu;
    char *sim_name = (char*)FLAGS_sim_name.c_str();

    assert(P <= 4); // very important
    assert(Nlayers == 3); // very important
    srand(FLAGS_seed);

// print command line arguments
#ifdef _OPENMP
    printf("Compiled with OpenMP and running with %s threads.\n", getenv("OMP_NUM_THREADS"));
#endif
    if(verbosity >= 2) {
        printf("Minit_file = %s \n", FLAGS_Minit_file.c_str());
        printf("P = %d \n", P);
        printf("finaltime = %g \n", finaltime);
        printf("timestep = %g \n", timestep);
        printf("Nlayers = %d \n", Nlayers);
        printf("demag = %d \n", demag);
        printf("exchange = %d \n", exchange);
        printf("external = %d \n", external);
        printf("use_fmm = %d \n", use_fmm);
        printf("use_gpu = %d \n", use_gpu);
        printf("sim_name = %s \n", sim_name);
        printf("SEED = %d \n", FLAGS_seed);
    }

// create directory to hold results
    int err = mkdir(sim_name, 0755);
    // if(err == EEXIST)
        // printf("ERROR%d: Directory already exists\n", err);
    // else
        // printf("ERROR%d: Some other error\n", err);
    // sprintf(command,"copy %s %s",source_file,destination_file);
    // system(command);
    // if(status) return EXIT_FAILURE;

// Material parameters
// ================================================
    const fptype mu_0 = 4 * M_PI * 1e-7; // permeability of vacuum
    const fptype Ms = 8.6e5;             // saturation magnetization (permalloy)
    const fptype Aexch = 1.3e-11;        // exchange constant (permalloy)
    const fptype alfa = 0.008;             // damping coefficient (permalloy)
    const fptype gamma = 2.21e5;         // gyromagnetic ratio (permalloy)

// Mask configuration for magnetic material
    int xdim = Nx;
    int ydim = Ny;
    int zdim = Nlayers; // +2
    assert(zdim == 3); // very important
    assert(cellSize <= 5e-9);
    const fptype dx = cellSize;
    const fptype dy = cellSize;
    const fptype dz = cellSize;
    const fptype sample_width =  dx * xdim;
    const fptype sample_height =  dy * ydim;
    const fptype sample_depth =  dz * zdim;
    // assert(dx == dy);
    // assert(dy == dz);
    printf("(xdim, ydim, zdim) = (%d, %d, %d)\n", xdim, ydim, zdim);
    printf("(sample_width, sample_height, sample_depth) = (%g, %g, %g)\n", sample_width, sample_height, sample_depth);
    printf("(dx, dy, dz) = (%g, %g, %g)\n", dx, dy, dz);

// load Minitial from file
    Vector3 *Minit = new Vector3[ydim*xdim]();  // magnetization matrix
    if(Minit == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    status |= load_Minit(FLAGS_Minit_file.c_str(), ydim, xdim, Minit, verbosity);
    if(status) return EXIT_FAILURE;


// generate the initial magnetization distribution
    byte *material = new byte[zdim*ydim*xdim]();  // material matrix
    Vector3 *M = new Vector3[zdim*ydim*xdim]();  // magnetization matrix
    if(material == NULL || M == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

    assert(zdim == 3); // assert(zdim >= 3);
    for(int z = 1; z < zdim-1; z++) {
        for(int y = 0; y < ydim; y++) {
            for(int x = 0; x < xdim; x++) {
                if(1) // mask[y*xdim + x])
                {
                    Vector3 M1 = Minit[0*ydim*xdim + y*xdim + x];
                    M[z*ydim*xdim + y*xdim + x] = M1;
                    if(M1.x || M1.y || M1.z)
                        material[z*ydim*xdim + y*xdim + x] = 1;
                }
            }
        }
    }

// write M field to file
    char filename[1000];
    sprintf(filename, "%s/%s", sim_name, "Minit.dat");
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


    // return EXIT_FAILURE;

// magnetization dynamics
// ===================================================================
    // return EXIT_FAILURE;
    status |= time_marching(    material, M,
                                finaltime, timestep,
                                xdim, ydim, zdim, dx, dy, dz,
                                P, mu_0, Ms, Aexch, alfa, gamma,
                                demag, exchange, external, use_fmm,
                                use_gpu, sim_name, verbosity );

// closing
    delete []Minit;
    delete []M;
    delete []material;

    printf("SEED = %d\n", FLAGS_seed);
    // printf("%s\n", status ? "failed to complete" : "successfuly completed");

    status |= gettimeofday(&time2, NULL);
    double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    printf("Simulation completed in %f seconds.\n", deltatime);

    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}




/*
#ifdef USE_FREEIMAGE
    BYTE *mask = NULL; // mask matrix
    if(strlen(maskImg) != 0)
        load_mask(maskImg, &mask, &xdim, &ydim);
    else
        // load_mask_txt(maskTxt, &mask, &xdim, &ydim);
#else
    printf("Please compile with USE_FREEIMAGE directive\n");
    return EXIT_FAILURE;
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
*/

/*
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

    assert(zdim >= 3);
    for(int z = 1; z < zdim-1; z++) {
        for(int y = 0; y < ydim; y++) {
            for(int x = 0; x < xdim; x++) {
                if(mask[y*xdim + x])
                {
                    fptype theta = M_PI/2;
                    fptype phi   = 0;
                    // fptype theta = frand_atob(0, 180) * M_PI/180;
                    // fptype phi   = frand_atob(0, 360) * M_PI/180;
                    // fptype theta = M_PI/2 + frand_atob(-10, 10) * M_PI/180;
                    // fptype phi   = 0 + frand_atob(-90, 90) * M_PI/180;
                    if(IC_singledomain) {
                        theta = 90 * M_PI/180;
                        phi   = 45 * M_PI/180;
                        // theta = 0 + frand_atob(-1, 1) * M_PI/180;
                        // phi   = frand_atob(0, 360) * M_PI/180;
                        // theta = (90 + frand_atob(-1, 1)) * M_PI/180;
                        // phi = (90 + frand_atob(-1, 1)) * M_PI/180;
                    }
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
*/
