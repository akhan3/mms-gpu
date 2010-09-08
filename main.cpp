#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <assert.h>
#include <iostream>
#include <cmath>
#include <time.h>
#include <sys/time.h>
#include "Box.hpp"
#include "Cmpx.hpp"
#include "Vector3.hpp"
#include "helper_functions.hpp"
#include "vector_functions.hpp"
// #include "ode_functions.hpp"
using std::cout;
using std::endl;

int Hfield (    const Vector3 *M, Vector3 *H, fptype *charge, fptype *potential,
                const int xdim, const int ydim, const int zdim, const fptype meshwidth,
                const fptype mu_0, const fptype Ms, const fptype Aexch,
                const int demag, const int exchange, const int external, const int use_fmm, const int P,
                const int verbose_level );

//*************************************************************************//
//******************** Main function **************************************//
//*************************************************************************//
int main(int argc, char **argv)
{
    int status = 0;
    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);
    const int verbose_level = 2;

// intial command line arguments
    char filename_arg[1000] = "verysmall_16x16.png";
    unsigned int P = 3;
    fptype finaltime = 1e-9;
    fptype timestep = 1e-14;
    fptype meshwidth = 1e-9;
    unsigned int xdim = 16;
    unsigned int ydim = 16;
    unsigned int zdim = 1;
    int demag = true;
    int exchange = true;
    int external = false;
    int use_fmm = false;
    char sim_name[1000] = "sim_untitled";
    unsigned int seed = time(NULL);
// read command line arguments
    if(argc >= 2)
        sscanf(argv[1], "%s", filename_arg);
    if(argc >= 3) {
        sscanf(argv[2], "%u", &P);
        assert(P <= 4);
    }
    if(argc >= 4)
        sscanf(argv[3], "%lf", &finaltime);
    if(argc >= 5)
        sscanf(argv[4], "%lf", &timestep);
    if(argc >= 6)
        sscanf(argv[5], "%lf", &meshwidth);
    if(argc >= 7)
        sscanf(argv[6], "%d", &xdim);
    if(argc >= 8)
        sscanf(argv[7], "%d", &ydim);
    if(argc >= 9)
        sscanf(argv[8], "%d", &zdim);
    if(argc >= 10)
        sscanf(argv[9], "%d", &demag);
    if(argc >= 11)
        sscanf(argv[10], "%d", &exchange);
    if(argc >= 12)
        sscanf(argv[11], "%d", &external);
    if(argc >= 13)
        sscanf(argv[12], "%d", &use_fmm);
    if(argc >= 14)
        sscanf(argv[13], "%s", sim_name);
    if(argc >= 15)
        sscanf(argv[14], "%u", &seed);
    srand(seed);
// print command line arguments
    printf("imagefile = %s \n", filename_arg);
    printf("P = %d \n", P);
    printf("finaltime = %g \n", finaltime);
    printf("timestep = %g \n", timestep);
    printf("meshwidth = %g \n", meshwidth);
    printf("xdim = %d \n", xdim);
    printf("ydim = %d \n", ydim);
    printf("zdim = %d \n", zdim);
    printf("demag = %d \n", demag);
    printf("exchange = %d \n", exchange);
    printf("external = %d \n", external);
    printf("use_fmm = %d \n", use_fmm);
    printf("sim_name = %s \n", sim_name);
    printf("SEED = %d \n", seed);

// Material parameters
// ================================================
    const fptype mu_0 = 4 * M_PI * 1e-7; // permeability of vacuum
    const fptype Ms = 8.6e5;             // saturation magnetization (permalloy)
    const fptype Aexch = 1.3e-11;        // exchange constant (permalloy)
    const fptype alfa = 0.5;             // damping coefficient (permalloy)
    const fptype gamma = 2.21e5;         // gyromagnetic ratio (permalloy)

// Mask configuration for magnetic material
#ifdef USE_FREEIMAGE
    BYTE *mask = NULL; // mask matrix
    char filename[1000];
    sprintf(filename, "%s", filename_arg);
// read the mask from file
    load_mask(filename, &mask, &xdim, &ydim);
    zdim = 1;
#else
    byte *mask = new byte[ydim*xdim](); // mask matrix
// specimen magnet 20x20x20
    for(unsigned int y = 0; y < ydim; y++)
        for(unsigned int x = 0; x < xdim; x++)
            mask[y*xdim + x] = 1;   // all white (no material)
    // for(unsigned int y = 21; y <= 41; y++)
        // for(unsigned int x = 21; x <= 41; x++)
    for(unsigned int y = 1; y < ydim-1; y++)
        for(unsigned int x = 1; x < xdim-1; x++)
            mask[y*xdim + x] = 0;   // selected black (material)
#endif
    // assert(xdim == ydim);
    printf("(xdim, ydim, zdim) = (%d, %d, %d)\n", xdim, ydim, zdim);

// generate the initial magnetization distribution
    byte *material = new byte[zdim*ydim*xdim]();  // material matrix
    Vector3 *M = new Vector3[zdim*ydim*xdim]();  // magnetization matrix
    if(material == NULL || M == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

    assert(zdim == 1);
    for(unsigned int z = 0; z < zdim; z++) {
        for(unsigned int y = 0; y < ydim; y++) {
            for(unsigned int x = 0; x < xdim; x++) {
                if(!mask[y*xdim + x])
                {
                    fptype theta = frand_atob(0, 180) * M_PI/180;
                    fptype phi   = frand_atob(0, 360) * M_PI/180;
                    // fptype theta = M_PI/2;
                    // fptype phi   = 0;
                    M[z*ydim*xdim + y*xdim + x] = Ms * Vector3(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
                    material[z*ydim*xdim + y*xdim + x] = 1;
                }
            }
        }
    }
    delete []mask;

// write M field to file
    status |= save_vector3d(M, zdim, ydim, xdim, "M.dat", verbose_level);
    if(status) return EXIT_FAILURE;

// write material field to file
    fptype *m = new fptype[xdim*ydim];
    int z = 0;
    for(unsigned int y = 0; y < ydim; y++)
        for(unsigned int x = 0; x < xdim; x++)
            m[y*xdim + x] = (fptype)material[z*ydim*xdim + y*xdim + x];
    status |= matrix2file(m, ydim, xdim, "material.dat", 100);
    if(status) return EXIT_FAILURE;
    delete []m;

// magnetization dynamics
// ===================================================================
    const int xyzdim = zdim*ydim*xdim;
    Vector3 *H_fmm      = new Vector3[xyzdim]();
    Vector3 *H_exact    = new Vector3[xyzdim]();
    fptype *charge      = new fptype[xyzdim]();
    fptype *potential   = new fptype[xyzdim]();
    if(H_fmm == NULL && H_exact == NULL && charge == NULL && potential == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

    status |= Hfield(M, H_fmm,   charge, potential, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, demag, exchange, external, 1||use_fmm, P, verbose_level+200);
    if(status) return EXIT_FAILURE;
    status |= save_vector3d(H_fmm, zdim, ydim, xdim, "H_fmm.dat", verbose_level);

    // status |= Hfield(M, H_exact, charge, potential, xdim, ydim, zdim, meshwidth, mu_0, Ms, Aexch, demag, exchange, external, 0&&use_fmm, P, verbose_level+200);
    // if(status) return EXIT_FAILURE;
    // status |= save_vector3d(H_exact, zdim, ydim, xdim, "H_exact.dat", verbose_level);


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
