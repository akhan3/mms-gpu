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
DEFINE_string   (simName, "simUntitled", "name of this simulation");
DEFINE_double   (Ms, 8.6e5, "saturation magnetization [A/m]");
DEFINE_double   (Aexch, 1.3e-11, "exchange constant [?]");
DEFINE_double   (alpha, 0.008, "damping coefficient []");
DEFINE_double   (gamma, 2.21e5, "gyromagnetic ratio [1/(A/m/s)]");
// DEFINE_int32    (Nx, 13, "number of mesh points in X-direction");
// DEFINE_int32    (Ny, 13, "number of mesh points in Y-direction");
DEFINE_int32    (Nlayers, 3, "number of layers of magnetic materials");
DEFINE_double   (cellSize, 5e-9, "mesh cell size (dx,dy,dz) [m]");
DEFINE_string   (MinitFile, "disc_13x13.txt", "Text file in matrix format for intial state of M");
DEFINE_double   (finaltime, 10e-9, "final time for simulation to terminate [s]");
DEFINE_double   (timestep, 1e-12, "simulation time step [s]");
DEFINE_int32    (subsample_demag, 1, "Subsampling of Hdemag calculation");
DEFINE_bool     (silent_stdout, false, "No printing on stdout. Useful for batch runs");
DEFINE_bool     (demag, true, "account for demagnetization field");
DEFINE_bool     (exchange, true, "account for exchange field");
DEFINE_bool     (external, false, "account for external field");
DEFINE_bool     (useGPU, false, "use Graphics Processor for potential calculation");
DEFINE_int32    (cudaDevice, 1, "Tesla card number to use");
DEFINE_bool     (useFMM, false, "use Fast-Multipole-Method for potential calculation");
DEFINE_int32    (fmmP, 3, "order of multipole expansion");
DEFINE_int32    (seed, time(NULL), "seed for random number generator");
DEFINE_int32    (verbosity, 4, "Verbosity level of the simulator");
DEFINE_bool     (printArgsAndExit, false, "Print command line arguments and exit");

//*************************************************************************//
//******************** Main function **************************************//
//*************************************************************************//
int main(int argc, char **argv)
{
    int status = 0;
    // const int verbosity = 4;

    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);

// read command line arguments
    google::ParseCommandLineFlags(&argc, &argv, true);

// sanity check
    assert(FLAGS_fmmP <= 4); // very important
    assert(FLAGS_Nlayers == 3); // very important
    assert(FLAGS_cellSize <= 5e-9);

// seed the random generator
    srand(FLAGS_seed);

// convert the strings
    char *simName = (char*)FLAGS_simName.c_str();
    char *MinitFile = (char*)FLAGS_MinitFile.c_str();

// print command line arguments
#ifdef _OPENMP
    printf("Compiled with OpenMP and running with %s threads.\n", getenv("OMP_NUM_THREADS"));
#endif
    if(FLAGS_verbosity >= 2 || FLAGS_printArgsAndExit) {
        printf("sim_name = %s \n", simName);

        // printf("Nx = %d \n", FLAGS_Nx);
        // printf("Ny = %d \n", FLAGS_Ny);
        printf("Nlayers = %d \n", FLAGS_Nlayers);
        printf("cellSize = %g \n", FLAGS_cellSize);
        printf("timestep = %g \n", FLAGS_timestep);
        printf("finaltime = %g \n", FLAGS_finaltime);
        printf("terminatingTorque = %g \n", FLAGS_terminatingTorque);
        printf("adjust_step = %d \n", FLAGS_adjust_step);
        NEWLINE;
        printf("Minit_file = %s \n", MinitFile);
        printf("Ms = %g \n", FLAGS_Ms);
        printf("alpha = %g \n", FLAGS_alpha);
        printf("gamma = %g \n", FLAGS_gamma);
        printf("Aexch = %g \n", FLAGS_Aexch);
        NEWLINE;
        printf("Bext = %s \n", FLAGS_Bext.c_str());
        printf("STO_JFile = %s \n", FLAGS_STO_JFile.c_str());
        printf("STO_I = %g \n", FLAGS_STO_I);
        printf("STO_Pdir = %s \n", FLAGS_STO_Pdir.c_str());
        // printf("STO_A = %g \n", FLAGS_STO_A);
        printf("STO_P = %g \n", FLAGS_STO_P);
        printf("STO_Lambda = %g \n", FLAGS_STO_Lambda);
        printf("STO_t0 = %g \n", FLAGS_STO_t0);
        NEWLINE;
        printf("demag = %d \n", FLAGS_demag);
        printf("subsample_demag = %d \n", FLAGS_subsample_demag);
        printf("exchange = %d \n", FLAGS_exchange);
        printf("external = %d \n", FLAGS_external);
        printf("useGPU = %d \n", FLAGS_useGPU);
        printf("cudaDevice = %d \n", FLAGS_cudaDevice);
        printf("useFMM = %d \n", FLAGS_useFMM);
        printf("fmmP = %d \n", FLAGS_fmmP);
        NEWLINE;
        printf("log_Mfield = %d \n", FLAGS_log_Mfield);
        printf("subsample = %d \n", FLAGS_subsample);
        NEWLINE;
        printf("SEED = %d \n", FLAGS_seed);
        printf("verbosity = %d \n", FLAGS_verbosity);
        printf("printArgsAndExit = %d \n", FLAGS_printArgsAndExit);
        NEWLINE;
    }

    if(FLAGS_printArgsAndExit)
        return EXIT_SUCCESS;

// create directory to hold results
    int err = mkdir(simName, 0755);
    // if(err == EEXIST)
        // printf("ERROR%d: Directory already exists\n", err);
    // else
        // printf("ERROR%d: Some other error\n", err);
    // sprintf(command,"copy %s %s",source_file,destination_file);
    // system(command);
    // if(status) return EXIT_FAILURE;

// Material parameters
// ================================================
    // const fptype Ms = FLAGS_Ms;
    // const fptype Aexch = FLAGS_Aexch;
    // const fptype alpha = FLAGS_alpha;
    // const fptype gamma = FLAGS_gamma;

// Mask configuration for magnetic material
    int xdim;   // = FLAGS_Nx;
    int ydim;   // = FLAGS_Ny;
    int zdim = FLAGS_Nlayers; // +2
    const fptype dx = FLAGS_cellSize;
    const fptype dy = FLAGS_cellSize;
    const fptype dz = FLAGS_cellSize;

// load Minitial from file
    Vector3 *Minit = NULL;  // magnetization matrix. will be allocated inside the following function
    status |= load_Minit(MinitFile, &ydim, &xdim, &Minit, FLAGS_verbosity);
    if(status) return EXIT_FAILURE;

// print some info about geometry
    const fptype sample_width =  dx * xdim;
    const fptype sample_height =  dy * ydim;
    const fptype sample_depth =  dz * zdim;
    // assert(dx == dy);
    // assert(dy == dz);
    printf("(xdim, ydim, zdim) = (%d, %d, %d)\n", xdim, ydim, zdim);
    printf("(sample_width, sample_height, sample_depth) = (%g, %g, %g) nm\n", sample_width/1e-9, sample_height/1e-9, sample_depth/1e-9);
    printf("(dx, dy, dz) = (%g, %g, %g) nm\n", dx/1e-9, dy/1e-9, dz/1e-9);

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
    sprintf(filename, "%s/%s", simName, "Minit.dat");
    // status |= save_vector3d(M, zdim, ydim, xdim, filename, verbosity);
    status |= save_vector3d(&M[ydim*xdim], 1, ydim, xdim, filename, FLAGS_verbosity);
    if(status) return EXIT_FAILURE;

// write material field to file
    fptype *m = new fptype[xdim*ydim];
    int z = 1;
    for(int y = 0; y < ydim; y++)
        for(int x = 0; x < xdim; x++)
            m[y*xdim + x] = (fptype)material[z*ydim*xdim + y*xdim + x];
    sprintf(filename, "%s/%s", simName, "material.dat");
    status |= matrix2file(m, ydim, xdim, filename, FLAGS_verbosity);
    if(status) return EXIT_FAILURE;
    delete []m;


    // return EXIT_FAILURE;

// magnetization dynamics
// ===================================================================
    status |= time_marching(    material, M,
                                FLAGS_finaltime, FLAGS_timestep,
                                xdim, ydim, zdim, dx, dy, dz,
                                FLAGS_fmmP, mu_0, FLAGS_Ms, FLAGS_Aexch, FLAGS_alpha, FLAGS_gamma,
                                FLAGS_demag, FLAGS_exchange, FLAGS_external, FLAGS_useFMM,
                                FLAGS_useGPU, simName, FLAGS_verbosity );

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
