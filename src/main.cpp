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
DEFINE_string   (simName, "simUntitled", "Name of this simulation");
DEFINE_string   (simDesc, "Untitled simulation", "Description about this simulation");
DEFINE_double   (Ms, 8.6e5, "Saturation magnetization [A/m]");
DEFINE_double   (Aexch, 1.3e-11, "Exchange constant [?]");
DEFINE_double   (alpha, 0.008, "Damping coefficient []");
DEFINE_double   (gamma, 2.21e5, "Gyromagnetic ratio [1/(A/m/s)]");
DEFINE_int32    (Nlayers, 3, "Number of layers of magnetic materials");
DEFINE_double   (cellSize, 5e-9, "Mesh cell size (dx,dy,dz) [m]");
DEFINE_string   (MinitFile, "disc_13x13.txt", "Text file in matrix format for intial state of M");
DEFINE_double   (finaltime, 10e-9, "Final time for simulation to terminate [s]");
DEFINE_double   (timestep, 1e-12, "Simulation time step [s]");
DEFINE_int32    (subsample_demag, 1, "Subsampling of Hdemag calculation");
DEFINE_bool     (silent_stdout, false, "No printing on stdout. Useful for batch runs");
DEFINE_bool     (demag, true, "Account for demagnetization field");
DEFINE_bool     (exchange, true, "Account for exchange field");
DEFINE_bool     (external, false, "Account for external field");
DEFINE_bool     (useGPU, false, "Use Graphics Processor for potential calculation");
DEFINE_int32    (cudaDevice, 1, "Tesla card number to use");
DEFINE_bool     (useFMM, false, "Use Fast-Multipole-Method for potential calculation");
DEFINE_int32    (fmmP, 3, "Order of multipole expansion");
DEFINE_int32    (seed, time(NULL), "Seed for random number generator. Default value is current time.");
DEFINE_int32    (verbosity, 4, "Verbosity level of the simulator");
DEFINE_bool     (printArgsAndExit, false, "Print command line arguments and exit");


//*************************************************************************//
//******************** Main function **************************************//
//*************************************************************************//
int main(int argc, char **argv)
{
    int status = 0;
    char filename[1000];
    // const int verbosity = 4;

    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);

// read command line arguments
    NEWLINE;
    printf("INFO: Parsing command line flags... ");
    google::ParseCommandLineFlags(&argc, &argv, true);
    printf("done!\n");

// sanity check
    assert(FLAGS_fmmP <= 4); // very important
    assert(FLAGS_Nlayers == 3); // very important
    assert(FLAGS_cellSize <= 5e-9);

// seed the random generator
    srand(FLAGS_seed);

// convert the strings
    char *simName = (char*)FLAGS_simName.c_str();
    char *MinitFile = (char*)FLAGS_MinitFile.c_str();
    time_t now;
    time(&now);
    printf("#===============================================================================\n");
    printf("# Simulation title : %s\n", simName);
    printf("# Description      : %s\n", FLAGS_simDesc.c_str());
    printf("# Time started     : %s", ctime(&now));
    printf("#===============================================================================\n");

// create directory to hold results
    {
        int status = mkdir(simName, S_IRWXU);
        if(!status)
            printf("INFO: Directory \"%s\" created.\n", simName);
        else {
            // printf("mkdir caused errno = %d\n", errno);
            if(errno == EEXIST)
                printf("INFO: Directory \"%s\" already exists and will be overwritten!\n", simName);
            else if(errno == EACCES) {
                printf("ERROR: Write permission is denied for the parent directory in which the new directory is to be added.\n");
                return EXIT_FAILURE;
            }
            else if(errno == ENOSPC) {
                printf("ERROR: The file system doesn't have enough room to create the new directory.\n");
                return EXIT_FAILURE;
            }
            else if(errno == EROFS) {
                printf("ERROR: The parent directory of the directory being created is on a read-only file system and cannot be modified.\n");
                return EXIT_FAILURE;
            }
            else if(errno == EMLINK) {
                printf("ERROR: The parent directory has too many links (entries).\n");
                return EXIT_FAILURE;
            }
        }
    }

// write flagfile.mif
    sprintf(filename, "%s/parameters-%s.mif", simName, simName);
    FILE *ff = fopen(filename, "w");
    if(ff == NULL) {
        fprintf(stderr, "FATAL ERROR: Error opening file %s\n", filename);
        return EXIT_FAILURE;
    }
    fprintf(ff, "#===============================================================================\n");
    fprintf(ff, "# Simulation title : %s\n", simName);
    fprintf(ff, "# Description      : %s\n", FLAGS_simDesc.c_str());
    fprintf(ff, "# Time started     : %s", ctime(&now));
    fprintf(ff, "#===============================================================================\n");
    fprintf(ff, "\n\n");
    std::vector<google::CommandLineFlagInfo> allFlags;
    google::GetAllFlags(&allFlags);
    for (std::vector<google::CommandLineFlagInfo>::iterator flag = allFlags.begin(); flag != allFlags.end(); ++flag) {
        if(flag->name == "flagfile")
            continue;
        if(flag->name == "fromenv") // built-in flags from this point
            break;
        fprintf(ff, "# %s ", flag->description.c_str());
        fprintf(ff, "# \"%s\"\n", flag->default_value.c_str());
        fprintf(ff, "\t-%s=%s\n\n", flag->name.c_str(), flag->current_value.c_str());
    }
    // Don't close here! Keep opened for appending finishing time!
    // fclose(ff);

    if(FLAGS_printArgsAndExit)
        return EXIT_SUCCESS;

#ifdef _OPENMP
    printf("INFO: Compiled with OpenMP and running with %s threads.\n", getenv("OMP_NUM_THREADS"));
#endif

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
    printf("INFO: (dx, dy, dz) = (%g, %g, %g) nm\n", dx/1e-9, dy/1e-9, dz/1e-9);
    printf("INFO: (sample_width, sample_height, sample_depth) = (%g, %g, %g) nm\n", sample_width/1e-9, sample_height/1e-9, sample_depth/1e-9);
    printf("INFO: (xdim, ydim, zdim) = (%d, %d, %d)\n", xdim, ydim, zdim);
    printf("INFO: (timestep, finaltime, steps) = (%g, %g). %d steps\n", FLAGS_timestep, FLAGS_finaltime, (int)(FLAGS_finaltime/FLAGS_timestep));

// Check to see if xdim and ydim are in powers of 2
    if(FLAGS_useFMM) {
        assert(xdim == ydim);
        double logN = log2f(xdim);
        printf("pow(2,logN) = %d\n", (int)powf(2,logN));
        printf("pow(2,ceil(logN)) = %d\n", (int)powf(2,ceil(logN)));
        assert(xdim == pow(2,logN));
    }

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

    status |= gettimeofday(&time2, NULL);
    double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    int hours = deltatime / 3600;
    int mins = (deltatime - 3600*hours) / 60;
    int secs = deltatime - 3600*hours - 60*mins;
    time(&now);
    printf("#===============================================================================\n");
    printf("# Time finished     : %s", ctime(&now));
    printf("# Time taken        : %d:%d:%d (h:mm:ss)\n", hours, mins, secs);
    printf("#===============================================================================\n");

    fprintf(ff, "\n");
    fprintf(ff, "#===============================================================================\n");
    fprintf(ff, "# Time finished     : %s", ctime(&now));
    fprintf(ff, "# Time taken        : %d:%d:%d (h:mm:ss)\n", hours, mins, secs);
    fprintf(ff, "#===============================================================================\n");
    fclose(ff);

    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}
