#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <cmath>
// #include <cutil_inline.h>
#include "mydefs.hpp"
#include "Box.hpp"
// #include "Queue.hpp"
#include "Cmpx.hpp"
#include "Vector3.hpp"
#include "numerics.hpp"
#include "helper_functions.hpp"

#include "potential_calc.cu"

// FMM interaction evaluation kernel
// ===================================
__global__
void fmm_kernel(    const Box *const n_gmem,
                    const Box *const ni_gmem,
                    const Cmpx *const mpc_gmem,
                    fptype *potential_gmem,
                    const unsigned int limit,
                    const int P,
                    const int xdim,
                    const int ydim,
                    const int zdim,
                    const int zc
                )
{
    int bi = blockIdx.x;
    int ti = threadIdx.x;
    if(ti == 0) {
        const Box n = *n_gmem;
        fptype width = powf(2, limit - n.level);
        // block index as interaction loop
        if(n.interaction[bi] != NULL) {
            Box ni = ni_gmem[bi];
            for(int yy = ceilf(ni.cy-width/2); yy <= floorf(ni.cy+width/2); yy++) {
                for(int xx = ceilf(ni.cx-width/2); xx <= floorf(ni.cx+width/2); xx++) {
                    for (int zp = 0; zp < zdim; zp++) { // for each potential layer in zdim
                        Vector3 r(xx-n.cx, yy-n.cy, zp-zc);
                        Cmpx sum_over_lm;
                        for(int l=0; l<=P; l++) {
                            Cmpx sum_over_m;
                            for(int m=-l; m<=l; m++) {
                                Cmpx sph = spherical_harmonic(l, m, r.colatitude(), r.azimuth());
                                sph *= (1.0*factorial(l-fabsf(m))) / factorial(l+fabsf(m));
                                sph *= mpc_gmem[l*l+l+m];
                                sum_over_m += sph;
                            }
                            sum_over_m *= 1 / powf(r.magnitude(), l+1);
                            sum_over_lm += sum_over_m;
                        }
                        potential_gmem[zp*ydim*xdim + yy*xdim + xx] += sum_over_lm.get_re();
                    }
                }
            }
        } // if(ni != NULL)
    } // if thread 0
}



// wrapper function for FMM kernel
// =======================================================
int fmm_gpu(        const Box *const n_ptr,
                    const Cmpx *const mpc,
                    fptype *potential,
                    const unsigned int limit,
                    const int P,    // multipole series truncation (l = 0...P)
                    const int xdim, const int ydim, const int zdim,
                    const int zc,   // charge layer
                    const int use_gpu,
                    const int verbose_level
                )
{
    int status = 0;
    static int first_time = 1;

    // concatenate 27 interaction boxes in CPU memory
    static Box *ni_cpumem = NULL;
    if(first_time) {
        ni_cpumem = (Box*)malloc(27 * sizeof(Box));
        if(ni_cpumem == NULL) {
            fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
            return EXIT_FAILURE;
        }
    }
    for(int bi = 0; bi < 27; bi++) {
        if(n_ptr->interaction[bi] != NULL)
            memcpy(ni_cpumem + bi, n_ptr->interaction[bi], sizeof(Box));
    }

    // set up device memory pointers
    static Box *n_gmem = NULL;
    static Box *ni_gmem = NULL;
    static Cmpx *mpc_gmem = NULL;
    static fptype *potential_gmem = NULL;

    if(first_time) {
        // select device to use
        // cudaSetDevice(1);
        // allocate memory on device
        cudaMalloc((void**)&n_gmem, sizeof(Box));         checkCUDAError("Allocate n_gmem");
        cudaMalloc((void**)&ni_gmem, 27 * sizeof(Box));
        cudaMalloc((void**)&mpc_gmem, (P+1)*(P+1) * sizeof(Cmpx));
        cudaMalloc((void**)&potential_gmem, zdim*ydim*xdim * sizeof(fptype));
        if(n_gmem == NULL || ni_gmem == NULL || mpc_gmem == NULL || potential_gmem == NULL) {
            fprintf(stderr, "%s:%d Error allocating memory on GPU\n", __FILE__, __LINE__);
            return EXIT_FAILURE;
        }
    }

    // copy required data to GPU global memory
    cudaMemcpy(n_gmem, n_ptr, sizeof(Box), cudaMemcpyHostToDevice);
    checkCUDAError("Copying n_gmem");
    cudaMemcpy(ni_gmem, ni_cpumem, 27 * sizeof(Box), cudaMemcpyHostToDevice);
    checkCUDAError("Copying ni_cpumem");
    cudaMemcpy(mpc_gmem, mpc, (P+1)*(P+1) * sizeof(Cmpx), cudaMemcpyHostToDevice);
    checkCUDAError("Copying mpc_gmem");
    cudaMemcpy(potential_gmem, potential, zdim*ydim*xdim * sizeof(fptype), cudaMemcpyHostToDevice);
    checkCUDAError("Copying potential_gmem");

    // int currentDevice;
    // cudaGetDevice(&currentDevice);
    // printf("using device %d\n", currentDevice);

    // set up kernel parameters
    #define MAXTHREADSPERBLOCK    1024
    int problem_size = xdim*ydim;
    // dim3 grid = ceil(total_threads / (fptype)MAXTHREADSPERBLOCK);
    // dim3 threads(MAXTHREADSPERBLOCK, 1, 1);
    dim3 grid = 27;
    const int stride = ceil(problem_size / (fptype)MAXTHREADSPERBLOCK);
    dim3 threads = 32;
    assert(threads.x <= MAXTHREADSPERBLOCK);    // max_threads_per_block

    if(first_time) {
        // printf("x=%u, y=%u, threads.x=%u, threads.y=%u, threads.z=%u, stride=%u, grid.x=%u, grid.y=%u, grid.z=%u\n",
                // xdim, ydim, threads.x,    threads.y,    threads.z,    stride,    grid.x,    grid.y,    grid.z);
        printf("\n\nlaunching kernel with %u blocks and %u threads...\n\n",
                    grid.x*grid.y*grid.z, threads.x*threads.y*threads.z);
        fflush(NULL);
    }

    // start timer
    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);

    // launch the kernel
    fmm_kernel <<<grid, threads>>>
        (n_gmem, ni_gmem, mpc_gmem, potential_gmem, limit, P, xdim, ydim, zdim, zc);
    checkCUDAError("Exeuting Kernel fmm_kernel()");
    cudaThreadSynchronize();

    // read the timer
    status |= gettimeofday(&time2, NULL);
    // double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    // printf("  Kernel completed in %f seconds.\n", deltatime);

    // copy potential (the result of kernel) to host main memory
    cudaMemcpy(potential, potential_gmem, zdim*ydim*xdim * sizeof(fptype), cudaMemcpyDeviceToHost);
    checkCUDAError("Copying potential_gmem");
    // cudaFree(charge_gmem);
    // cudaFree(potential_gmem);
    cudaThreadSynchronize();

    first_time = 0;
    return EXIT_SUCCESS;
}
