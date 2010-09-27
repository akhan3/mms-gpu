#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <cmath>
// #include <cutil_inline.h>
#include "mydefs.hpp"
#include "Box.hpp"
#include "Cmpx.hpp"
#include "Vector3.hpp"
#include "numerics.hpp"
#include "helper_functions.hpp"

#include "potential_calc.cu"
#include "fmm_calc.cu"


// FMM interaction evaluation kernel
// ===================================
__global__
void fmm_kernel
    (   const Box *const n_gmem,
        const Box *const ni_gmem,
        const Cmpx *const mpc_gmem,
        fptype *potential_gmem,
        const unsigned int limit,
        const int P,
        const int xdim,
        const int ydim,
        const int zdim,
        const int zc    )
{
    int bi = blockIdx.x;
    int ti = threadIdx.x;

    __shared__ Box n;
    __shared__ Box ni;
    __shared__ int x1;
    __shared__ int y1;
    __shared__ int width;
    __shared__ int stride;

    if(ti == 0) {
        n = *n_gmem;
    }
    __syncthreads();

    if(n.interaction[bi] == NULL)
        return;

    if(ti == 0) {
        // block index as interaction loop
        ni = ni_gmem[bi];
        width = powf(2, limit - n.level);
        x1 = ceilf(ni.cx - width/2);
        y1 = ceilf(ni.cy - width/2);
        stride = ceilf(width*width / (fptype)blockDim.x);
    }
    __syncthreads();

    // int ti1 = ti * stride;
    if(ti * stride < width*width) {
        for(int i = ti * stride; i < ti * stride + stride; i++) {
            if(i < width*width) {
                int x = x1 + i % width;
                int y = y1 + (i - (x - x1)) / width;
                for (int zp = 0; zp < zdim; zp++) { // for each potential layer in zdim
                    Vector3 r(x-n.cx, y-n.cy, zp-zc);
                    Cmpx sum_over_lm;
                    for(int l=0; l<=P; l++) {
                        Cmpx sum_over_m;
                        for(int m=-l; m<=l; m++) {
                            Cmpx sph = spherical_harmonic(l, m, r.colatitude(), r.azimuth());
                            sph *= (1.0*factorial(l-fabsf(m))) / factorial(l+fabsf(m));
                            sph *= mpc_gmem[l*l+l+m];
                            sum_over_m += sph;
                        }
                        sum_over_m *= powf(r.magnitude(), -(l+1));
                        sum_over_lm += sum_over_m;
                    }
                    potential_gmem[zp*ydim*xdim + y*xdim + x] += sum_over_lm.get_re();
                }
            }
        } // thread loop
    } // if thread within limits
}




// wrapper function for FMM kernel
// =======================================================
int fmm_gpu(        const Box *const n_ptr,
                    const Cmpx *const mpc,
                    fptype *potential_gmem,
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

    if(first_time) {
        // select device to use
        // cudaSetDevice(1);
        // allocate memory on device
        cudaMalloc((void**)&n_gmem, sizeof(Box));         checkCUDAError("Allocate n_gmem");
        cudaMalloc((void**)&ni_gmem, 27 * sizeof(Box));
        cudaMalloc((void**)&mpc_gmem, (P+1)*(P+1) * sizeof(Cmpx));
        if(n_gmem == NULL || ni_gmem == NULL || mpc_gmem == NULL) {
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

    // int currentDevice;
    // cudaGetDevice(&currentDevice);
    // printf("using device %d\n", currentDevice);

    // set up kernel parameters
    #define MAXTHREADSPERBLOCK    1024
    int width = (int)powf(2, limit - n_ptr->level);
    int problem_size = width * width;
    dim3 grid = 27;
    dim3 threads = (problem_size <= 256) ? problem_size : 256;
    assert(threads.x <= MAXTHREADSPERBLOCK);    // max_threads_per_block

    if(first_time)
        printf("launching kernel with %u blocks and %u threads...\n",
                    grid.x*grid.y*grid.z, threads.x*threads.y*threads.z);

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

    first_time = 0;
    return EXIT_SUCCESS;
}





















// **************************************************
// attempt to pin 28 boxes and mpc
// **************************************************

// #include <stdio.h>
// #include <stdlib.h>
// #include <assert.h>
// #include <time.h>
// #include <sys/time.h>
// #include <cmath>
// // #include <cutil_inline.h>
// #include "mydefs.hpp"
// #include "Box.hpp"
// // #include "Queue.hpp"
// #include "Cmpx.hpp"
// #include "Vector3.hpp"
// #include "numerics.hpp"
// #include "helper_functions.hpp"

// #include "potential_calc.cu"
// #include "fmm_calc.cu"


// // FMM interaction evaluation kernel
// // ===================================
// __global__
// void fmm_kernel
    // (   const Box *const n_gmem,
        // const Box *const ni_gmem,
        // const Cmpx *const mpc_gmem,
        // fptype *potential_gmem,
        // const unsigned int limit,
        // const int P,
        // const int xdim,
        // const int ydim,
        // const int zdim,
        // const int zc    )
// {
    // int bi = blockIdx.x;
    // int ti = threadIdx.x;

    // __shared__ Box n;
    // __shared__ Box ni;
    // __shared__ int x1;
    // __shared__ int y1;
    // __shared__ int width;
    // __shared__ int stride;

    // if(ti == 0) {
        // n = *n_gmem;
    // }
    // __syncthreads();

    // if(n.interaction[bi] == NULL)
        // return;

    // if(ti == 0) {
        // // block index as interaction loop
        // ni = ni_gmem[bi];
        // width = powf(2, limit - n.level);
        // x1 = ceilf(ni.cx - width/2);
        // y1 = ceilf(ni.cy - width/2);
        // stride = ceilf(width*width / (fptype)blockDim.x);
    // }
    // __syncthreads();

    // // int ti1 = ti * stride;
    // if(ti * stride < width*width) {
        // for(int i = ti * stride; i < ti * stride + stride; i++) {
            // if(i < width*width) {
                // int x = x1 + i % width;
                // int y = y1 + (i - (x - x1)) / width;
                // for (int zp = 0; zp < zdim; zp++) { // for each potential layer in zdim
                    // Vector3 r(x-n.cx, y-n.cy, zp-zc);
                    // Cmpx sum_over_lm;
                    // for(int l=0; l<=P; l++) {
                        // Cmpx sum_over_m;
                        // for(int m=-l; m<=l; m++) {
                            // Cmpx sph = spherical_harmonic(l, m, r.colatitude(), r.azimuth());
                            // sph *= (1.0*factorial(l-fabsf(m))) / factorial(l+fabsf(m));
                            // sph *= mpc_gmem[l*l+l+m];
                            // sum_over_m += sph;
                        // }
                        // sum_over_m *= powf(r.magnitude(), -(l+1));
                        // sum_over_lm += sum_over_m;
                    // }
                    // potential_gmem[zp*ydim*xdim + y*xdim + x] += sum_over_lm.get_re();
                // }
            // }
        // } // thread loop
    // } // if thread within limits
// }




// // wrapper function for FMM kernel
// // =======================================================
// int fmm_gpu(        const Box *const n_ptr,
                    // const Cmpx *const mpc,
                    // fptype *potential_gmem,
                    // const unsigned int limit,
                    // const int P,    // multipole series truncation (l = 0...P)
                    // const int xdim, const int ydim, const int zdim,
                    // const int zc,   // charge layer
                    // const int use_gpu,
                    // const int verbose_level
                // )
// {
    // int status = 0;
    // static int first_time = 1;

    // // allocate pinned memory for nbox, 27 interaction boxes and mpc
    // static void *cpumem_pinned = NULL;
    // static void *gpumem_pinned = NULL;
    // if(first_time) {
        // cudaHostAlloc((void **)&cpumem_pinned, 28*sizeof(Box) + (P+1)*(P+1)*sizeof(Cmpx), cudaHostAllocMapped);
        // checkCUDAError("cudaHostAllocMapped");
        // // Get the device pointers to the mapped memory
        // cudaHostGetDevicePointer((void **)&gpumem_pinned, (void *)cpumem_pinned, 0);
        // checkCUDAError("cudaHostGetDevicePointer");
    // }

    // // concatenate nbox, 27 interaction boxes and mpc in CPU pinned memory
    // memcpy(cpumem_pinned, n_ptr, sizeof(Box));
    // for(int bi = 0; bi < 27; bi++) {
        // if(n_ptr->interaction[bi] != NULL)
            // memcpy(((char*)cpumem_pinned) + (1+bi)*sizeof(Box), n_ptr->interaction[bi], sizeof(Box));
    // }
    // memcpy(((char*)cpumem_pinned) + 28*sizeof(Box), mpc, sizeof(Cmpx));

    // // set up device memory pointers
    // static Box *n_gmem = (Box*)gpumem_pinned;
    // static Box *ni_gmem = (Box*)(((char*)gpumem_pinned) + sizeof(Box));
    // static Cmpx *mpc_gmem = (Cmpx*)(((char*)gpumem_pinned) + 28*sizeof(Box));

    // // // allocate memory on device
    // // if(first_time) {
        // // cudaMalloc((void**)&n_gmem, sizeof(Box));         checkCUDAError("Allocate n_gmem");
        // // cudaMalloc((void**)&ni_gmem, 27 * sizeof(Box));
        // // cudaMalloc((void**)&mpc_gmem, (P+1)*(P+1) * sizeof(Cmpx));
        // // if(n_gmem == NULL || ni_gmem == NULL || mpc_gmem == NULL) {
            // // fprintf(stderr, "%s:%d Error allocating memory on GPU\n", __FILE__, __LINE__);
            // // return EXIT_FAILURE;
        // // }
    // // }

    // // copy required data to GPU global memory
    // // cudaMemcpy(n_gmem, n_ptr, sizeof(Box), cudaMemcpyHostToDevice);
    // // checkCUDAError("Copying n_gmem");
    // // cudaMemcpy(ni_gmem, ni_cpumem, 27 * sizeof(Box), cudaMemcpyHostToDevice);
    // // checkCUDAError("Copying ni_cpumem");
    // // cudaMemcpy(mpc_gmem, mpc, (P+1)*(P+1) * sizeof(Cmpx), cudaMemcpyHostToDevice);
    // // checkCUDAError("Copying mpc_gmem");

    // // int currentDevice;
    // // cudaGetDevice(&currentDevice);
    // // printf("using device %d\n", currentDevice);

    // // set up kernel parameters
    // #define MAXTHREADSPERBLOCK    1024
    // int width = (int)powf(2, limit - n_ptr->level);
    // int problem_size = width * width;
    // dim3 grid = 27;
    // dim3 threads = (problem_size <= 256) ? problem_size : 256;
    // assert(threads.x <= MAXTHREADSPERBLOCK);    // max_threads_per_block

    // if(first_time)
        // printf("launching kernel with %u blocks and %u threads...\n",
                    // grid.x*grid.y*grid.z, threads.x*threads.y*threads.z);

    // // start timer
    // timeval time1, time2;
    // status |= gettimeofday(&time1, NULL);

    // if(first_time) {
        // // select device to use
        // // cudaSetDevice(1);
    // }

    // // launch the kernel
    // fmm_kernel <<<grid, threads>>>
        // (n_gmem, ni_gmem, mpc_gmem, potential_gmem, limit, P, xdim, ydim, zdim, zc);
    // checkCUDAError("Exeuting Kernel fmm_kernel()");
    // cudaThreadSynchronize();

    // // read the timer
    // status |= gettimeofday(&time2, NULL);
    // // double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    // // printf("  Kernel completed in %f seconds.\n", deltatime);

    // first_time = 0;
    // return EXIT_SUCCESS;
// }
