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
    (   const Cmpx n_cxy,
        const int width,
        const int num_interactions,
        const Cmpx *const ni_cxy_gmem,
        const Cmpx *const mpc_gmem,
        fptype *potential_gmem,
        const int P,
        const int xdim, const int ydim, const int zdim,
        const fptype dx, const fptype dy, const fptype dz,
        const int zc    )
{
    int bi = blockIdx.x;
    int ti = threadIdx.x;

    if(bi >= num_interactions)
        return;

    __shared__ fptype ni_cx;
    __shared__ fptype ni_cy;
    __shared__ int x1;
    __shared__ int y1;
    __shared__ int stride;

    if(ti == 0) {
        // block index as interaction loop
        ni_cx = ni_cxy_gmem[bi].get_re();
        ni_cy = ni_cxy_gmem[bi].get_im();
        x1 = ceilf(ni_cx - width/2);
        y1 = ceilf(ni_cy - width/2);
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
                    Vector3 r((x-n_cxy.get_re())*dx, (y-n_cxy.get_im())*dy, (zp-zc)*dz);
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
                    const Cmpx *const mpc_gmem,
                    fptype *potential_gmem,
                    const unsigned int limit,
                    const int P,    // multipole series truncation (l = 0...P)
                    const int xdim, const int ydim, const int zdim,
                    const fptype dx, const fptype dy, const fptype dz,
                    const int zc,   // charge layer
                    const int use_gpu,
                    const int verbosity
                )
{
    // int status = 0;
    static int first_time = 1;

// setup device
    if(first_time) {
        // select device to use
        cudaSetDevice(1); // device-0 is tied to display output

        int currentDevice;
        cudaGetDevice(&currentDevice);
        printf("using device %d\n", currentDevice);
    }

// concatenate coordinates of 27 interaction boxes in CPU pinned memory
    static Cmpx *ni_cxy = NULL;      // pinned memory on CPU
    static Cmpx *ni_cxy_gmem = NULL; // mapped pointer on device memory
    if(first_time) {
        cudaHostAlloc((void**)&ni_cxy, 27*sizeof(Cmpx), cudaHostAllocMapped);
        checkCUDAError("cudaHostAllocMapped");
        // Get the device pointer to the mapped memory
        cudaHostGetDevicePointer((void**)&ni_cxy_gmem, (void*)ni_cxy, 0);
        checkCUDAError("cudaHostGetDevicePointer");
    }
    int num_interactions = 0;
    for(int bi = 0; bi < 27; bi++) {
        if(n_ptr->interaction[bi] != NULL) {
            ni_cxy[num_interactions] = Cmpx(n_ptr->interaction[bi]->cx, n_ptr->interaction[bi]->cy);
            num_interactions++;
        }
    }

// prepare other arguments
    Cmpx n_cxy(n_ptr->cx, n_ptr->cy);
    int width = (int)powf(2, limit - n_ptr->level);

// set up kernel parameters
    int problem_size = width * width;
    dim3 grid = 27;
    dim3 threads = (problem_size <= 256) ? problem_size : 256;
    if(first_time) {
        assert(grid.x <= MAXBLOCKS);
        assert(threads.x <= MAXTHREADSPERBLOCK);
        printf("launching kernel with %u blocks and %u threads...\n",
                    grid.x*grid.y*grid.z, threads.x*threads.y*threads.z);
    }

// launch the kernel
    fmm_kernel <<<grid, threads>>>
        (n_cxy, width, num_interactions, ni_cxy_gmem, mpc_gmem, potential_gmem, P, xdim, ydim, zdim, dx, dy, dz, zc);
    checkCUDAError("Exeuting Kernel fmm_kernel()");
    cudaThreadSynchronize();


    first_time = 0;
    return EXIT_SUCCESS;
}
