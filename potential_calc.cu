#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <cutil_inline.h>
#include "Vector3.hpp"
#include "mydefs.hpp"

// #define SDATA(index)      cutilBankChecker(sdata, index)
// #define SDATA(index)      sdata[index]


// Kernel definition (2nd version)
__global__ void
calc_potential_exact_kernel(   const fptype *charge_gmem,
                                const int xdim, const int ydim, const int zdim,
                                const int stride,
                                fptype *potential_gmem  )
{
    const int N = xdim*ydim*zdim;
    // const int n = blockIdx.x * blockDim.x + threadIdx.x;
    int bi = blockIdx.x;
    int ti1 = threadIdx.x * stride;

    if(bi >= N) // if block exceeds, don't proceed (unnecessary check)
        return;

// reset shared memory
    extern __shared__ fptype sdata[];
    sdata[threadIdx.x] = 0;
    // __syncthreads;

    if(ti1 < N) // calculate potential only if start of thread doesn't exceed
    {
        __shared__ int xs, ys, zs; // target coords per block
        if(threadIdx.x == 0) {
            xs = bi % xdim;
            ys = ((bi - xs) / xdim) % ydim;
            zs = (bi - xs - ys*xdim) / (xdim*ydim);
        }
        __syncthreads();
        int x = xs;
        int y = ys;
        int z = zs;

        fptype pot = 0;
        // int i = ti1;
        for(int i = ti1; i < ti1 + stride; i++)
        {
            if(i < N) { // calculate potential only if this point doesn't exceed
                fptype q = charge_gmem[i];
                // __syncthreads();
                // if(q != 0)
                {
                    int x_ = i % xdim; // source coords per thread
                    int y_ = ((i - x_) / xdim) % ydim;
                    int z_ = (i - x_ - y_*xdim) / (xdim*ydim);
                    // potential due this thread's charge
                    if(bi != i) { // skip on itself to avoid div by zero
                        // Vector3 V(x-x_, y-y_, z-z_);
                        // fptype dist = V.magnitude();
                        fptype R = sqrtf( (x-x_)*(x-x_) + (y-y_)*(y-y_) + (z-z_)*(z-z_) );
                        pot += q / R;
                    }
                }
            }
        }
        sdata[threadIdx.x] = pot;
    }
    __syncthreads();

// parallel reduction to sum up potential from threads
// (must use all 1024 threads)
    if (threadIdx.x < 512)
        sdata[threadIdx.x] += sdata[threadIdx.x + 512];
    __syncthreads();
    if (threadIdx.x < 256)
        sdata[threadIdx.x] += sdata[threadIdx.x + 256];
    __syncthreads();
    if (threadIdx.x < 128)
        sdata[threadIdx.x] += sdata[threadIdx.x + 128];
    __syncthreads();
    if (threadIdx.x < 64)
        sdata[threadIdx.x] += sdata[threadIdx.x + 64];
    __syncthreads();
    if (threadIdx.x < 32)
        sdata[threadIdx.x] += sdata[threadIdx.x + 32];
    __syncthreads();
    if (threadIdx.x < 16)
        sdata[threadIdx.x] += sdata[threadIdx.x + 16];
    __syncthreads();
    if (threadIdx.x < 8)
        sdata[threadIdx.x] += sdata[threadIdx.x + 8];
    __syncthreads();
    if (threadIdx.x < 4)
        sdata[threadIdx.x] += sdata[threadIdx.x + 4];
    __syncthreads();
    if (threadIdx.x < 2)
        sdata[threadIdx.x] += sdata[threadIdx.x + 2];
    __syncthreads();
    if (threadIdx.x < 1)
        sdata[threadIdx.x] += sdata[threadIdx.x + 1];
    __syncthreads();
// write summed potential to global memory
    if(threadIdx.x == 0)
        potential_gmem[bi] = sdata[0];
}



// Exact O(N^2) calculation of potential
int calc_potential_exact_gpu( const fptype *charge,
                        const int xdim, const int ydim, const int zdim,
                        fptype *potential)
{

    int status = 0;
    static int first_time = 1;

    // set up device memory pointers
    static fptype *charge_d = NULL;
    static fptype *potential_d = NULL;

    if(first_time) {
        // select device to use
        cudaSetDevice(1);
        // allocate memory on device
        cutilSafeCall( cudaMalloc( (void**)&charge_d,    zdim*ydim*xdim * sizeof(fptype) ) );
        cutilSafeCall( cudaMalloc( (void**)&potential_d, zdim*ydim*xdim * sizeof(fptype) ) );
        if(charge_d == NULL || potential_d == NULL) {
            fprintf(stderr, "%s:%d Error allocating memory on GPU\n", __FILE__, __LINE__);
            return EXIT_FAILURE;
        }
    }

    // int currentDevice;
    // cudaGetDevice(&currentDevice);
    // printf("using device %d\n", currentDevice);

    // copy charge array to device global memory
    cutilSafeCall( cudaMemcpy( charge_d, charge, zdim*ydim*xdim * sizeof(fptype), cudaMemcpyHostToDevice ) );


    // set up kernel parameters
    #ifdef __DEVICE_EMULATION__
        #define MAXTHREADSPERBLOCK    64
    #else
        #define MAXTHREADSPERBLOCK    1024
    #endif
    int problem_size = zdim*ydim*xdim;
    // dim3 grid = ceil(total_threads / (fptype)MAXTHREADSPERBLOCK);
    // dim3 threads(MAXTHREADSPERBLOCK, 1, 1);
    dim3 grid = problem_size;
    const int stride = ceil(problem_size / (fptype)MAXTHREADSPERBLOCK);
    dim3 threads(MAXTHREADSPERBLOCK, 1, 1);
    assert(threads.x <= MAXTHREADSPERBLOCK);    // max_threads_per_block

    // if(first_time) {
        // printf("x=%u, y=%u, z=%u, threads.x=%u, threads.y=%u, threads.z=%u, stride=%u, grid.x=%u, grid.y=%u, grid.z=%u\n",
                // xdim, ydim, zdim, threads.x,    threads.y,    threads.z,    stride,    grid.x,    grid.y,    grid.z);
        // printf("launching kernel with %u blocks and %u threads...\n",
                    // grid.x*grid.y*grid.z, threads.x*threads.y*threads.z);
    // }

    // start timer
    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);

    // launch the kernel
    calc_potential_exact_kernel <<<grid, threads, 1024 * sizeof(fptype)>>>
        (charge_d, xdim, ydim, zdim, stride, potential_d);
    cutilCheckMsg("Kernel execution failed");
    cudaThreadSynchronize();

    // read the timer
    status |= gettimeofday(&time2, NULL);
    // double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    // printf("  Kernel completed in %f seconds.\n", deltatime);

    // copy potential (the result of kernel) to host main memory
    cutilSafeCall( cudaMemcpy( potential, potential_d, zdim*ydim*xdim * sizeof(fptype), cudaMemcpyDeviceToHost ) );
    // cutilSafeCall( cudaFree(charge_d) );
    // cutilSafeCall( cudaFree(potential_d) );
    cudaThreadSynchronize();

    first_time = 0;
    return EXIT_SUCCESS;
}
