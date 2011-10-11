#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include "Vector3.hpp"
#include "mydefs.hpp"

// Print a message if a CUDA error occurred
void checkCUDAError(const char *msg) {
    cudaError_t err = cudaGetLastError();
    if(cudaSuccess != err) {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }
}

// Kernel definition (2nd version)
__global__ void
calc_potential_exact_kernel(    const fptype *charge_gmem,
                                const int xdim, const int ydim, const int zdim,
                                const fptype dx, const fptype dy, const fptype dz,
                                fptype *potential_gmem  )
{
    const int N = xdim*ydim*zdim;
    int bstride = ceilf(N / (fptype)gridDim.x);
    int tstride = ceilf(N / (fptype)blockDim.x);
    int bi1 = blockIdx.x * bstride;
    int ti1 = threadIdx.x * tstride;

    if(bi1 >= N) // if start of block exceeds, don't proceed
        return;

    for(int bi = bi1; bi < bi1 + bstride; bi++) // block subindex
    {
        if(bi >= N)  // consider this target point only if this point doesn't exceed
            break;

        // reset shared memory
        extern __shared__ fptype sdata[];
        sdata[threadIdx.x] = 0;

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
            for(int i = ti1; i < ti1 + tstride; i++)
            {
                if(i < N) { // calculate potential only if this point doesn't exceed
                    fptype dV = dx*dy*dz;
                    fptype q = charge_gmem[i] * dV;
                    // __syncthreads();
                    // if(q != 0)
                    {
                        int x_ = i % xdim; // source coords per thread
                        int y_ = ((i - x_) / xdim) % ydim;
                        int z_ = (i - x_ - y_*xdim) / (xdim*ydim);
                        // potential due this thread's charge
                        if(bi != i) { // skip on itself to avoid div by zero
                            Vector3 R((x-x_)*dx, (y-y_)*dy, (z-z_)*dz);
                            pot += q / R.magnitude();
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
    } // block subindex loop
}


// Exact O(N^2) calculation of potential
int calc_potential_exact_gpu( const fptype *charge,
                        const int xdim, const int ydim, const int zdim,
                        const fptype dx, const fptype dy, const fptype dz,
                        fptype *potential)
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

// set up device memory pointers
    static fptype *charge_gmem = NULL;
    static fptype *potential_gmem = NULL;

// allocate memory on device
    if(first_time) {
        cudaMalloc((void**)&charge_gmem,    zdim*ydim*xdim * sizeof(fptype));
        checkCUDAError("Allocate charge_gmem");
        cudaMalloc((void**)&potential_gmem, zdim*ydim*xdim * sizeof(fptype));
        checkCUDAError("Allocate potential_gmem");
        if(charge_gmem == NULL || potential_gmem == NULL) {
            fprintf(stderr, "%s:%d Error allocating memory on GPU\n", __FILE__, __LINE__);
            return EXIT_FAILURE;
        }
    }

// copy charge array to device global memory
    cudaMemcpy(charge_gmem, charge, zdim*ydim*xdim * sizeof(fptype), cudaMemcpyHostToDevice);
    checkCUDAError("Copying charge_gmem");


// set up kernel parameters
    int problem_size = zdim*ydim*xdim;
    dim3 grid = ceil(problem_size / ceil((fptype)problem_size/MAXBLOCKS));
    // dim3 grid = (problem_size < MAXBLOCKS) ? problem_size : MAXBLOCKS;
    dim3 threads = 1024;
    if(first_time) {
        assert(grid.x <= MAXBLOCKS);
        assert(threads.x <= MAXTHREADSPERBLOCK);
        printf("launching kernel with %u blocks and %u threads...\n",
                    grid.x*grid.y*grid.z, threads.x*threads.y*threads.z);
    }

// launch the kernel
    calc_potential_exact_kernel <<<grid, threads, 1024 * sizeof(fptype)>>>
        (charge_gmem, xdim, ydim, zdim, dx, dy, dz, potential_gmem);
    checkCUDAError("Exeuting Kernel calc_potential_exact_kernel()");
    // cudaThreadSynchronize();

// copy potential (the result of kernel) to host main memory
    cudaMemcpy(potential, potential_gmem, zdim*ydim*xdim * sizeof(fptype), cudaMemcpyDeviceToHost);
    checkCUDAError("Copying potential_gmem");
    // cudaThreadSynchronize();

    // cudaFree(charge_gmem);
    // cudaFree(potential_gmem);

    first_time = 0;
    return EXIT_SUCCESS;
}
