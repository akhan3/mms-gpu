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

// #define SDATA(index)      cutilBankChecker(sdata, index)
// #define SDATA(index)      sdata[index]


__global__
void fmm_kernel0(const Cmpx sph, const Cmpx mpc, int l, int m, Vector3 r, Cmpx *ans_d)
{
    if(threadIdx.x == 0) {
        Cmpx product(sph);
        product *= (1.0*factorial(l-fabsf(m))) / factorial(l+fabsf(m));
        product *= mpc;
        product *= powf(r.magnitude(), -(l+1));
        *ans_d = product;
    }
}

__global__
void fmm_kernel1(const Cmpx mpc, int l, int m, Vector3 r, Cmpx *ans_d)
{
    if(threadIdx.x == 0) {
        Cmpx product = spherical_harmonic(l, m, r.colatitude(), r.azimuth());
        product *= (1.0*factorial(l-fabsf(m))) / factorial(l+fabsf(m));
        product *= mpc;
        product *= powf(r.magnitude(), -(l+1));
        *ans_d = product;
    }
}

__global__
void fmm_kernel2(Cmpx *mpc_gmem, int l, int m, int P, Vector3 r, Cmpx *ans_d)
{
    if(threadIdx.x == 0) {
        Cmpx product = spherical_harmonic(l, m, r.colatitude(), r.azimuth());
        product *= (1.0*factorial(l-fabsf(m))) / factorial(l+fabsf(m));
        product *= mpc_gmem[l*(2*P+1) + m+l];
        product *= powf(r.magnitude(), -(l+1));
        *ans_d = product;
    }
}


__global__
void fmm_kernel22(  const Cmpx *const multipole_coeff_gmem,
                    fptype *potential_gmem,
                    const unsigned int limit,
                    const fptype width,
                    const int P,    // multipole series truncation (l = 0...P)
                    const fptype n_cx, const fptype n_cy,
                    const int xx, const int yy,
                    const int xdim,
                    const int ydim,
                    const int zdim,
                    const int zp, const int zc
                )
{
    if(threadIdx.x == 0) {
                        Vector3 r(xx-n_cx, yy-n_cy, zp-zc);
                        Cmpx sum_over_lm;
                        for(int l=0; l<=P; l++) {
                            Cmpx sum_over_m;
                            for(int m=-l; m<=l; m++) {
                                Cmpx sph = spherical_harmonic(l, m, r.colatitude(), r.azimuth());
                                sph *= (1.0*factorial(l-fabsf(m))) / factorial(l+fabsf(m));
                                sph *= multipole_coeff_gmem[l*(2*P+1) + m+l];
                                sum_over_m += sph;
                            }
                            sum_over_m *= 1 / powf(r.magnitude(), l+1);
                            sum_over_lm += sum_over_m;
                        }
                        // atomicAdd(&potential_gmem[zp*ydim*xdim + yy*xdim + xx], sum_over_lm.get_re());
                        potential_gmem[zp*ydim*xdim + yy*xdim + xx] += sum_over_lm.get_re();
    } // if thread 0
}


// FMM algorithm in BFS
// ===============================
__global__
void fmm_kernel(    const Box *const n_gmem,
                    const Box *const ni_gmem, const int nbi,
                    const Cmpx *const multipole_coeff_gmem,
                    fptype *potential_gmem,
                    const unsigned int limit,
                    const int P,    // multipole series truncation (l = 0...P)
                    const int xdim,
                    const int ydim,
                    const int zdim,
                    const int zc   // charge layer
                )
{
    // int bi = blockIdx.x;
    int ti = threadIdx.x;
    if(ti == 0) {
        const Box n = *n_gmem;
        fptype width = powf(2, limit - n.level);
        // block index as interaction loop
        if(n.interaction[nbi] != NULL) {
            Box ni = ni_gmem[nbi];
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
                                sph *= multipole_coeff_gmem[l*(2*P+1) + m+l];
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


// Move the vector to the host and sum
float sumOnHost(const float *v_d, int n)
{
    float sum = 0.f;
    int i;
    // create space on the host for the device data
    float *v_h = (float*)malloc(n*sizeof(float));
    assert(v_h != NULL); // check if the malloc succeeded
    // copy the vector from the device to the host
    cudaMemcpy(v_h, v_d, n*sizeof(float), cudaMemcpyDeviceToHost);
    for(i=0; i<n; i++) sum += v_h[i];
    free(v_h); // free the vector on host
    return(sum);
}


int fmm_gpu(        const Box *const n_ptr,
                    const Cmpx *const multipole_coeff,
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
    static Cmpx *multipole_coeff_gmem = NULL;
    static fptype *potential_gmem = NULL;

    if(first_time) {
        // select device to use
        // cudaSetDevice(1);
        // allocate memory on device
        cudaMalloc((void**)&n_gmem, sizeof(Box));         checkCUDAError("Allocate n_gmem");
        cudaMalloc((void**)&ni_gmem, 27 * sizeof(Box));
        cudaMalloc((void**)&multipole_coeff_gmem, (P+1)*(2*P+1) * sizeof(Cmpx));
        cudaMalloc((void**)&potential_gmem, zdim*ydim*xdim * sizeof(fptype));
        if(n_gmem == NULL || ni_gmem == NULL || multipole_coeff_gmem == NULL || potential_gmem == NULL) {
            fprintf(stderr, "%s:%d Error allocating memory on GPU\n", __FILE__, __LINE__);
            return EXIT_FAILURE;
        }
    }

    // copy required data to GPU global memory
    cudaMemcpy(n_gmem, n_ptr, sizeof(Box), cudaMemcpyHostToDevice);
    checkCUDAError("Copying n_gmem");
    cudaMemcpy(ni_gmem, ni_cpumem, 27 * sizeof(Box), cudaMemcpyHostToDevice);
    checkCUDAError("Copying ni_cpumem");
    cudaMemcpy(multipole_coeff_gmem, multipole_coeff, (P+1)*(2*P+1) * sizeof(Cmpx), cudaMemcpyHostToDevice);
    checkCUDAError("Copying multipole_coeff_gmem");
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

    // if(0)
    if(use_gpu == 2 || use_gpu == 3) {
        if(first_time)
            printf("\n\nEmulating GPU...\n\n");
        first_time = 0;

        const Box n = *n_ptr;
        fptype width = powf(2, limit - n.level);
        for(int bi = 0; bi < 27; bi++) {
            if(n.interaction[bi] != NULL) {
                Box ni = *(n.interaction[bi]);
                for(int yy = ceilf(ni.cy-width/2); yy <= floorf(ni.cy+width/2); yy++) {
                    for(int xx = ceilf(ni.cx-width/2); xx <= floorf(ni.cx+width/2); xx++) {
                        for (int zp = 0; zp < zdim; zp++) { // for each potential layer in zdim
                            Vector3 r(xx-n.cx, yy-n.cy, zp-zc);
                            Cmpx sum_over_lm;
                            for(int l=0; l<=P; l++) {
                                for(int m=-l; m<=l; m++) {
                                    Cmpx sph;
                                    if(use_gpu != 3) {
                                        Cmpx mpc = multipole_coeff[l*(2*P+1) + m+l];
                                        sph = spherical_harmonic(l, m, r.colatitude(), r.azimuth());
                                        sph *= (1.0*factorial(l-fabsf(m))) / factorial(l+fabsf(m));
                                        sph *= mpc;
                                        sph *= powf(r.magnitude(), -(l+1));
                                    }
                                    else {
                                        Cmpx *ans_d = NULL;
                                        cudaMalloc((void**)&ans_d, sizeof(Cmpx)); checkCUDAError("Allocate ans_d");
                                        // fmm_kernel0 <<<1, 1>>> (sph, mpc, l, m, r, ans_d);
                                        // fmm_kernel1 <<<1, 1>>> (mpc, l, m, r, ans_d);
                                        fmm_kernel2 <<<1, 1>>> (multipole_coeff_gmem, l, m, P, r, ans_d);
                                        cudaMemcpy(&sph, ans_d, sizeof(Cmpx), cudaMemcpyDeviceToHost);
                                        checkCUDAError("Copying summation");
                                    }
                                    sum_over_lm += sph;
                                }
                            }
                            potential[zp*ydim*xdim + yy*xdim + xx] += sum_over_lm.get_re();
                        }
                    }
                }
            } // if(ni != NULL)
        } // interaction loop
        return 0;
    }



// debug MPC
    Cmpx *mpc = (Cmpx*)malloc((P+1)*(2*P+1) * sizeof(Cmpx));
    cudaMemcpy(mpc, multipole_coeff_gmem, (P+1)*(2*P+1) * sizeof(Cmpx), cudaMemcpyDeviceToHost);
    checkCUDAError("Copying multipole_coeff_gmem");
    printf("Box(%.1f,%.1f):\n", n_ptr->cx, n_ptr->cy);
    for(int i = 0; i < (P+1)*(2*P+1); i++) {
        printf("MPC[%d] = %s = %s\n", i, multipole_coeff[i].cartesian(), mpc[i].cartesian());
    }
    printf("\n");
    return 0;


    // launch the kernel
    // fmm_kernel <<<grid, threads>>>
        // (n_gmem, ni_gmem, multipole_coeff_gmem, potential_gmem, limit, P, xdim, ydim, zdim, zc);

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
