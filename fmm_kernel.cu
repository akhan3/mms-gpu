#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include <cutil_inline.h>
#include "mydefs.hpp"
#include "Box.hpp"
#include "Queue.hpp"
#include "Cmpx.hpp"
#include "Vector3.hpp"
#include "numerics.hpp"
#include "helper_functions.hpp"

#include "potential_calc.cu"

// #define SDATA(index)      cutilBankChecker(sdata, index)
// #define SDATA(index)      sdata[index]


// FMM algorithm in BFS
// ===============================
__global__
void fmm_kernel(    const fptype *charge_gmem,
                    fptype *potential_gmem,
                    void **queue_gmem,
                    const Box *root,
                    const unsigned int limit,
                    const int P,    // multipole series truncation (l = 0...P)
                    const int xdim, const int ydim, const int zdim,
                    const int zc   // charge layer
                )
{
    int ti = threadIdx.x;

    if(threadIdx.x == 0)
    {
        const unsigned int N = (unsigned int)powf(4, limit);
        Queue Q_tree(N, queue_gmem);
        Q_tree.enqueue((void*)root);

    // iterate over all the boxes in tree
        while(!Q_tree.isEmpty()) {
            Box *n = (Box*)Q_tree.dequeue();
            // populate queue with children nodes
            if(n->level < limit)
                for(int i=0; i<=3; i++)
                    Q_tree.enqueue(n->child[i]);

            if(n->level <= 1)   // no FMM steps for Level-0 and Level-1
                continue;

    // function to perform on node
            if(n->is_pruned()) {
                continue;
            }

        // Calculate multipole coefficients for the source box
            fptype q = 0;
            Cmpx multipole_coeff[3+1][2*3+1];
            // checking for source charges in the source box
            fptype charge_found = 0;
            fptype width = powf(2, limit-n->level);
            int yy1 = ceil(n->cy-width/2);
            int yy2 = floor(n->cy+width/2);
            for(int yy=yy1; yy<=yy2; yy++) {
            // for(int yy=ceil(n->cy-width/2); yy<=floor(n->cy+width/2); yy++) {
                for(int xx=ceil(n->cx-width/2); xx<=floor(n->cx+width/2); xx++) {
                    q = charge_gmem[yy*xdim + xx];
                    if(q != 0) { // if charge found
                        charge_found = 1;
                        Cmpx r_(xx - n->cx, yy - n->cy);
                        for(int l=0; l<=P; l++) {
                            for(int m=-l; m<=l; m++) {
                                Cmpx sph = spherical_harmonic(l, m, M_PI/2, r_.get_ang()).conjugate();
                                sph *= q * pow(r_.get_mag(), l);
                                multipole_coeff[l][m+l] += sph;
                                // multipole_coeff[l][m+l] += q * pow(r_.get_mag(), l) * spherical_harmonic(l, m, M_PI/2, r_.get_ang()).conjugate();
                            } // m loop
                        } // l loop
                    } // if(q != 0)
                } // source charge loop
            } // source charge loop
            // NEWLINE;


            if(! charge_found) {
                n->prune();
                continue;
            }

            if(charge_found)
            {
                for (int zp = 0; zp < zdim; zp++) // for each potential layer in zdim
                {
                // calculation of potential at the boxes in interaction list
                    for(int i=0; i<27; i++) {
                        Box *ni = n->interaction[i];
                        if(ni != NULL) {
                            for(int yy=ceil(ni->cy-width/2); yy<=floor(ni->cy+width/2); yy++) {
                                for(int xx=ceil(ni->cx-width/2); xx<=floor(ni->cx+width/2); xx++) {
                                    Vector3 r(xx - n->cx, yy - n->cy, zp - zc);
                                    Cmpx sum_over_lm;
                                    for(int l=0; l<=P; l++) {
                                        Cmpx sum_over_m;
                                        for(int m=-l; m<=l; m++) {
                                            Cmpx sph = spherical_harmonic(l, m, r.colatitude(), r.azimuth());
                                            sph *= (1.0*factorial(l-abs(m))) / factorial(l+abs(m));
                                            sph *= multipole_coeff[l][m+l];
                                            sum_over_m += sph;
                                            // sum_over_m += (1.0*factorial(l-abs(m))) / factorial(l+abs(m)) * multipole_coeff[l][m+l] * spherical_harmonic(l, m, r.colatitude(), r.azimuth());
                                        }
                                        sum_over_m *= 1 / pow(r.magnitude(), l+1);
                                        sum_over_lm += sum_over_m;
                                        // sum_over_lm += 1 / pow(r.magnitude(), l+1) * sum_over_m;
                                    }
                                    potential_gmem[zp*ydim*xdim + yy*xdim + xx] += sum_over_lm.get_re();
                                    // potential[yy*xdim+xx] += (sum_over_lm.get_re() > 0) ? sum_over_lm.get_mag() : -sum_over_lm.get_mag();

                                    // const fptype threshold = 1e-2;
                                    // fptype modangle = fabs(sum_over_lm.get_ang());
                                    // modangle = (modangle < M_PI-modangle) ? modangle : M_PI-modangle;
                                    // if(modangle > threshold) {
                                        // if(verbose_level >= 0)
                                            // printf("PANIC!! L%d   R=%g   angle=%g\n", n->level, r.magnitude(), modangle);
                                        // fprintf(paniclog, "%d   %g   %g\n", n->level, r.magnitude(), modangle);
                                    // }
                                }
                            }
                        } // if(ni != NULL)
                    } // interaction loop

                // calculation with neighbor list at the deepest level
                    if(n->level == limit) {
                        if(zp != zc) { // neighbor on other layers at self position
                            Vector3 r(0, 0, zp - zc);
                            potential_gmem[zp*ydim*xdim + (int)(n->cy*xdim + n->cx)] += q / r.magnitude();
                        }
                        for(int i=0; i<8; i++) {
                            Box *nb = n->neighbor[i];
                            if(nb != NULL) {
                                Vector3 r(nb->cx - n->cx, nb->cy - n->cy, zp - zc);
                                potential_gmem[zp*ydim*xdim + (int)(nb->cy*xdim + nb->cx)] += q / r.magnitude();
                            }
                        } // neighbor loop
                    } // if deepest level
                } // for each potential layer in zdim
            } // if(charge_found)
        } // while(!Q_tree.isEmpty())
    } // if(threadIdx.x == 0)
    return;
}




int fmm_gpu(        const fptype *charge,
                    fptype *potential,
                    const Box *root,
                    const unsigned int limit,
                    const unsigned int actual_limit,
                    const int P,    // multipole series truncation (l = 0...P)
                    const int xdim, const int ydim, const int zdim,
                    const int zc,   // charge layer
                    FILE *paniclog,
                    const int verbose_level
                )
{
    assert(zdim == 1);
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
    dim3 grid = 1;
    const int stride = ceil(problem_size / (fptype)MAXTHREADSPERBLOCK);
    dim3 threads = 32;
    assert(threads.x <= MAXTHREADSPERBLOCK);    // max_threads_per_block

    // if(first_time) {
        // printf("x=%u, y=%u, z=%u, threads.x=%u, threads.y=%u, threads.z=%u, stride=%u, grid.x=%u, grid.y=%u, grid.z=%u\n",
                // xdim, ydim, zdim, threads.x,    threads.y,    threads.z,    stride,    grid.x,    grid.y,    grid.z);
        // printf("launching kernel with %u blocks and %u threads...\n",
                    // grid.x*grid.y*grid.z, threads.x*threads.y*threads.z);
    // }

// allocate memory for FMM Queue
    void **queue_mem = (void**)malloc(xdim*ydim * sizeof(void*));
    if(queue_mem == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    // Queue Q_tree(xdim*ydim, queue_mem);

    // start timer
    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);

    // launch the kernel
    fmm_kernel <<<grid, threads, 1024 * sizeof(fptype)>>>
        (charge_d, potential_d, queue_mem, root,
                    limit, P, xdim, ydim, zdim, zc);

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
