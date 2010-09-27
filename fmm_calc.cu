#ifdef _OPENMP
#include <omp.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include "Box.hpp"
#include "Queue.hpp"
#include "Cmpx.hpp"
#include "Vector3.hpp"
#include "numerics.hpp"
#include "helper_functions.hpp"


// FMM algorithm in BFS
// ===============================
int fmm_bfs(        const fptype *charge,
                    fptype *potential,
                    fptype *potential_gmem,
                    Box *const root,
                    const unsigned int limit,
                    const unsigned int actual_limit,
                    const int P,    // multipole series truncation (l = 0...P)
                    const int xdim, const int ydim, const int zdim,
                    const int zc,   // charge layer
                    FILE *paniclog,
                    const int use_gpu,
                    const int verbose_level
                )
{
    int status = 0;
    assert(limit == actual_limit);  // unable to support arbitrary depth calculations.
    assert(limit <= actual_limit);
    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);
    if(verbose_level >= 6)
        printf("    Executing FMM algorithm...\n");
    unsigned int prev_level = 0;

    const unsigned int N = (unsigned int)pow(4, limit);
    void **queue_mem = (void**)malloc(N * sizeof(void*));
    if(queue_mem == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    Queue Q_tree(N, queue_mem);
    Q_tree.enqueue((void*)root);

// timers for profiling
    double t_coeff = 0;
    double t_potential = 0;
    double t_potential_nearest = 0;
    double deltatime = 0;
    timeval t1, t2;

// iterate over all the boxes in tree
    // int total_boxes = (4*xdim*ydim - 1) / 3;
    // for(int b = 0; b < total_boxes; b++)
    // {
        // Box *n = root + b;
    while(!Q_tree.isEmpty())
    {
        Box *n = (Box*)Q_tree.dequeue();
        if(n->level < limit)
            for(int i=0; i<=3; i++)
                Q_tree.enqueue(n->child[i]); // populate queue with children nodes

        if(n->level <= 1)   // no FMM steps for Level-0 and Level-1
            continue;

// function to perform on node
        if(prev_level != n->level) {
            if(prev_level >= 2) {
                status |= gettimeofday(&time2, NULL);
                double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
                status |= gettimeofday(&time1, NULL);
                if(verbose_level >= 10)
                    printf("done in %f seconds.\n", deltatime); fflush(NULL);
            // saving this level potential
                // char filename_pot[200];
                // if     (use_gpu == 0) sprintf(filename_pot, "potential_cpu_L%d.dat", n->level - 1);
                // else if(use_gpu == 1) sprintf(filename_pot, "potential_gpu_L%d.dat", n->level - 1);
                // else if(use_gpu == 2) sprintf(filename_pot, "potential_gpuemu_L%d.dat", n->level - 1);
                // else                  sprintf(filename_pot, "potential_gpugpu_L%d.dat", n->level - 1);
                // status |= save_scalar3d(potential, zdim, ydim, xdim, filename_pot, 100);
                // if(status) return EXIT_FAILURE;
            }
            prev_level = n->level;
            if(verbose_level >= 6) {
                int width = pow(2, actual_limit-n->level);
                printf("    Level%d (%dx%d boxes, size=%dx%d)... ",
                    n->level, (int)pow(2, n->level), (int)pow(2, n->level), width, width);
                fflush(NULL);
            }
        }

        // if(n->is_pruned()) {
            // continue;
        // }

        // char idstring[100];
        // n->get_idstring(idstring);
        // printf("L%d%s(%d,%d)=L%d(%.1f,%.1f) \n", n->level, idstring, n->x, n->y, actual_limit, n->cx, n->cy);


        gettimeofday(&t1, NULL);

        fptype q = 0;

    // Calculate multipole coefficients for the source box
        Cmpx mpc[(P+1)*(P+1)];
        // checking for source charges in the source box
        fptype charge_found = 0;
        fptype width = pow(2, actual_limit-n->level);
        int yy1 = ceil(n->cy-width/2);
        int yy2 = floor(n->cy+width/2);
        #ifdef _OPENMP
        // #pragma omp parallel for
        #endif
        for(int yy=yy1; yy<=yy2; yy++) {
        // for(int yy=ceil(n->cy-width/2); yy<=floor(n->cy+width/2); yy++) {
            for(int xx=ceil(n->cx-width/2); xx<=floor(n->cx+width/2); xx++) {
                q = charge[yy*xdim + xx];
                if(q != 0) { // if charge found
                    charge_found = 1;
                    Cmpx r_(xx - n->cx, yy - n->cy);
                    for(int l=0; l<=P; l++) {
                        for(int m=-l; m<=l; m++) {
                            Cmpx sph = spherical_harmonic(l, m, M_PI/2, r_.get_ang()).conjugate();
                            sph *= q * pow(r_.get_mag(), l);
                            mpc[l*l+l+m] += sph;
                            // mpc[l*l+l+m] += q * pow(r_.get_mag(), l) * spherical_harmonic(l, m, M_PI/2, r_.get_ang()).conjugate();
                        } // m loop
                    } // l loop
                } // if(q != 0)
            } // source charge loop
        } // source charge loop
        // NEWLINE;

        gettimeofday(&t2, NULL);
        deltatime = (t2.tv_sec + t2.tv_usec/1e6) - (t1.tv_sec + t1.tv_usec/1e6);
        t_coeff += deltatime;

        if(! charge_found) {
            // n->prune();
            continue;
        }

        // gettimeofday(&t1, NULL);

        if(charge_found)
        {
            gettimeofday(&t1, NULL);

            // calculation of potential at the boxes in 27 boxes of interaction list
            if(use_gpu && width >= 4) {
                status |= fmm_gpu(  n,
                                    mpc,
                                    potential_gmem, limit, P,
                                    xdim, ydim, zdim, zc,
                                    use_gpu, verbose_level);
            }
            else {
                #ifdef _OPENMP
                // #pragma omp parallel for
                #endif
                for(int i=0; i<27; i++) {
                    Box *ni = n->interaction[i];
                    if(ni != NULL) {
                        for(int yy=ceil(ni->cy-width/2); yy<=floor(ni->cy+width/2); yy++) {
                            for(int xx=ceil(ni->cx-width/2); xx<=floor(ni->cx+width/2); xx++) {
                                for (int zp = 0; zp < zdim; zp++) { // for each potential layer in zdim
                                    Vector3 r(xx - n->cx, yy - n->cy, zp - zc);
                                    Cmpx sum_over_lm;
                                    for(int l=0; l<=P; l++) {
                                        Cmpx sum_over_m;
                                        for(int m=-l; m<=l; m++) {
                                            Cmpx sph = spherical_harmonic(l, m, r.colatitude(), r.azimuth());
                                            sph *= (1.0*factorial(l-abs(m))) / factorial(l+abs(m));
                                            sph *= mpc[l*l+l+m];
                                            sum_over_m += sph;
                                            // sum_over_m += (1.0*factorial(l-abs(m))) / factorial(l+abs(m)) * mpc[l*l+l+m] * spherical_harmonic(l, m, r.colatitude(), r.azimuth());
                                        }
                                        sum_over_m *= 1 / pow(r.magnitude(), l+1);
                                        sum_over_lm += sum_over_m;
                                        // sum_over_lm += 1 / pow(r.magnitude(), l+1) * sum_over_m;
                                    }
                                    potential[zp*ydim*xdim + yy*xdim + xx] += sum_over_lm.get_re();
                                    // potential[yy*xdim+xx] += (sum_over_lm.get_re() > 0) ? sum_over_lm.get_mag() : -sum_over_lm.get_mag();

                                    // const fptype threshold = 1e-2;
                                    // fptype modangle = fabs(sum_over_lm.get_ang());
                                    // modangle = (modangle < M_PI-modangle) ? modangle : M_PI-modangle;
                                    // if(modangle > threshold) {
                                        // if(verbose_level >= 0)
                                            // printf("PANIC!! L%d   R=%g   angle=%g\n", n->level, r.magnitude(), modangle);
                                        // fprintf(paniclog, "%d   %g   %g\n", n->level, r.magnitude(), modangle);
                                    // }
                                } // potential layers
                            }
                        }
                    } // if(ni != NULL)
                } // interaction loop
            } // if(! use_gpu)

            gettimeofday(&t2, NULL);
            deltatime = (t2.tv_sec + t2.tv_usec/1e6) - (t1.tv_sec + t1.tv_usec/1e6);
            t_potential += deltatime;

            // calculation with neighbor list at the deepest level
            if(n->level == actual_limit) {
                gettimeofday(&t1, NULL);
                // printf("nearest potential calulcation.\n");
                // assert(n->cx == n->x && n->cy == n->y);
                // fptype q_prev = q;
                // q = charge[(int)(n->cy*xdim + n->cx)];
                // assert(q == q_prev);

                for (int zp = 0; zp < zdim; zp++) { // for each potential layer in zdim
                    if(zp != zc) { // neighbor on other layers at self position
                        Vector3 r(0, 0, zp - zc);
                        potential[zp*ydim*xdim + (int)(n->cy*xdim + n->cx)] += q / r.magnitude();
                    }
                    for(int i=0; i<8; i++) {
                        Box *nb = n->neighbor[i];
                        if(nb != NULL) {
                            Vector3 r(nb->cx - n->cx, nb->cy - n->cy, zp - zc);
                            potential[zp*ydim*xdim + (int)(nb->cy*xdim + nb->cx)] += q / r.magnitude();
                        }
                    } // neighbor loop

                    gettimeofday(&t2, NULL);
                    deltatime = (t2.tv_sec + t2.tv_usec/1e6) - (t1.tv_sec + t1.tv_usec/1e6);
                    t_potential_nearest += deltatime;
                    // printf("nearest potential calulcation took %f seconds so far.\n", t_potential_nearest);

                } // for each potential layer in zdim
            } // if deepest level

        } // if(charge_found)

        // gettimeofday(&t2, NULL);
        // deltatime = (t2.tv_sec + t2.tv_usec/1e6) - (t1.tv_sec + t1.tv_usec/1e6);
        // t_potential += deltatime;

    } // while(!Q_tree.isEmpty())

    status |= gettimeofday(&time2, NULL);
    deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    if(verbose_level >= 10)
        printf("done in %f seconds.\n", deltatime); fflush(NULL);
    // saving this level potential
    // char filename_pot[200];
    // if     (use_gpu == 0) sprintf(filename_pot, "potential_cpu_L%d.dat", limit);
    // else if(use_gpu == 1) sprintf(filename_pot, "potential_gpu_L%d.dat", limit);
    // else if(use_gpu == 2) sprintf(filename_pot, "potential_gpuemu_L%d.dat", limit);
    // else                  sprintf(filename_pot, "potential_gpugpu_L%d.dat", limit);
    // status |= save_scalar3d(potential, zdim, ydim, xdim, filename_pot, 100);
    // if(status) return EXIT_FAILURE;
    // if(verbose_level >= 10) {
        // printf("done in %f seconds.\n", deltatime);
        // printf("FMM coeff calulcation took %f seconds.\n", t_coeff);
        // printf("FMM potential calulcation took %f seconds.\n", t_potential);
        // printf("nearest potential calulcation took %f seconds.\n", t_potential_nearest);
    // }

    free(queue_mem);
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}





int fmm_calc(   const fptype *charge,
                fptype *potential,
                const int xdim, const int ydim, const int zdim,
                const int P,    // multipole series truncation (l = 0...P)
                const int use_gpu,
                const int verbose_level )
{
    static int first_time = 1;
    int status = 0;
    const unsigned int logN = ceil(log2f(xdim * ydim) / log2f(4));
    // FILE *paniclog = fopen("paniclog.dat", "w");
    FILE *paniclog = fopen("paniclog.dat", "a");
    fprintf(paniclog, "# FMM: New run\n");

    timeval time1, time2;
    double deltatime;
    status |= gettimeofday(&time1, NULL);

// allocate memory for the Tree and its associated BFS Queue
    int total_boxes = (4*xdim*ydim - 1) / 3;
    printf("sizeof(Box) = %lu\n", sizeof(Box));
    printf("total Boxes in the tree = %d\n", total_boxes);
    printf("memory required for the tree = %lu Bytes\n", total_boxes * sizeof(Box));
    // printf("memory required for the tree = %.0f KB\n", total_boxes * sizeof(Box) / 1024.0);
    Box *tree = (Box*)malloc(total_boxes * sizeof(Box));
    int len = xdim * ydim;
    // void **contents_ = new void*[len]();
    void **queue_mem = (void**)malloc(len * sizeof(void*));
    if(tree == NULL || queue_mem == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

// generate the tree
    tree[0] = Box(0, 0, logN);
    Box *root = &tree[0];
    // Box *root = new Box(0, 0, logN);
    // root->create_tree_recurse(logN);
    root->create_tree_bfs(logN, queue_mem);
    root->find_neighbors_recurse(root, logN);
    free(queue_mem);

    status |= gettimeofday(&time2, NULL);
    deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    if(verbose_level >= 0)
        printf("Tree: took %f seconds\n", deltatime);
    fflush(NULL);


    // timeval time1, time2;
    status |= gettimeofday(&time1, NULL);

// pinned host memory and associated device memory
    static fptype *potential_pinned = NULL;
    static fptype *potential_gmem = NULL;
    if(first_time) {
        cudaHostAlloc((void **)&potential_pinned, zdim*ydim*xdim * sizeof(fptype), cudaHostAllocMapped);
        checkCUDAError("cudaHostAllocMapped");
        // Get the device pointers to the mapped memory
        cudaHostGetDevicePointer((void **)&potential_gmem, (void *)potential_pinned, 0);
        checkCUDAError("cudaHostGetDevicePointer");
    }
// reset potential before beginning
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < zdim*ydim*xdim; i++)
        potential_pinned[i] = 0;

// reset potential before beginning
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for(int i = 0; i < zdim*ydim*xdim; i++)
        potential[i] = 0;

// for each charge layer in zdim
    for (int zc = 0; zc < zdim; zc++) {
        if(verbose_level >= 3)
            printf("  FMM: charge layer %d\n", zc);
        fprintf(paniclog, "# FMM:   charge layer %d\n", zc);
        fflush(NULL);
        // call the actual function
        status |= fmm_bfs(charge+zc*ydim*xdim, potential_pinned, potential_gmem, root, logN, logN, P, xdim, ydim, zdim, zc, paniclog, use_gpu, verbose_level);

        if(status) return EXIT_FAILURE;
        // root->grow();
    }

    memcpy(potential, potential_pinned, zdim*ydim*xdim*sizeof(fptype));

    if(status) return EXIT_FAILURE;
    status |= gettimeofday(&time2, NULL);
    deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    if(verbose_level >= 0)
        printf("FMM: took %f seconds\n", deltatime);
    fflush(NULL);

// closing
    fclose(paniclog);
    // delete root;
    free(tree);
    first_time = 0;
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}


// Exact O(N^2) calculation of potential
void calc_potential_exact( const fptype *charge,
                        const int xdim, const int ydim, const int zdim,
                        fptype *potential, int use_gpu)
{
    int status = 0;
    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);
    if(use_gpu)
        calc_potential_exact_gpu(charge, xdim, ydim, zdim, potential);
    else
    {
    // reset potential before beginning
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for(int i = 0; i < zdim*ydim*xdim; i++)
            potential[i] = 0;

        for(int z_ = 0; z_ < zdim; z_++) {  // source loop
            for(int y_ = 0; y_ < ydim; y_++) {
                for(int x_ = 0; x_ < xdim; x_++) {
                    fptype q = charge[z_*ydim*xdim + y_*xdim + x_];
                    if(q == 0) continue;
                    for(int z = 0; z < zdim; z++) { // observation point loop
                        #ifdef _OPENMP
                        #pragma omp parallel for
                        #endif
                        for(int y = 0; y < ydim; y++) {
                            for(int x = 0; x < xdim; x++) {
                                if(z == z_ && y == y_ && x == x_) continue;    // skip on itself
                                Vector3 R(x-x_, y-y_, z-z_);
                                potential[z*ydim*xdim + y*xdim + x] += q / R.magnitude();
                            }
                        }
                    }
                }
            }
        }
    }
    status |= gettimeofday(&time2, NULL);
    double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    if(0)
        printf("Exact: took %f seconds\n", deltatime);
    fflush(NULL);
}


// H field based on nearest neighbor coupling only
void calc_H_nearest_neighbor(   const Vector3 *M, Vector3 *H,
                                const int xdim, const int ydim, const int zdim )
{
    for(int i = 0; i < zdim*ydim*xdim; i++)
        H[i] = -0.2 * (   ((i-xdim >= 0)        ? M[i-xdim] : Vector3(0,0,0))     // top
                        + ((i+xdim < ydim*xdim) ? M[i+xdim] : Vector3(0,0,0))     // bottom
                        + ((i%xdim != 0)        ? M[i-1]    : Vector3(0,0,0))     // left
                        + (((i+1)%xdim != 0)    ? M[i+1]    : Vector3(0,0,0)) );  // right
}
