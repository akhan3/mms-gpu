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
#include "helper_functions.hpp"


// FMM algorithm in BFS
// ===============================
int fmm_bfs(        const fptype *charge,
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
    int status = 0;
    assert(limit == actual_limit);  // unable to support arbitrary depth calculations.
    // assert(limit <= actual_limit);
    // timeval time1, time2;
    // status |= gettimeofday(&time1, NULL);
    if(verbose_level >= 6)
        printf("Executing FMM algorithm...\n");
    // unsigned int prev_level = 0;

    const unsigned int N = (unsigned int)pow(4, limit);
    Queue Q_tree(N);
    Q_tree.enqueue((void*)root);
    while(!Q_tree.isEmpty()) {
        Box *n = (Box*)Q_tree.dequeue();
        // populate queue with children nodes
        if (n->level < limit)
            for(int i=0; i<=3; i++)
                Q_tree.enqueue(n->child[i]);
        // function to perform on node
        if(n->level >= 2)   // no FMM steps for Level-0 and Level-1
        {
            // if(prev_level != n->level) {
                // if(prev_level >= 2) {
                    // status |= gettimeofday(&time2, NULL);
                    // double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
                    // status |= gettimeofday(&time1, NULL);
                    // if(verbose_level >= 10)
                        // printf("done in %f seconds.\n", deltatime); fflush(NULL);
                // }
                // prev_level = n->level;
                // if(verbose_level >= 6)
                    // printf("Level%d (%d boxes)... ", n->level, (int)pow(4, n->level)); fflush(NULL);
            // }

            if(n->is_pruned()) {
                continue;
            }

            // char idstring[100];
            // n->get_idstring(idstring);
            // printf("L%d%s(%d,%d)=L%d(%.1f,%.1f) \n", n->level, idstring, n->x, n->y, actual_limit, n->cx, n->cy);

        // Calculate multipole moments for the source box
            Cmpx multipole_coeff[P+1][2*P+1];
            // checking for source charges in the source box
            fptype charge_found = 0;
            fptype width = pow(2, actual_limit-n->level);
            for(int l=0; l<=P; l++) {
                for(int m=-l; m<=l; m++) {
                    for(int yy=ceil(n->cy-width/2); yy<=floor(n->cy+width/2); yy++) {
                        for(int xx=ceil(n->cx-width/2); xx<=floor(n->cx+width/2); xx++) {
                            fptype q = charge[yy*xdim + xx];
                            if (q != 0) { // if charge found
                                charge_found = 1;
                                Cmpx r_(xx - n->cx, yy - n->cy, 0);
                                multipole_coeff[l][m+l] += q * pow(r_.get_mag(), l) * spherical_harmonic(l, m, M_PI/2, r_.get_ang()).conjugate();
                            } // if (q != 0)
                        } // source charge loop
                    } // source charge loop
                } // m loop
            } // l loop
            // NEWLINE;

            if(! charge_found) {
                n->prune();
                continue;
            }

            if(charge_found)
            {
                for (int zp = 0; zp < zdim; zp++) // for each potential layer in zdim
                {
                    // printf("FMM:   charge layer %d, potential layer %d\n", zc, zp);
                // calculation of potential at the boxes in interaction list
                    for(int i=0; i<27; i++) {
                        Box *ni = n->interaction[i];
                        if (ni != NULL) {
                            for(int yy=ceil(ni->cy-width/2); yy<=floor(ni->cy+width/2); yy++) {
                                for(int xx=ceil(ni->cx-width/2); xx<=floor(ni->cx+width/2); xx++) {
                                    // Cmpx r(xx - n->cx, yy - n->cy, 0);
                                    Vector3 r(xx - n->cx, yy - n->cy, zp - zc);
                                    Cmpx sum_over_lm;
                                    for(int l=0; l<=P; l++) {
                                        Cmpx sum_over_m;
                                        for(int m=-l; m<=l; m++) {
                                            sum_over_m += (1.0*factorial(l-abs(m))) / factorial(l+abs(m)) * multipole_coeff[l][m+l] * spherical_harmonic(l, m, r.colatitude(), r.azimuth());
                                        }
                                        sum_over_lm += 1 / pow(r.magnitude(), l+1) * sum_over_m;
                                    }
                                    potential[zp*ydim*xdim + yy*xdim + xx] += sum_over_lm.get_re();
                                    // potential[yy*xdim+xx] += (sum_over_lm.get_re() > 0) ? sum_over_lm.get_mag() : -sum_over_lm.get_mag();

                                    const fptype threshold = 1e-2;
                                    fptype modangle = fabs(sum_over_lm.get_ang());
                                    modangle = (modangle < M_PI-modangle) ? modangle : M_PI-modangle;
                                    if(modangle > threshold) {
                                        if(verbose_level >= 0)
                                            printf("PANIC!! L%d   R=%g   angle=%g\n", n->level, r.magnitude(), modangle);
                                        fprintf(paniclog, "%d   %g   %g\n", n->level, r.magnitude(), modangle);
                                    }
                                }
                            }
                        } // if (ni != NULL)
                    } // interaction loop

                // calculation with neighbor list at the deepest level
                    if(n->level == actual_limit) {
                        assert(n->cx == n->x && n->cy == n->y);
                        fptype q = charge[(int)(n->cy*xdim + n->cx)];
                        if(zp != zc) { // neighbor on other layers at self position
                            Vector3 r(0, 0, zp - zc);
                            potential[zp*ydim*xdim + (int)(n->cy*xdim + n->cx)] += q / r.magnitude();
                        }
                        for(int i=0; i<8; i++) {
                            Box *nb = n->neighbor[i];
                            if (nb != NULL) {
                                Vector3 r(nb->cx - n->cx, nb->cy - n->cy, zp - zc);
                                potential[zp*ydim*xdim + (int)(nb->cy*xdim + nb->cx)] += q / r.magnitude();
                            }
                        } // neighbor loop
                    } // if deepest level
                } // for each potential layer in zdim
            } // if(charge_found)
        }   // if(n->level >= 2)
    } // while(!Q_tree.isEmpty())

    // status |= gettimeofday(&time2, NULL);
    // double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    // if(verbose_level >= 10)
        // printf("done in %f seconds.\n", deltatime);
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}


int fmm_calc(   const fptype *charge,
                fptype *potential,
                const int xdim, const int ydim, const int zdim,
                const int P,    // multipole series truncation (l = 0...P)
                const int verbose_level )
{
    int status = 0;
    const unsigned int logN = ceil(log2(xdim * ydim) / log2(4));
    // FILE *paniclog = fopen("paniclog.dat", "w");
    FILE *paniclog = fopen("paniclog.dat", "a");
    fprintf(paniclog, "# FMM: New run\n");

// generate the tree
    Box *root = new Box(0, logN);
    create_tree_recurse(root, logN);
    find_neighbors_recurse(root, root, logN);

// for each charge layer in zdim
    for (int zc = 0; zc < zdim; zc++) {
        if(verbose_level >= 3)
            printf("FMM: charge layer %d\n", zc);
        fprintf(paniclog, "# FMM:   charge layer %d\n", zc);
        fflush(NULL);
        timeval time1, time2;
        status |= gettimeofday(&time1, NULL);
        // call the actual function
        status |= fmm_bfs(charge+zc*ydim*xdim, potential, root, logN, logN, P, xdim, ydim, zdim, zc, paniclog, verbose_level);
        status |= gettimeofday(&time2, NULL);
        double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
        if(verbose_level >= 5)
            printf("FMM: charge layer %d took %f seconds\n", zc, deltatime);
        if(status) return EXIT_FAILURE;
        root->grow();
    }

// closing
    fclose(paniclog);
    delete root;
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}


// Exact O(N^2) calculation of potential
void calc_potential_exact( const fptype *charge,
                        const int xdim, const int ydim, const int zdim,
                        fptype *potential)
{
    for(int z_ = 0; z_ < zdim; z_++) {  // source loop
        for(int y_ = 0; y_ < ydim; y_++) {
            for(int x_ = 0; x_ < xdim; x_++) {
                fptype q = charge[z_*ydim*xdim + y_*xdim + x_];
                if (q == 0) continue;
                for(int z = 0; z < zdim; z++) { // observation point loop
                    #ifdef _OPENMP
                    #pragma omp parallel for
                    #endif
                    for(int y = 0; y < ydim; y++) {
                        for(int x = 0; x < xdim; x++) {
                            if (z == z_ && y == y_ && x == x_) continue;    // skip on itself
                            Vector3 R(x-x_, y-y_, z-z_);
                            potential[z*ydim*xdim + y*xdim + x] += q / R.magnitude();
                        }
                    }
                }
            }
        }
    }
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
