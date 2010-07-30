#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include "Box.hpp"
#include "Queue.hpp"
#include "Cmpx.hpp"
#include "helper_functions.hpp"
#define NEWLINE printf("\n");


// FMM algorithm in BFS
// ===============================
int fmm_bfs(        const Box *root,
                    const unsigned int limit,
                    const unsigned int actual_limit,
                    const unsigned int H,
                    const float *charge,
                    const int P,    // multipole series (l = 0...P)
                    float *potential
                )
{
    int status = 0;
    assert(limit == actual_limit);  // unable to currently support arbitrary depth calculations.
    // assert(limit <= actual_limit);
    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);
    printf("Executing FMM algorithm...\n");
    unsigned int prev_level = 0;

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
            if(prev_level != n->level) {
                if(prev_level >= 2) {
                    status |= gettimeofday(&time2, NULL);
                    double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
                    status |= gettimeofday(&time1, NULL);
                    printf("done in %f seconds.\n", deltatime); fflush(NULL);
                    // collect successive stages of FMM in matrices and write them to files for pedagogical purpose only
                    // char filename[100];
                    // sprintf(filename, "potential_L%d.dat", prev_level);
                    // status |= matrix2file(potential, H, H, filename);
                    // if(status) return EXIT_FAILURE;
                }
                prev_level = n->level;
                printf("Level%d (%d boxes)... ", n->level, (int)pow(4, n->level)); fflush(NULL);
            }

            if(n->is_pruned()) {
                continue;
            }

            // char idstring[100];
            // n->get_idstring(idstring);
            // printf("L%d%s(%d,%d)=L%d(%.1f,%.1f) \n", n->level, idstring, n->x, n->y, actual_limit, n->cx, n->cy);

        // Calculate multipole moments for the source box
            Cmpx multipole_coeff[P+1][2*P+1];
            // checking for source charges in the source box
            float charge_found = 0;
            float width = pow(2, actual_limit-n->level);
            for(int l=0; l<=P; l++) {
                for(int m=-l; m<=l; m++) {
                    for(int yy=ceil(n->cy-width/2); yy<=floor(n->cy+width/2); yy++) {
                        for(int xx=ceil(n->cx-width/2); xx<=floor(n->cx+width/2); xx++) {
                            float q = charge[yy*H+xx];
                            if (q != 0) { // if charge found
                                charge_found = 1;
                                Cmpx r_ = Cmpx(xx - n->cx, yy - n->cy, 0);
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
            // calculation of potential at the boxes in interaction list
                for(int i=0; i<27; i++) {
                    Box *ni = n->interaction[i];
                    if (ni != NULL) {
                        for(int yy=ceil(ni->cy-width/2); yy<=floor(ni->cy+width/2); yy++) {
                            for(int xx=ceil(ni->cx-width/2); xx<=floor(ni->cx+width/2); xx++) {
                                Cmpx r = Cmpx(xx - n->cx, yy - n->cy, 0);
                                Cmpx sum_over_lm;
                                for(int l=0; l<=P; l++) {
                                    Cmpx sum_over_m;
                                    for(int m=-l; m<=l; m++) {
                                        sum_over_m += (1.0*factorial(l-abs(m))) / factorial(l+abs(m)) * multipole_coeff[l][m+l] * spherical_harmonic(l, m, M_PI/2, r.get_ang());
                                    }
                                    sum_over_lm += 1 / pow(r.get_mag(), l+1) * sum_over_m;
                                }
                                potential[yy*H+xx] += sum_over_lm.get_re();

                                const float threshold = 1e-2;
                                float modangle = fabs(sum_over_lm.get_ang());
                                modangle = (modangle < M_PI-modangle) ? modangle : M_PI-modangle;
                                if(modangle > threshold)
                                    printf("PANIC!! sum_over_lm.ang=%g , modangle=%g \n", sum_over_lm.get_ang(), modangle);
                                // assert(fabs(sum_over_lm.get_ang()) <= threshold || (M_PI-fabs(sum_over_lm.get_ang())) <= threshold);   // make sure there is no imaginary part remaining
                            }
                        }
                    } // if (ni != NULL)
                } // interaction loop

            // calculation with neighbor list at the deepest level
                if(n->level == actual_limit) {
                    for(int i=0; i<8; i++) {
                        Box *nb = n->neighbor[i];
                        if (nb != NULL)
                        {
                            assert(n->cx == n->x && n->cy == n->y);
                            float q = charge[(int)(n->cy*H + n->cx)];
                            Cmpx r = Cmpx(nb->cx - n->cx, nb->cy - n->cy, 0);
                            potential[(int)(nb->cy*H+nb->cx)] += q / r.get_mag();
                        }
                    } // neighbor loop
                } // if deepest level
            } // if(charge_found)
        }   // if(n->level >= 2)
    } // while(!Q_tree.isEmpty())

    status |= gettimeofday(&time2, NULL);
    double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    printf("done in %f seconds.\n", deltatime);
    // collect successive stages of FMM in matrices and write them to files for pedagogical purpose only
    // char filename[100];
    // sprintf(filename, "potential_L%d.dat", prev_level);
    // status |= matrix2file(potential, H, H, filename);
    // if(status) return EXIT_FAILURE;
    return status;
}


int fmm_calc(const Box *root, const unsigned int limit, const float *charge, float *potential) {
    const unsigned int N = (unsigned int)pow(4, limit);
    const unsigned int H = (unsigned int)sqrt(N);
    const          int P = 3;   // multipole series truncation (k = 1...P)
    printf("=================\n");
    printf("|  fmm_calc()   |\n");
    printf("=================\n");
    printf("N = %d, log4(N) = %d, sqrt(N) = %d, P = %d\n", N, limit, H, P);
    int status = 0;

// Run the FMM algorithm on tree
    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);
        status |= fmm_bfs(root, limit, limit, H, charge, P, potential);
    status |= gettimeofday(&time2, NULL);
    double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    printf("fmm_bfs() took %f seconds\n", deltatime);

// closing
    return status;
}
