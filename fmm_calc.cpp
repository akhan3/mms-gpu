#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <cmath>
#include "Box.hpp"
#include "Queue.hpp"
#include "helper_functions.hpp"
#define NEWLINE printf("\n");


// Collect potential data BFS
// ===============================
void collect_potential_bfs( Box *n,
                            const unsigned int limit,
                            const unsigned int actual_limit,
                            const unsigned int H,
                            float *potential
                          )
{
    assert(limit <= actual_limit);
    printf("Collecting potential data at L%d and smearing it down to L%d...\n", limit, actual_limit);
    const unsigned int N = (unsigned int)pow(4, limit);
    Queue Q(N);
    Q.enqueue(n);
    while(!Q.isEmpty()) {
        Box *n = (Box*)Q.dequeue();
        // populate queue with children nodes
        if (n->level < limit)
            for(int i=0; i<=3; i++)
                Q.enqueue(n->child[i]);
        // function to perform on node
        if (n->level == limit) {
            char idstring[100];
            n->get_idstring(idstring);
            // smear this potential value to all deeper levels
            int width = (int)pow(2, actual_limit - limit);
            // printf("L%d%s(%d,%d)=L%d(%.1f,%.1f) ", n->level, idstring, n->x, n->y, actual_limit, n->cx, n->cy);
            // printf("width=%d :: ", (int)width);
            for(unsigned int y = width*n->y; y < width*(n->y+1); y++) {
                for(unsigned int x = width*n->x; x < width*(n->x+1); x++) {
                    potential[(int)(y*H+x)] = n->potential;
                    // printf("(%d,%d)", x,y);
                }
            }
            // NEWLINE;
        }
    }
}

// FMM algorithm in BFS
// ===============================
int fmm_bfs(   Box *root,
                    const unsigned int limit,
                    const unsigned int actual_limit,
                    const unsigned int H,
                    const float *charge,
                    const unsigned int P    // multipole series (k = 1...P)
                )
{
    int status = 0;
    assert(limit == actual_limit);  // unable to currently support arbitrary depth calculations.
    // assert(limit <= actual_limit);
    printf("Executing FMM algorithm...\n");
    unsigned int prev_level = 0;

    const unsigned int N = (unsigned int)pow(4, limit);
    const float lengthBox = 10;
    const float meshwidth = lengthBox / sqrt(N);
    Queue Q_tree(N);
    Q_tree.enqueue(root);
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
                printf("   level%d, %d boxes...\n", n->level, (int)pow(4, n->level));
                prev_level = n->level;
            }
            float width = pow(2, actual_limit-n->level);
            // char idstring[100];
            // n->get_idstring(idstring);
            // printf("L%d%s(%d,%d)=L%d(%.1f,%.1f) ", n->level, idstring, n->x, n->y, actual_limit, n->cx, n->cy);

        // calculation with interaction list
            for(int i=0; i<27; i++) {
                Box *ni = n->interaction[i];
                if (ni != NULL) {
                    Queue Q_nic((int)(N / pow(4,n->level)));
                    Q_nic.enqueue(ni);
                    while(!Q_nic.isEmpty()) {
                        Box *nic = (Box*)Q_nic.dequeue();
                        // populate queue with children nodes
                        if (nic->level < limit)
                            for(int i=0; i<=3; i++)
                                Q_nic.enqueue(nic->child[i]);
                        // function to perform on node
                        if (nic->level == limit) {
                            float rx = nic->cx - n->cx;
                            float ry = nic->cy - n->cy;
                            float potential_tmp = 0;
                            if(n->level < limit) // if (except last level), do the multipole expansion
                            {
                                // checking for source charges in the source box
                                for(int yy=ceil(n->cy-width/2); yy<=floor(n->cy+width/2); yy++) {
                                    for(int xx=ceil(n->cx-width/2); xx<=floor(n->cx+width/2); xx++) {
                                        float q = charge[yy*H+xx];
                                        if (q != 0) { // if charge found
                                            float rx_ = xx - n->cx;
                                            float ry_ = yy - n->cy;
                                            float cos_dtheta = cmpx_costheta_between(rx,ry,rx_,ry_);
                                            // printf("cos_theta=%f ", cos_dtheta);
                                            float multipole_series = 0;
                                            // multipole series
                                            for(unsigned int k=0; k<=P; k++) {
                                                float num = pow(cmpx_magnitude(rx_, ry_), k);
                                                float den = pow(cmpx_magnitude(rx, ry), k+1);
                                                float leg = legendre(k, cos_dtheta);
                                                multipole_series += num / den * leg;
                                            }
                                            potential_tmp += q * multipole_series / meshwidth;
                                        } // if (q != 0)
                                    } // source charge loop
                                } // source charge loop
                            } // if (except last level)
                            else // if (last level), multipole expansions is not possible.
                                 // do exact calculation instead
                            {
                                assert(n->cx == n->x && n->cy == n->y);
                                float q = charge[(int)(n->cy*H + n->cx)];
                                potential_tmp += q / cmpx_magnitude(rx, ry) / meshwidth;
                            } // if (last level)
                            nic->potential += potential_tmp;
                        } // if (nic->level == limit)
                    } // while(!Q_nic.isEmpty())
                } // if (ni != NULL)
            } // interaction loop
            // NEWLINE;

        // calculation with neighbor list at the deepest level
            if(n->level == limit) {
                for(int i=0; i<8; i++) {
                    Box *nb = n->neighbor[i];
                    if (nb != NULL)
                    {
                        assert(n->cx == n->x && n->cy == n->y);
                        float q = charge[(int)(n->cy*H + n->cx)];
                        float rx = nb->cx - n->cx;
                        float ry = nb->cy - n->cy;
                        nb->potential += q / cmpx_magnitude(rx, ry) / meshwidth;
                    }
                } // neighbor loop
            } // if deepest level
        }   // if(n->level >= 2)
    } // while(!Q_tree.isEmpty())
    return status;
}


int fmm_calc(Box *root, const unsigned int limit) {
    const unsigned int N = (unsigned int)pow(4, limit);
    const unsigned int H = (unsigned int)sqrt(N);
    const unsigned int P = 3;   // multipole series truncation (k = 1...P)
    printf("=================\n");
    printf("|  fmm_calc()   |\n");
    printf("=================\n");
    printf("N = %d, log4(N) = %d, sqrt(N) = %d, P = %d\n", N, limit, H, P);
    int status = 0;

// charge matrix
    float *charge = new float[N]();

// Charge configuration
// ===========================
// // dipole
    // int l = H/8;
    // int x0 = H/2;
    // int y1 = (H-l)/2;
    // int y2 = (H+l)/2;
    // charge[y1*H + x0] = 1;
    // charge[y2*H + x0] = -1;

// place a random charge at random location
    // int xx = H/2;
    // int yy = H/2;
    // // charge[yy*H + xx] = 1;
    // charge[yy*H + xx] = 1/4.0;
    // charge[(yy-1)*H + xx] = 1/4.0;
    // charge[yy*H + xx-1] = 1/4.0;
    // charge[(yy-1)*H + xx-1] = 1/4.0;

// Few random charges at random locations
    const float charge_probability = 0.005;
    for(unsigned int yy=0; yy<H; yy++)
        for(unsigned int xx=0; xx<H; xx++)
            if(frand_atob(0, 1) < charge_probability)
                charge[yy*H + xx] = frand_atob(1, 10);

// write charge matrix to file
    status |= matrix2file(charge, H, H, "charge.dat");
    // status |= matrix2stdout(charge, H, H, "charge");


// Run the FMM algorithm on tree
    timeval time1, time2;
    status |= gettimeofday(&time1, NULL);
        status |= fmm_bfs(root, limit, limit, H, charge, P);
    status |= gettimeofday(&time2, NULL);
    double deltatime = (time2.tv_sec + time2.tv_usec/1e6) - (time1.tv_sec + time1.tv_usec/1e6);
    printf("fmm_bfs() took %f seconds\n", deltatime);



// collect result in a matrix
    float *potential = new float[N]();    // potential matrix
    collect_potential_bfs(root, limit-0, limit, H, potential);
// write potential matrix to file
    status |= matrix2file(potential, H, H, "potential.dat");

// // collect successive stages of FMM in matrices and write them to files
// // for pedagogical purpose only
    // for(unsigned int l=2; l<=limit; l++) {
        // char filename[100];
        // collect_potential_bfs(root, l, limit, H, potential);
        // sprintf(filename, "potential_L%d.dat", l);
        // status |= matrix2file(potential, H, H, filename);
    // }

// closing
    delete []potential;
    delete []charge;
    return status;
}
