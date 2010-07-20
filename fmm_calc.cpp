#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
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
        // populate queue with children nodes        
        if (n->level < limit)
            for(int i=0; i<=3; i++)
                Q.enqueue(n->child[i]);
    }
}

// FMM algorithm in BFS
// ===============================
void fmm_bfs(   Box *n, 
                    const unsigned int limit, 
                    const unsigned int actual_limit, 
                    const unsigned int H, 
                    const float *charge, 
                    const unsigned int P    // multipole series (k = 1...P)
                ) 
{
    assert(limit <= actual_limit);
    printf("Executing FMM algorithm...\n");
    unsigned int prev_level = 0;

    const unsigned int N = (unsigned int)pow(4, limit);
    Queue Q(N);
    Q.enqueue(n);
    while(!Q.isEmpty()) {
        Box *n = (Box*)Q.dequeue();
        // function to perform on node
        if(n->level >= 2)   // no FMM steps for Level-0 and Level-1
        {
            // inheret with parent's potential
            n->potential += n->parent->potential;            
            float Q = 0;
            float *a = new float[P]();
            float width = pow(2, actual_limit-n->level);

            if(prev_level != n->level) {
                printf("   level%d, %d boxes...\n", n->level, (int)pow(4, n->level));
                prev_level = n->level;
            }
            // char idstring[100];
            // n->get_idstring(idstring);
            // printf("L%d%s(%d,%d)=L%d(%.1f,%.1f) ", n->level, idstring, n->x, n->y, actual_limit, n->cx, n->cy);
            // printf("width=%d ", (int)width);
            // printf("Xrange=[%.0f,%.0f] Yrange=[%.0f,%.0f]", ceil(n->cx-width/2), floor(n->cx+width/2), 
                                                            // ceil(n->cy-width/2), floor(n->cy+width/2));
            for(int yy=ceil(n->cy-width/2); yy<=floor(n->cy+width/2); yy++) {
                for(int xx=ceil(n->cx-width/2); xx<=floor(n->cx+width/2); xx++) {
                    float q = charge[yy*H+xx];
                    if (q != 0) {
                        Q += q;
                        float rx_ = xx - n->cx;
                        float ry_ = yy - n->cy;
                        // printf("{q=%g@(%d,%d)=[%.1f%s%.1fi]} ", q, xx,yy, rx_, ry_>=0 ? "+" : "", ry_);
                        // printf("ak = {");
                        for(unsigned int k=1; k<=P; k++) {
                            a[k-1] += -q * pow(magnitude(rx_, ry_), k);
                            // printf("%g, ", a[k-1]);
                        }
                        // printf("}");
                    }
                }
            }
            // if(0 || Q !=0) {
                // printf("Q=%g ", Q);
                // printf("ak={");
                // for(unsigned int k=1; k<=P; k++)
                    // printf("%g, ", a[k-1]);
                // printf("} ");
            // }

        // calculation with interaction list 
            float rx;
            float ry;
            // printf("interacts with ");
            for(int i=0; i<27; i++) {
                Box *ni = n->interaction[i];
                if (ni != NULL) 
                {
                    float potential_tmp = 0;
                    rx = ni->cx - n->cx;
                    ry = ni->cy - n->cy;
                    // ni->get_idstring(idstring);
                    // NEWLINE;
                    // printf("    ");
                    // printf("L%d%s: ", ni->level, idstring);
                    // printf("[%.0f%s%.0fi]", rx, ry>=0 ? "+" : "", ry);
                    potential_tmp += Q * log(magnitude(rx, ry));
                    for(unsigned int k=1; k<=P; k++) {
                        potential_tmp += a[k-1] / pow(magnitude(rx, ry), k) / k;
                    }
                    // printf("%.3g + ", potential_tmp);
                    // printf("%.3g + ", ni->parent->potential);
                    // printf("%.3g = ", ni->potential);
                    // inheret potential value from the parent
                    ni->potential += potential_tmp;
                    // ni->potential += potential_tmp + ni->parent->potential;
                    // printf("%g", ni->potential);
                }
            }
            // NEWLINE;

        // calculation with neighbor list at the deepest level
            if(n->level == limit) {
                for(int i=0; i<8; i++) {
                    Box *nb = n->neighbor[i];
                    if (nb != NULL) 
                    {
                        rx = nb->cx - n->cx;
                        ry = nb->cy - n->cy;
                        nb->potential += Q * log(magnitude(rx, ry));
                    }
                }
            }
        }   // if(n->level >= 2)

        // populate queue with children nodes
        if (n->level < limit)
            for(int i=0; i<=3; i++)
                Q.enqueue(n->child[i]);
    } // while(!Q.isEmpty())
}


int fmm_calc(Box *root, const unsigned int limit) {
    const unsigned int N = (unsigned int)pow(4, limit);
    const unsigned int H = (unsigned int)sqrt(N);
    const unsigned int P = 3;   // multipole series truncation (k = 1...P)
    const float charge_probability = 0.001;
    printf("=================\n");    
    printf("|  fmm_calc()   |\n");
    printf("=================\n");    
    printf("N = %d, log4(N) = %d, sqrt(N) = %d, P = %d\n", N, limit, H, P);
    int status = 0;

    // printf("#################\n");    
    // printf("#      BFS      #\n");
    // printf("#################\n");    
    // traverse_tree_bfs(root, limit);
    // printf("#################\n");    
    // printf("#      DFS      #\n");
    // printf("#################\n");    
    // traverse_tree_dfs(root, limit);
    // return 0;



// charge matrix
    float *charge = new float[N]();
    

// Charge configuration
// ===========================
// dipole
    // int l = H/8;
    // int x0 = H/2;
    // int y1 = (H-l)/2;
    // int y2 = (H+l)/2;
    // charge[y1*H + x0] = 1;
    // charge[y2*H + x0] = -1;

// place a random charge at random location
    int xx = H/2;
    int yy = H/2;
    charge[yy*H + xx] = 1;
    charge[rand_atob(0,H)*H + rand_atob(0,H)] = 1;
    charge[rand_atob(0,H)*H + rand_atob(0,H)] = 1;
    // for(int i=0; i<(int)(N*charge_probability); i++)
        // charge[rand_atob(0,H)*H + rand_atob(0,H)] = 1;
        // charge[rand_atob(0,H)*H + rand_atob(0,H)] = frand_atob(1, 10);

// Few random charges at random locations
    // for(unsigned int yy=0; yy<H; yy++) {
        // for(unsigned int xx=0; xx<H; xx++) {
            // if(frand_atob(0, 1) < charge_probability)
                // charge[yy*H + xx] = 1;
                // // charge[yy*H + xx] = frand_atob(1, 10);
        // }
    // }

// write potential matrix to file
    status |= matrix2file(charge, H, H, "charge.dat");

// Run the FMM algorithm on tree
    fmm_bfs(root, limit, limit, H, charge, P);
        
// collect result in a matrix
    float *potential = new float[N]();    // potential matrix
    collect_potential_bfs(root, limit-0, limit, H, potential);
// write potential matrix to file
    status |= matrix2file(potential, H, H, "potential.dat");

// collect successive stages of FMM in matrices and write them to files
// for pedagogical purpose only
    for(unsigned int l=2; l<=limit; l++) {
        char filename[100];
        collect_potential_bfs(root, l, limit, H, potential);
        sprintf(filename, "potential_L%d.dat", l);
        status |= matrix2file(potential, H, H, filename);
    }


// closing
    delete []potential;
    delete []charge;
    return status;
}
