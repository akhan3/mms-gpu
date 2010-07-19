#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <cmath>
#include "Box.hpp"
#define NEWLINE printf("\n");

int rand_atob(const int a, const int b);
float frand_atob(const float a, const float b);

using namespace std;
using std::cout;
using std::endl;

// =========================================================================
// ==== Complex number functions ===========================================
// =========================================================================
void cmpx_pow(  const float x, const float y, const unsigned int p,
                float *xo, float *yo)
{
    *xo = 1;
    *yo = 0;
    for(unsigned int i=1; i<=p; i++) {
        float x1;
        float y1;
        x1 = x*(*xo) - y*(*yo);
        y1 = y*(*xo) + x*(*yo);
        *xo = x1;
        *yo = y1;
    }
}

void cmpx_ln(  const float x, const float y,
                float *xo, float *yo)
{
    float r = sqrt(x*x + y*y);
    float theta = atan2(y, x);
    *xo = log(r);
    *yo = theta;
}

void cmpx_div(  const float a, const float b,
                const float c, const float d,
                float *xo, float *yo)
{
    float denominator = c*c + d*d;
    *xo = (a*c + b*d) / denominator;
    *yo = (b*c - a*d) / denominator;
}
// =========================================================================


// Breadth-first tarversal
// ===============================
// void traverse_tree_bfs( Box *n, 
                            // const unsigned int limit, 
                            // const unsigned int H, 
                            // float *potential
                          // ) 
// {
    
    // // function to perform on node
    // if (n->level == limit)
        // potential[(int)(n->cy*H+n->cx)] = n->potential;
    // // Recursion
    // if (n->level < limit) {
        // for(int i=0; i<=3; i++)
            // traverse_tree_recurse(n->child[i], limit, H, potential);
    // }
// }

// to collect potential data
// ===============================
void traverse_tree_recurse( Box *n, 
                            const unsigned int limit, 
                            const unsigned int H, 
                            float *potential
                          ) 
{
    // function to perform on node
    if (n->level == limit)
        potential[(int)(n->cy*H+n->cx)] = n->potential;
    // Recursion
    if (n->level < limit) {
        for(int i=0; i<=3; i++)
            traverse_tree_recurse(n->child[i], limit, H, potential);
    }
}

// to run FMM algo
// ===============================
void fmm_recurse(   Box *n, 
                    const unsigned int limit, 
                    const unsigned int H, 
                    const unsigned int N,
                    const float *charge, 
                    const unsigned int P    // multipole series (k = 1...P)
                ) 
{
// no FMM steps for Level-0 and Level-1
    if(n->level >= 2)   // but this condition is not necessary
    {
    // function to perform on node
        float Q = 0;
        float *a_r = new float[P]();
        float *a_i = new float[P]();
        float width = pow(2, limit-n->level);
        
        char idstring[100];
        n->get_idstring(idstring);
        printf("L%d%s(%d,%d)=L%d(%.1f,%.1f) ", n->level, idstring, n->x, n->y, limit, n->cx, n->cy);
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
                    printf("{q=%g@(%d,%d)=(%.1f,%.1f)} ", q, xx,yy,rx_,ry_);
                    for(unsigned int k=1; k<=P; k++) {
                        float rx__k, ry__k;
                        cmpx_pow(rx_, ry_, k, &rx__k, &ry__k);
                        a_r[k-1] += -q * rx__k;
                        a_i[k-1] += -q * ry__k;
                        // printf("{a_%d=[%g,%g]} ", k, a_r[k-1],a_i[k-1]);
                    }
                }
            }
        }
        if(Q !=0) {
            printf("Q=%g ", Q);
            printf("a_k={");
            for(unsigned int k=1; k<=P; k++)
                printf("[%g,%g] ", a_r[k-1],a_i[k-1]);
            printf("}");
        }

        float rx;
        float ry;
        // printf("interacts with ");
        for(int i=0; i<27; i++) {
            Box *ni = n->interaction[i];
            if (ni == NULL) 
                break;
            rx = ni->cx - n->cx;
            ry = ni->cy - n->cy;
            printf("[%.0f%s%.0fi] ", rx, ry>=0 ? "+" : "", ry);

            ni->potential += ni->parent->potential;
            float ln_rx, ln_ry;
            cmpx_ln(rx,ry, &ln_rx, &ln_ry);
            ni->potential += Q * ln_rx;
            // float p_real = 0, p_float = 0;
            for(unsigned int k=1; k<=P; k++) {
                float rx_k, ry_k;
                cmpx_pow(rx, ry, k, &rx_k, &ry_k);
                float div_r, div_i;
                cmpx_div(a_r[k-1], a_i[k-1], rx, ry, &div_r, &div_i);
                ni->potential += div_r / k;
            }
            printf("pot = %f ", ni->potential);
        }
        NEWLINE;
    }
    
    // Recursion
    if (n->level < limit)
        for(int i=0; i<=3; i++)
            fmm_recurse(n->child[i], limit, H, N, charge, P);
}


int fmm_calc(Box *root, const unsigned int limit) {
    const unsigned int N = (unsigned int)pow(4, limit);
    const unsigned int H = (unsigned int)sqrt(N);
    const unsigned int P = 3;   // multipole series truncation (k = 1...P)
    printf("#################\n");    
    printf("#  fmm_calc()   #\n");
    printf("#################\n");    
    printf("N = %d, log4(N) = %d, sqrt(N) = %d\n", N, limit, H);
    FILE *fh = NULL;
    
// charge matrix
    float *charge = new float[N]();
    
// place a random charge at random location
    charge[0*H + 1] = 1;
    // charge[rand_atob(0,H)*H + rand_atob(0,H)] = (float)rand_atob(1, 100);
// // place a few random charges at random locations
    // for(int yy=0; yy<H; yy++) {
        // for(int xx=0; xx<H; xx++) {
            // if(frand_atob(0, 1) < 0.1)
                // charge[yy*H + xx] = frand_atob(10.0, 100.0);
        // }
    // }
    

// display charge matrix
    printf("charge = [\n");
    fh = fopen("charge.dat", "w");
    for(unsigned int i=0; i<H; i++) {     // axis ij
        for(unsigned int j=0; j<H; j++) {
            printf("%g ", charge[i*H+j]);
            fprintf(fh, "%g ", charge[i*H+j]);
        }
        NEWLINE;
        fprintf(fh, "\n");
    }
    fclose(fh);
    printf("];\n");
    NEWLINE;

    
// do the calculations on tree
    fmm_recurse(root, limit, H, N, charge, P);


// collect result in a matrix
    float *potential = new float[N]();    // potential matrix
    traverse_tree_recurse(root, limit, H, potential);

// display potential matrix
    printf("potential = [\n");
    fh = fopen("potential.dat", "w");
    for(unsigned int i=0; i<H; i++) {     // axis ij
        for(unsigned int j=0; j<H; j++) {
            printf("%g ", potential[i*H+j]);
            fprintf(fh, "%g ", potential[i*H+j]);
        }
        NEWLINE;
        fprintf(fh, "\n");
    }
    fclose(fh);
    printf("];\n");
    NEWLINE;


    
    delete []charge ;
    return EXIT_SUCCESS;
}
