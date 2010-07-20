#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <cmath>
#include "Box.hpp"
#include "Queue.hpp"
#define NEWLINE printf("\n");

int rand_atob(const int a, const int b);
float frand_atob(const float a, const float b);

using namespace std;
using std::cout;
using std::endl;

// =========================================================================
// ==== Complex number functions ===========================================
// =========================================================================
inline float magnitude(const float x, const float y)
{
    return sqrt(x*x + y*y);
}

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


// Depth-first tarversal
// ===============================
void traverse_tree_dfs(Box *n, const unsigned int limit) 
{
    // function to perform on node
    char idstring[100];
    n->get_idstring(idstring);
    // printf("this Box is at L%d%s(%d,%d) = L%d(%.1f,%.1f)", n->level, idstring, n->x, n->y, limit, n->cx, n->cy);
    // NEWLINE;
    if (n->level < limit)
        for(int i=0; i<=3; i++)
            traverse_tree_dfs(n->child[i], limit);
}

// Breadth-first tarversal
// ===============================
void traverse_tree_bfs(Box *root, const unsigned int limit) 
{
    const unsigned int N = (unsigned int)pow(4, limit);
    Queue Q(N);
    Q.enqueue(root);
    while(!Q.isEmpty()) {
        Box *n = (Box*)Q.dequeue();
        // function to perform on node
        char idstring[100];
        n->get_idstring(idstring);
        // printf("this Box is at L%d%s(%d,%d) = L%d(%.1f,%.1f)", n->level, idstring, n->x, n->y, limit, n->cx, n->cy);
        // NEWLINE;
        // populate queue with children nodes
        if (n->level < limit)
            for(int i=0; i<=3; i++)
                Q.enqueue(n->child[i]);
    }
}

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
    const unsigned int N = (unsigned int)pow(4, limit);
    Queue Q(N);
    Q.enqueue(n);
    while(!Q.isEmpty()) {
        Box *n = (Box*)Q.dequeue();
        // function to perform on node
        // if(n->level >= 2)   // no FMM steps for Level-0 and Level-1
                            // but this condition is not necessary
        {
            float Q = 0;
            float *a = new float[P]();
            float width = pow(2, actual_limit-n->level);
            
            char idstring[100];
            n->get_idstring(idstring);
            printf("L%d%s(%d,%d)=L%d(%.1f,%.1f) ", n->level, idstring, n->x, n->y, actual_limit, n->cx, n->cy);
            printf("width=%d ", (int)width);
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
            if(0 || Q !=0) {
                printf("Q=%g ", Q);
                printf("ak={");
                for(unsigned int k=1; k<=P; k++)
                    printf("%g, ", a[k-1]);
                printf("} ");
            }

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
                    // printf("[%.0f%s%.0fi]", rx, ry>=0 ? "+" : "", ry);
                    potential_tmp += Q * log(magnitude(rx, ry));
                    for(unsigned int k=1; k<=P; k++) {
                        potential_tmp += a[k-1] / pow(magnitude(rx, ry), k) / k;
                    }
                    potential_tmp += ni->parent->potential;     // inheret potential value from the parent
                    // printf("%g>", potential_tmp);                    
                    ni->potential += potential_tmp;
                    // printf("%g ", ni->potential);
                }
            }
            NEWLINE;
        }

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
    printf("=================\n");    
    printf("|  fmm_calc()   |\n");
    printf("=================\n");    
    printf("N = %d, log4(N) = %d, sqrt(N) = %d, P = %d\n", N, limit, H, P);
    FILE *fh = NULL;

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
    
// place a random charge at random location
    int xx = H/2;
    int yy = H/2;
    charge[yy*H + xx] = 1;
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
    fmm_bfs(root, limit, limit, H, charge, P);


// collect result in a matrix
    float *potential = new float[N]();    // potential matrix
    collect_potential_bfs(root, limit, limit, H, potential);

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
