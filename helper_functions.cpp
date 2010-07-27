#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cmath>
#include "helper_functions.hpp"

#define NEWLINE printf("\n");


int factorial(int x)
{
    // int fac = 1;
    // for(int i = 2; i <= x; i++)
        // fac *= i;
    // return fac;
    // return (x == 0) ? 1 : x*factorial(x-1);
    // assert(x >= 0 && x <= 8);
    switch (x) {
        case 0: return     1;
        case 1: return     1;
        case 2: return     2;
        case 3: return     6;
        case 4: return    24;
        case 5: return   120;
        case 6: return   720;
        case 7: return  5040;
        case 8: return 40320;
        default:
            return (x == 0) ? 1 : x*factorial(x-1);
    }
}

// Legendre function (recursive implementation)
float legendre(int k, float x)
{
    assert(k >= 0);
    assert(x >= -1 && x <= 1);
    switch (k) {
        case 0:
            return 1;
        case 1:
            return x;
        case 2:
            return 1/2.0 * (3 * x*x - 1);
        case 3:
            return 1/2.0 * (5 * x*x*x - 3 * x);
        default:
            return ((2*k-1) * x * legendre(k-1, x) - (k-1) * legendre(k-2, x)) / k;
    }
}

// Associated Legendre function (Lookup table for x = 0)
float associated_legendre(int l, int m, float x)
{
    assert(l >= 0 && l <= 4);
    assert(abs(m) <= l);
    assert(m == abs(m));
    assert(fabs(x) <= 1e-7);    // assert(fabs(x) == 0.0);
    switch (l) {
        case 0:
            switch (m) {
                case 0:  return 1;
                default: return 0;
            }
        case 1:
            switch (m) {
                case 0:  return 0;
                case 1:  return -1;
                default: return 0;
            }
        case 2:
            switch (m) {
                case 0:  return -1/2.0;
                case 1:  return 0;
                case 2:  return 3;
                default: return 0;
            }
        case 3:
            switch (m) {
                case 0:  return 0;
                case 1:  return 3/2.0;
                case 2:  return 0;
                case 3:  return -15;
                default: return 0;
            }
        case 4:
            switch (m) {
                case 0:  return 3/8.0;
                case 1:  return 0;
                case 2:  return -15/2.0;
                case 3:  return 0;
                case 4:  return 105;
                default: return 0;
            }
        default:
            printf("FATAL ERROR؛ associated_legendre(l=%d, m=%d, x=%f) is not implemented\n", l,m,x);
            return 0.0/0.0; // NaN
    }
}

// Spherical harmonics
Cmpx spherical_harmonic(int l, int m, float theta, float phi)
{
    // return Cmpx(sqrt((1.0*factorial(l-abs(m))) / factorial(l+abs(m))) * associated_legendre(l,abs(m),cos(theta)), m*phi, 1);
    return Cmpx(associated_legendre(l,abs(m),cos(theta)), m*phi, 1);
}

// Write matrix to file
// ===============================
int matrix2file(const float* matrix, const int rows, const int cols, const char* filename) {
    FILE *fh = fopen(filename, "w");
    if(fh == NULL) {
        printf("FATAL ERROR: Erro opening file %s\n", filename);
        return EXIT_FAILURE;
    }
    for(int r=0; r<rows; r++) {     // axis ij
        for(int c=0; c<cols; c++)
            fprintf(fh, "%g ", matrix[r*cols+c]);
        fprintf(fh, "\n");
    }
    fclose(fh);
    return EXIT_SUCCESS;
}

int matrix2stdout(const float* matrix, const int rows, const int cols, const char* matrixname) {
    printf("%s = [\n", matrixname);
    for(int r=0; r<rows; r++) {     // axis ij
        for(int c=0; c<cols; c++)
            printf("%g ", matrix[r*cols+c]);
        printf("\n");
    }
    printf("];\n");
    return EXIT_SUCCESS;
}

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
        if (n->level < limit)
            for(int i=0; i<=3; i++)
                Q.enqueue(n->child[i]);
        // function to perform on node
        char idstring[100];
        n->get_idstring(idstring);
        // printf("this Box is at L%d%s(%d,%d) = L%d(%.1f,%.1f)", n->level, idstring, n->x, n->y, limit, n->cx, n->cy);
        // NEWLINE;
        // populate queue with children nodes
    }
}


void create_tree_recurse(Box *thisBox, const unsigned int limit) {
    // Recursion
    if (thisBox->level < limit) {
        // function to perform on node
        thisBox->split(limit);
        for(int i=0; i<=3; i++)
            create_tree_recurse(thisBox->child[i], limit);
    }
}

void find_neighbors_recurse(Box *thisBox, Box *root, const unsigned int limit) {
    // function to perform on node
    thisBox->find_neighbors(root);
    // Recursion
    if (thisBox->level < limit) {
        for(int i=0; i<=3; i++)
            find_neighbors_recurse(thisBox->child[i], root, limit);
    }
}

int rand_atob(const int a, const int b) {
    double r = rand() / (double)RAND_MAX;
    r = a + (b-a) * r;
    return (int)r;
}

float frand_atob(const float a, const float b) {
    double r = rand() / (double)RAND_MAX;
    r = a + (b-a) * r;
    return (float)r;
}

// void cmpx_pow(  const float x, const float y, const unsigned int p,
                // float *xo, float *yo)
// {
    // *xo = 1;
    // *yo = 0;
    // for(unsigned int i=1; i<=p; i++) {
        // float x1;
        // float y1;
        // x1 = x*(*xo) - y*(*yo);
        // y1 = y*(*xo) + x*(*yo);
        // *xo = x1;
        // *yo = y1;
    // }
// }

// void cmpx_ln(  const float x, const float y,
                // float *xo, float *yo)
// {
    // float r = sqrt(x*x + y*y);
    // float theta = atan2(y, x);
    // *xo = log(r);
    // *yo = theta;
// }
