#ifndef _HELPER_FUNCTIONS_H_
#define _HELPER_FUNCTIONS_H_

#include "Box.hpp"
#include "Queue.hpp"
#include "Cmpx.hpp"

int rand_atob(const int a, const int b);
float frand_atob(const float a, const float b);
int factorial(int x);

float legendre(int k, float x);
float associated_legendre(int l, int m, float x);
Cmpx spherical_harmonic(int l, int m, float theta, float phi);

int matrix4mfile(const char* filename, const int rows, const int cols, int* matrix);
int matrix2file(const float* matrix, const int rows, const int cols, const char* filename);
int matrix2stdout(const float* matrix, const int rows, const int cols, const char* matrixname);

// Depth-first tarversal
void traverse_tree_dfs(Box *n, const unsigned int limit);
// Breadth-first tarversal
void traverse_tree_bfs(Box *root, const unsigned int limit);



// prototype for FMM function
int fmm_calc(const Box *root, const unsigned int limit, const float *charge, float *potential);

void create_tree_recurse(Box *thisBox, const unsigned int limit);
void find_neighbors_recurse(Box *thisBox, Box *root, const unsigned int limit);

void divergence_2d( const Cmpx *V,
                    const int xdim, const int ydim, const float meshwidth,
                    float *S );
void gradient_2d(   const float *S,
                    int xdim, int ydim, float meshwidth,
                    Cmpx *V );

#endif // #ifndef  _HELPER_FUNCTIONS_H_
