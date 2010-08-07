#ifndef _HELPER_FUNCTIONS_H_
#define _HELPER_FUNCTIONS_H_

#include <FreeImage.h>
#include "Box.hpp"
#include "Queue.hpp"
#include "Cmpx.hpp"

int rand_atob(const int a, const int b);
float frand_atob(const float a, const float b);
int factorial(int x);

float legendre(int k, float x);
float associated_legendre(int l, int m, float x);
Cmpx spherical_harmonic(int l, int m, float theta, float phi);

// int matrix4mfile(const char* filename, const int rows, const int cols, int* matrix, int verbose_level);
int matrix2file(const float* matrix, const int rows, const int cols, const char* filename, int verbose_level);
int save_vector3d(const Vector3* vectorfield, const int zdim, const int ydim, const int xdim, const char* filename, int verbose_level);
int save_scalar3d(const float* scalarfield, const int zdim, const int ydim, const int xdim, const char* filename, int verbose_level);
// int matrix2stdout(const float* matrix, const int rows, const int cols, const char* matrixname);

// Depth-first tarversal
void traverse_tree_dfs(Box *n, const unsigned int limit);
// Breadth-first tarversal
void traverse_tree_bfs(Box *root, const unsigned int limit);



// prototype for FMM function
int fmm_calc(   const float *charge,
                float *potential,
                const int xdim, const int ydim, const int zdim,
                const int P,    // multipole series truncation (l = 0...P)
                const int verbose_level );

void calc_potential_exact( const float *charge,
                        const int xdim, const int ydim, const int zdim,
                        float *potential);

void calc_H_nearest_neighbor(   const Vector3 *M, Vector3 *H,
                                const int xdim, const int ydim, const int zdim );

void create_tree_recurse(Box *thisBox, const unsigned int limit);
void find_neighbors_recurse(Box *thisBox, Box *root, const unsigned int limit);

int load_mask(const char *filename, BYTE **mask, unsigned *xdim, unsigned *ydim);

#endif // #ifndef  _HELPER_FUNCTIONS_H_
