#ifndef _HELPER_FUNCTIONS_H_
#define _HELPER_FUNCTIONS_H_

#include "mydefs.hpp"
#include "Box.hpp"
#include "Queue.hpp"
#include "Cmpx.hpp"

int rand_atob(const int a, const int b);
fptype frand_atob(const fptype a, const fptype b);
int factorial(int x);

fptype legendre(int k, fptype x);
fptype associated_legendre(int l, int m, fptype x);
Cmpx spherical_harmonic(int l, int m, fptype theta, fptype phi);

// int matrix4mfile(const char* filename, const int rows, const int cols, int* matrix, int verbose_level);
int matrix2file(const fptype* matrix, const int rows, const int cols, const char* filename, int verbose_level);
int save_vector3d(const Vector3* vectorfield, const int zdim, const int ydim, const int xdim, const char* filename, int verbose_level);
int append_vector3d(const Vector3* vectorfield, const int zdim, const int ydim, const int xdim, const int tindex, const fptype time, FILE* fh, int verbose_level);
int save_scalar3d(const fptype* scalarfield, const int zdim, const int ydim, const int xdim, const char* filename, int verbose_level);
// int matrix2stdout(const fptype* matrix, const int rows, const int cols, const char* matrixname);

// Depth-first tarversal
void traverse_tree_dfs(Box *n, const unsigned int limit);
// Breadth-first tarversal
void traverse_tree_bfs(Box *root, const unsigned int limit);



// prototype for FMM function
int fmm_calc(   const fptype *charge,
                fptype *potential,
                const int xdim, const int ydim, const int zdim,
                const int P,    // multipole series truncation (l = 0...P)
                const int verbose_level );

void calc_potential_exact( const fptype *charge,
                        const int xdim, const int ydim, const int zdim,
                        fptype *potential);

void calc_H_nearest_neighbor(   const Vector3 *M, Vector3 *H,
                                const int xdim, const int ydim, const int zdim );

void create_tree_recurse(Box *thisBox, const unsigned int limit);
void find_neighbors_recurse(Box *thisBox, Box *root, const unsigned int limit);

#ifdef USE_FREEIMAGE
#include <FreeImage.h>
int load_mask(const char *filename, BYTE **mask, unsigned *xdim, unsigned *ydim);
#endif


#endif // #ifndef  _HELPER_FUNCTIONS_H_
