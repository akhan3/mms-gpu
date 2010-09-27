#ifndef _HELPER_FUNCTIONS_H_
#define _HELPER_FUNCTIONS_H_

#include "mydefs.hpp"
#include "Box.hpp"
#include "Queue.hpp"
#include "Cmpx.hpp"

// int matrix4mfile(const char* filename, const int rows, const int cols, int* matrix, int verbose_level);
int matrix2file(const fptype* matrix, const int rows, const int cols, const char* filename, int verbose_level);
int save_vector3d(const Vector3* vectorfield, const int zdim, const int ydim, const int xdim, const char* filename, int verbose_level);
int append_vector3d(const Vector3* vectorfield, const int zdim, const int ydim, const int xdim, FILE* fh, int verbose_level);
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
                const int use_gpu, const int verbose_level );

void calc_potential_exact( const fptype *charge,
                        const int xdim, const int ydim, const int zdim,
                        fptype *potential, int use_gpu);

int calc_potential_exact_gpu( const fptype *charge,
                        const int xdim, const int ydim, const int zdim,
                        fptype *potential);

int fmm_gpu(        const Box *const n,
                    const Cmpx *const multipole_coeff,
                    fptype *potential,
                    const unsigned int limit,
                    const int P,    // multipole series truncation (l = 0...P)
                    const int xdim, const int ydim, const int zdim,
                    const int zc,   // charge layer
                    const int use_gpu,
                    const int verbose_level
                );

void calc_H_nearest_neighbor(   const Vector3 *M, Vector3 *H,
                                const int xdim, const int ydim, const int zdim );

#ifdef USE_FREEIMAGE
#include <FreeImage.h>
int load_mask(const char *filename, BYTE **mask, unsigned *xdim, unsigned *ydim);
#endif


#endif // #ifndef  _HELPER_FUNCTIONS_H_
