#ifndef _HELPER_FUNCTIONS_H_
#define _HELPER_FUNCTIONS_H_

#include "mydefs.hpp"
#include "Box.hpp"
#include "Queue.hpp"
#include "Cmpx.hpp"

// int matrix4mfile(const char* filename, const int rows, const int cols, int* matrix, int verbosity);
int matrix2file(const fptype* matrix, const int rows, const int cols, const char* filename, int verbosity);
int save_vector3d(const Vector3* vectorfield, const int zdim, const int ydim, const int xdim, const char* filename, int verbosity);
int append_vector3d(const Vector3* vectorfield, const int zdim, const int ydim, const int xdim, FILE* fh, int verbosity);
int save_scalar3d(const fptype* scalarfield, const int zdim, const int ydim, const int xdim, const char* filename, int verbosity);
// int matrix2stdout(const fptype* matrix, const int rows, const int cols, const char* matrixname);

// Depth-first tarversal
void traverse_tree_dfs(Box *n, const unsigned int limit);
// Breadth-first tarversal
void traverse_tree_bfs(Box *root, const unsigned int limit);



// prototype for FMM function
int fmm_calc(   const fptype *charge,
                fptype *potential,
                const int xdim, const int ydim, const int zdim,
                const fptype dx, const fptype dy, const fptype dz,
                const int P,    // multipole series truncation (l = 0...P)
                const int use_gpu, const int verbosity );

int calc_potential_exact( const fptype *charge,
                        const int xdim, const int ydim, const int zdim,
                        const fptype dx, const fptype dy, const fptype dz,
                        fptype *potential, int use_gpu,
                        const int verbosity);

int calc_potential_exact_gpu( const fptype *charge,
                        const int xdim, const int ydim, const int zdim,
                        const fptype dx, const fptype dy, const fptype dz,
                        fptype *potential);

int fmm_gpu(        const Box *const n,
                    const Cmpx *const mpc_gmem,
                    fptype *potential_gmem,
                    const unsigned int limit,
                    const int P,    // multipole series truncation (l = 0...P)
                    const int xdim, const int ydim, const int zdim,
                    const fptype dx, const fptype dy, const fptype dz,
                    const int zc,   // charge layer
                    const int use_gpu,
                    const int verbosity
                );

void calc_H_nearest_neighbor(   const Vector3 *M, Vector3 *H,
                                const int xdim, const int ydim, const int zdim );

#ifdef USE_FREEIMAGE
#include <FreeImage.h>
int load_mask(const char *filename, BYTE **mask, int *xdim, int *ydim);
#endif


#endif // #ifndef  _HELPER_FUNCTIONS_H_
