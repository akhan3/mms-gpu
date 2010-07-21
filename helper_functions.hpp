#ifndef _HELPER_FUNCTIONS_H_
#define _HELPER_FUNCTIONS_H_

int rand_atob(const int a, const int b);
float frand_atob(const float a, const float b);
float legendre(unsigned int k, float x);

// ==== Complex number functions ===========================================
// =========================================================================
float cmpx_magnitude(const float x, const float y);
void cmpx_pow(  const float x, const float y, const unsigned int p,
                float *xo, float *yo);
void cmpx_ln(   const float x, const float y,
                float *xo, float *yo);
void cmpx_div(  const float a, const float b,
                const float c, const float d,
                float *xo, float *yo);
float cmpx_costheta_between(    const float a, const float b,
                                const float c, const float d    );

int matrix2file(const float* matrix, const int rows, const int cols, const char* filename);
int matrix2stdout(const float* matrix, const int rows, const int cols, const char* matrixname);

// Depth-first tarversal
void traverse_tree_dfs(Box *n, const unsigned int limit);
// Breadth-first tarversal
void traverse_tree_bfs(Box *root, const unsigned int limit);



// prototype for FMM function
int fmm_calc(Box *root, unsigned int logN);
void fmm_bfs(   Box *n, 
                    const unsigned int limit, 
                    const unsigned int actual_limit, 
                    const unsigned int H, 
                    const float *charge, 
                    const unsigned int P    // multipole series (k = 1...P)
                );
void collect_potential_bfs( Box *n,
                            const unsigned int limit, 
                            const unsigned int actual_limit, 
                            const unsigned int H, 
                            float *potential
                          );


void create_tree_recurse(Box *thisBox, const unsigned int limit);
void find_neighbors_recurse(Box *thisBox, Box *root, const unsigned int limit);


#endif // #ifndef  _HELPER_FUNCTIONS_H_
