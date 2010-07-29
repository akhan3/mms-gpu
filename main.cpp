#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <cmath>
#include <FreeImage.h>
#include "Box.hpp"
#include "helper_functions.hpp"
#define NEWLINE printf("\n");


//*************************************************************************//
//******************** Main function **************************************//
//*************************************************************************//
int main(int argc, char **argv)
{
    int status = 0;
    unsigned int seed = time(NULL);
    unsigned int N1 = 256*256;
    unsigned int P = 3;
    if(argc >= 2) {
        sscanf(argv[1], "%u", &N1);
        assert(N1 > 0);
    }
    if(argc >= 3) {
        sscanf(argv[2], "%u", &seed);
        assert(seed > 0);
    }
    if(argc >= 4) {
        sscanf(argv[3], "%u", &P);
        assert(seed > 0);
    }
    assert(P <= 4);
    srand(seed);
    const unsigned int logN = ceil(log2(N1) / log2(4));
    const unsigned int N = (unsigned int)pow(4, logN);
    const unsigned int h = (unsigned int)sqrt(N);
    printf("N = %d, log4(N) = %d\n", N, logN);
    // printf("sizeof(Box) = %ld\n", sizeof(Box));



    Box *root = new Box(0, logN);

// generate the tree
    create_tree_recurse(root, logN);
    find_neighbors_recurse(root, root, logN);



// Material parameters
// ================================================
    const float Ms = 8.6e5; // permalloy
    const float meshwidth = 1e-9;
    const float meshdepth = 1e-9; // depth of material in z-dimension

// Mask configuration for magnetic material
// ==========================================
    int *mask = new int[N]();    // mask matrix
// read the mask from file
    char filename[100];
    sprintf(filename, "barmagnet_60x90_%dx%d_mask.data", h,h);
    status |= matrix4mfile(filename, h, h, mask);
    if(status) return EXIT_FAILURE;
// dipole
    // int l = h/8;
    // int x0 = h/2;
    // int y1 = (h+l)/2;
    // int y2 = (h-l)/2;
    // mask[y1*h + x0] = 1;
    // mask[y2*h + x0] = 1;
// place a charge at the center
    // int xx = h/2;
    // int yy = h/2;
    // mask[yy*h + xx] = 1;
    // mask[(yy-1)*h + xx] = 1;
    // mask[yy*h + xx-1] = 1;
    // mask[(yy-1)*h + xx-1] = 1;
// Few random charges at random locations
    // const float mask_probability = .1;
    // for(unsigned int yy=0; yy<h; yy++)
        // for(unsigned int xx=0; xx<h; xx++)
            // if(frand_atob(0, 1) < mask_probability)
                // mask[yy*h + xx] = 1;


// generate the initial magnetization distribution
    Cmpx *M = new Cmpx[N]();    // magnetization matrix
    for(unsigned int y = 0; y < h; y++)
        for(unsigned int x = 0; x < h; x++)
            if (!mask[y*h + x])
                M[y*h + x].init(0, Ms, 0);

// magnetic volume charge density
    float *charge = new float[N]();
    divergence_2d(M, h, h, 1, charge);

// write charge matrix to file
    status |= matrix2file(charge, h, h, "charge.dat");
    if(status) return EXIT_FAILURE;

// magnetic scalar potential
    float *potential = new float[N]();

// call the function for main FMM calculation
    status |= fmm_calc(root, logN, charge, potential);
    if(status) return EXIT_FAILURE;

// write potential matrix to file
    status |= matrix2file(potential, h, h, "potential.dat");
    if(status) return EXIT_FAILURE;

// magnetic field
    Cmpx *H = new Cmpx[N]();    // magnetic field matrix
    const float constant_multiple = (meshdepth / meshwidth) / (4 * M_PI);
    printf("constant_multiple = %f\n", constant_multiple);
    gradient_2d(potential, h, h, 1, H);
    for(unsigned int i = 0; i < N; i++)
        H[i] *= constant_multiple;

// write H matrix to file
    float *Hx = new float[N]();
    float *Hy = new float[N]();
    for(unsigned int i = 0; i < N; i++) {
        Hx[i] = H[i].get_re();
        Hy[i] = H[i].get_im();
    }
    status |= matrix2file(Hx, h, h, "Hx.dat");
    status |= matrix2file(Hy, h, h, "Hy.dat");
    if(status) return EXIT_FAILURE;


// closing
    delete root;
    delete []potential;
    delete []Hx;
    delete []Hy;
    delete []H;
    delete []charge;

    printf("SEED = %d\n", seed);
    return EXIT_SUCCESS;
}
