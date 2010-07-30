#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <assert.h>
#include <iostream>
#include <cmath>
#include <FreeImage.h>
#include "Box.hpp"
#include "helper_functions.hpp"
#define NEWLINE printf("\n");
using std::cout;
using std::endl;


int load_mask(const char *filename, BYTE **mask, unsigned *dimx, unsigned *dimy) {
    assert( !strcmp(filename+strlen(filename)-3, "png") || !strcmp(filename+strlen(filename)-3, "PNG") );
    FIBITMAP *myimage = FreeImage_Load(FIF_PNG, filename, PNG_DEFAULT);
    assert(FreeImage_GetWidth(myimage) == FreeImage_GetHeight(myimage));
    assert(FreeImage_GetColorType(myimage) != FIC_RGBALPHA);

    // cout << "type = "       << FreeImage_GetImageType(myimage) << endl;
    // cout << "#colors = "    << FreeImage_GetColorsUsed(myimage) << endl;
    // cout << "bpp = "        << FreeImage_GetBPP(myimage) << endl;
    cout << "width = "        << FreeImage_GetWidth(myimage) << endl;
    cout << "height = "        << FreeImage_GetHeight(myimage) << endl;
    // cout << "color type = "        << FreeImage_GetColorType(myimage) << endl;
    // cout << "red mask = "        << FreeImage_GetRedMask(myimage) << endl;
    // cout << "green mask = "        << FreeImage_GetGreenMask(myimage) << endl;
    // cout << "blue mask = "        << FreeImage_GetBlueMask(myimage) << endl;
    // cout << "is transparent = "        << FreeImage_IsTransparent(myimage) << endl;
    // cout << "file type = "        << FreeImage_GetFileType(filename) << endl;

    *dimx = FreeImage_GetWidth(myimage);
    *dimy = FreeImage_GetHeight(myimage);
    *mask = new BYTE[FreeImage_GetHeight(myimage) * FreeImage_GetWidth(myimage)]();
    // Calculate the number of bytes per pixel (3 for 24-bit or 4 for 32-bit)
    int bytespp = FreeImage_GetLine(myimage) / FreeImage_GetWidth(myimage);
    for(unsigned y = 0; y < FreeImage_GetHeight(myimage); y++) {
        BYTE *bits = FreeImage_GetScanLine(myimage, y);
        for(unsigned x = 0; x < FreeImage_GetWidth(myimage); x++) {
            // read pixel color
            BYTE r = bits[FI_RGBA_RED];
            BYTE g = bits[FI_RGBA_GREEN];
            BYTE b = bits[FI_RGBA_BLUE];
            BYTE gray = (BYTE)(.3*r + .59*g + .11*b);
            (*mask)[y*FreeImage_GetWidth(myimage) + x] = gray;
            // jump to next pixel
            bits += bytespp;
        }
    }

    FreeImage_Unload(myimage);
    return EXIT_SUCCESS;
}

//*************************************************************************//
//******************** Main function **************************************//
//*************************************************************************//
int main(int argc, char **argv)
{
    int status = 0;
    char filename_arg[1000];
    unsigned int seed = time(NULL);
    unsigned int P = 3;
    assert(argc >= 2);
    if(argc >= 2) {
        sscanf(argv[1], "%s", filename_arg);
    }
    if(argc >= 3) {
        sscanf(argv[2], "%u", &P);
        assert(P <= 4);
    }
    if(argc >= 4) {
        sscanf(argv[3], "%u", &seed);
    }
    srand(seed);


// Material parameters
// ================================================
    const float Ms = 8.6e5; // permalloy
    const float meshwidth = 1e-9;
    const float meshdepth = 1e-9; // depth of material in z-dimension

// Mask configuration for magnetic material
    BYTE *mask = NULL; // mask matrix
    unsigned dimx = 0, dimy = 0;
    char filename[1000];
    sprintf(filename, "%s", filename_arg);
// read the mask from file
    load_mask(filename, &mask, &dimx, &dimy);
    assert(dimx == dimy);

// set up grid
    const unsigned int N = dimx * dimy;
    const unsigned int logN = ceil(log2(N) / log2(4));
    const unsigned int h = (unsigned int)sqrt(N);
    printf("N = %d, log4(N) = %d\n", N, logN);

// generate the tree
    Box *root = new Box(0, logN);
    create_tree_recurse(root, logN);
    find_neighbors_recurse(root, root, logN);

// generate the initial magnetization distribution
    Cmpx *M = new Cmpx[N]();    // magnetization matrix
    for(unsigned int y = 0; y < h; y++)
        for(unsigned int x = 0; x < h; x++)
            if (!mask[y*h + x])
                M[y*h + x].init(Ms, 1*M_PI/2, 1);
                // M[y*h + x].init(Ms, frand_atob(0, 2*M_PI), 1);

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
    delete []mask;

    printf("SEED = %d\n", seed);
    return EXIT_SUCCESS;
}
