#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cmath>
#include "Vector3.hpp"
#include "helper_functions.hpp"

// // Read matrix from file
// // ===============================
// int matrix4mfile(const char* filename, const int rows, const int cols, int* matrix, int verbosity)
// {
    // int status = 0;
    // FILE *fh = fopen(filename, "r");
    // if(fh == NULL) {
        // printf("FATAL ERROR: Error opening file %s\n", filename);
        // return EXIT_FAILURE;
    // }
    // int dummy;
    // for(int r=0; r<rows; r++) {     // axis ij
    // // for(int r=rows-1; r>=0; r--) {     // axis xy
        // for(int c=0; c<cols; c++)
            // dummy = fscanf(fh, "%d ", &matrix[r*cols+c]);
        // dummy = fscanf(fh, "\n");
    // }
    // fclose(fh);
    // if(verbosity >= 8)
        // printf("INFO: Read file %s with status=%d\n", filename, status);
    // return status ? EXIT_FAILURE : EXIT_SUCCESS;
// }

// Write matrix to file
// ===============================
int matrix2file(const fptype* matrix, const int rows, const int cols, const char* filename, int verbosity)
{
    int status = 0;
    FILE *fh = fopen(filename, "w");
    if(fh == NULL) {
        printf("FATAL ERROR: Error opening file %s\n", filename);
        return EXIT_FAILURE;
    }
    for(int r=0; r<rows; r++) {     // axis ij
    // for(int r=rows-1; r>=0; r--) {     // axis xy
        for(int c=0; c<cols; c++)
            fprintf(fh, "%g ", matrix[r*cols+c]);
        fprintf(fh, "\n");
    }
    fclose(fh);
    if(verbosity >= 8)
        printf("INFO: Written file %s with status=%d\n", filename, status);
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}

// Write 3D scalar field to file
// ===============================
int save_scalar3d(const fptype* scalarfield, const int zdim, const int ydim, const int xdim, const char* filename, int verbosity)
{
    int status = 0;
    FILE *fh = fopen(filename, "w");
    if(fh == NULL) {
        fprintf(stderr, "%s:%d FATAL ERROR: couldn't open file %s\n", __FILE__, __LINE__, filename);
        return EXIT_FAILURE;
    }
    for(int z = 0; z < zdim; z++) {
        for(int x = 0; x < xdim; x++) {     // mind it!! writing in column major order for MATLAB
            for(int y = 0; y < ydim; y++) {
                fprintf(fh, "%g ", scalarfield[z*ydim*xdim + y*xdim + x]);
            }
            fprintf(fh, "\n");
        }
        fprintf(fh, "\n");
    }
    fclose(fh);
    if(verbosity >= 8)
        printf("INFO: Written file %s with status=%d\n", filename, status);
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}

// Write 3D vector field to file
// ===============================
int save_vector3d(const Vector3* vectorfield, const int zdim, const int ydim, const int xdim, const char* filename, int verbosity)
{
    int status = 0;
    FILE *fh = fopen(filename, "w");
    if(fh == NULL) {
        fprintf(stderr, "%s:%d FATAL ERROR: couldn't open file %s\n", __FILE__, __LINE__, filename);
        return EXIT_FAILURE;
    }
    for(int z = 0; z < zdim; z++) {
        for(int x = 0; x < xdim; x++) {     // mind it!! writing in column major order for MATLAB
            for(int y = 0; y < ydim; y++) {
                fprintf(fh, "%g %g %g \n", vectorfield[z*ydim*xdim + y*xdim + x].x, vectorfield[z*ydim*xdim + y*xdim + x].y, vectorfield[z*ydim*xdim + y*xdim + x].z);
            }
            // fprintf(fh, "# y done\n");
        }
        // fprintf(fh, "# yx done\n");
    }
    // fprintf(fh, "# yxz done\n");
    fclose(fh);
    if(verbosity >= 8) {
        printf("INFO: Written file %s with status=%d\n", filename, status);
    }
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}

// Append 3D vector field to file
// ===============================
int append_vector3d(const Vector3* vectorfield, const int zdim, const int ydim, const int xdim, FILE* fh, int verbosity)
{
    if(verbosity) {}
    int status = 0;
    for(int z = 0; z < zdim; z++) {
        for(int x = 0; x < xdim; x++) {     // mind it!! writing in column major order for MATLAB
            for(int y = 0; y < ydim; y++) {
                fprintf(fh, "%g %g %g \n", vectorfield[z*ydim*xdim + y*xdim + x].x, vectorfield[z*ydim*xdim + y*xdim + x].y, vectorfield[z*ydim*xdim + y*xdim + x].z);
            }
            // fprintf(fh, "# y done\n");
        }
        // fprintf(fh, "# yx done\n");
    }
    // fprintf(fh, "# yxz done\n");
    // fprintf(fh, "# tindex=%d, time=%g done\n", tindex, time);
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}

int matrix2stdout(const fptype* matrix, const int rows, const int cols, const char* matrixname) {
    printf("%s = [\n", matrixname);
    for(int r=0; r<rows; r++) {     // axis ij
    // for(int r=rows-1; r>=0; r--) {     // axis xy
        for(int c=0; c<cols; c++)
            printf("%g ", matrix[r*cols+c]);
        printf("\n");
    }
    printf("];\n");
    return EXIT_SUCCESS;
}

// // Depth-first tarversal
// // ===============================
// void traverse_tree_dfs(Box *n, const unsigned int limit)
// {
    // // function to perform on node
    // char idstring[100];
    // n->get_idstring(idstring);
    // // printf("this Box is at L%d%s(%d,%d) = L%d(%.1f,%.1f)", n->level, idstring, n->x, n->y, limit, n->cx, n->cy);
    // // NEWLINE;
    // if(n->level < limit)
        // for(int i=0; i<=3; i++)
            // traverse_tree_dfs(n->child[i], limit);
// }

// // Breadth-first tarversal
// // ===============================
// void traverse_tree_bfs(Box *root, const unsigned int limit)
// {
    // const unsigned int N = (unsigned int)pow(4, limit);
    // Queue Q(N);
    // Q.enqueue(root);
    // while(!Q.isEmpty()) {
        // Box *n = (Box*)Q.dequeue();
        // if(n->level < limit)
            // for(int i=0; i<=3; i++)
                // Q.enqueue(n->child[i]);
        // // function to perform on node
        // char idstring[100];
        // n->get_idstring(idstring);
        // // printf("this Box is at L%d%s(%d,%d) = L%d(%.1f,%.1f)", n->level, idstring, n->x, n->y, limit, n->cx, n->cy);
        // // NEWLINE;
        // // populate queue with children nodes
    // }
// }


#ifdef USE_FREEIMAGE
int load_mask(const char *filename, BYTE **mask, unsigned *xdim, unsigned *ydim)
{
    assert( !strcmp(filename+strlen(filename)-3, "png") || !strcmp(filename+strlen(filename)-3, "PNG") );
    FIBITMAP *myimage = FreeImage_Load(FIF_PNG, filename, PNG_DEFAULT);
    assert(FreeImage_GetWidth(myimage) == FreeImage_GetHeight(myimage));
    assert(FreeImage_GetColorType(myimage) != FIC_RGBALPHA);

    // std::cout << "type = "       << FreeImage_GetImageType(myimage) << std::endl;
    // std::cout << "#colors = "    << FreeImage_GetColorsUsed(myimage) << std::endl;
    // std::cout << "bpp = "        << FreeImage_GetBPP(myimage) << std::endl;
    // std::cout << "width = "        << FreeImage_GetWidth(myimage) << std::endl;
    // std::cout << "height = "        << FreeImage_GetHeight(myimage) << std::endl;
    // std::cout << "color type = "        << FreeImage_GetColorType(myimage) << std::endl;
    // std::cout << "red mask = "        << FreeImage_GetRedMask(myimage) << std::endl;
    // std::cout << "green mask = "        << FreeImage_GetGreenMask(myimage) << std::endl;
    // std::cout << "blue mask = "        << FreeImage_GetBlueMask(myimage) << std::endl;
    // std::cout << "is transparent = "        << FreeImage_IsTransparent(myimage) << std::endl;
    // std::cout << "file type = "        << FreeImage_GetFileType(filename) << std::endl;

    *xdim = FreeImage_GetWidth(myimage);
    *ydim = FreeImage_GetHeight(myimage);
    *mask = new BYTE[FreeImage_GetHeight(myimage) * FreeImage_GetWidth(myimage)]();
    if(mask == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
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
#endif
