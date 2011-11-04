#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cmath>
#include "Vector3.hpp"
#include "helper_functions.hpp"

// Read 3d vector field from file
// ===============================
int load_Minit(const char* filename, int *rows, int *cols, Vector3 **M, int verbosity)
{
    int status = 0;
    FILE *fh = fopen(filename, "r");
    // FILE *fw = fopen("debugFile.txt", "w");
    if(fh == NULL) {
        printf("FATAL ERROR: Error opening file %s\n", filename);
        return EXIT_FAILURE;
    }
    int dummy;
    // read out Nx and Ny
    dummy = fscanf(fh, "Nx=%d\n", cols);
    dummy = fscanf(fh, "Ny=%d\n", rows);
    dummy = fscanf(fh, "# start field\n");
    // allocate M of size Nx*Ny
    *M = new Vector3[*rows * *cols]();  // magnetization matrix
    if(*M == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }

    Vector3 m;
    for(int r=0; r<*rows; r++) {     // axis ij
    // for(int r=rows-1; r>=0; r--) {   // axis xy
        for(int c=0; c<*cols; c++) {
            dummy = fscanf(fh, "%g %g %g", &m.x, &m.y, &m.z);
            (*M)[r* *cols+c] = m;
            if(dummy != 3) {
                printf("FATAL ERROR: Malformed file %s @ (%d,%d)\n", filename,r+1,c+1);
                return EXIT_FAILURE;
            }
        }
        dummy = fscanf(fh, "\n");
        // fprintf(fw, "\n");
    }
    fclose(fh);
    // fclose(fw);
    if(verbosity >= 8)
        printf("INFO: Read file %s with status=%d\n", filename, status);
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}

// Read 2D scalar field from file
// ===============================
int load_scalarField2D(const char* filename, int *rows, int *cols, fptype **scalar, int verbosity)
{
    int status = 0;
    FILE *fh = fopen(filename, "r");
    if(fh == NULL) {
        printf("FATAL ERROR: Error opening file %s\n", filename);
        return EXIT_FAILURE;
    }
    int dummy;
    // read out Nx and Ny
    dummy = fscanf(fh, "Nx=%d\n", cols);
    dummy = fscanf(fh, "Ny=%d\n", rows);
    dummy = fscanf(fh, "# start field\n");
    // allocate s of size Nx*Ny
    *scalar = new fptype[*rows * *cols]();  // scalar field matrix
    if(*scalar == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return EXIT_FAILURE;
    }
    for(int r=0; r<*rows; r++) {     // axis ij
        for(int c=0; c<*cols; c++) {
            dummy = fscanf(fh, "%g ", &(*scalar)[r* *cols+c]);
            if(dummy != 1) {
                printf("FATAL ERROR: Malformed file %s @ (%d,%d)\n", filename,r+1,c+1);
                return EXIT_FAILURE;
            }
        }
        dummy = fscanf(fh, "\n");
    }
    fclose(fh);
    if(verbosity >= 8)
        printf("INFO: Read file %s with status=%d\n", filename, status);
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
}

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
    // write Nx and Ny
    fprintf(fh, "Nx=%d\n", cols);
    fprintf(fh, "Ny=%d\n", rows);
    fprintf(fh, "# start field\n");
    for(int r=0; r<rows; r++) {     // axis ij
    // for(int r=rows-1; r>=0; r--) {     // axis xy
        for(int c=0; c<cols; c++)
            fprintf(fh, "%g ", matrix[r*cols+c]);
        fprintf(fh, "\n");
    }
    fprintf(fh, "# end field\n");
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
    // for(int z = 0; z < zdim; z++) {
        // for(int x = 0; x < xdim; x++) {     // mind it!! writing in column major order for MATLAB
            // for(int y = 0; y < ydim; y++) {
    for(int z = 0; z < zdim; z++) {
        for(int y = 0; y < ydim; y++) {
            for(int x = 0; x < xdim; x++) {     // mind it!! writing in row major order for C
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
    // write Nx and Ny
    fprintf(fh, "Nx=%d\n", xdim);
    fprintf(fh, "Ny=%d\n", ydim);
    fprintf(fh, "# start field\n");
    // for(int z = 0; z < zdim; z++) {
        // for(int x = 0; x < xdim; x++) {     // mind it!! writing in column major order for MATLAB
            // for(int y = 0; y < ydim; y++) {
    for(int z = 0; z < zdim; z++) {
        for(int y = 0; y < ydim; y++) {
            for(int x = 0; x < xdim; x++) {     // mind it!! writing in row major order for C
                fprintf(fh, "%g %g %g \n", vectorfield[z*ydim*xdim + y*xdim + x].x, vectorfield[z*ydim*xdim + y*xdim + x].y, vectorfield[z*ydim*xdim + y*xdim + x].z);
            }
            fprintf(fh, "\n");
            // fprintf(fh, "# y done\n");
        }
        fprintf(fh, "\n\n");
    }
    fprintf(fh, "# end field\n");
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
    // for(int z = 0; z < zdim; z++) {
        // for(int x = 0; x < xdim; x++) {     // mind it!! writing in column major order for MATLAB
            // for(int y = 0; y < ydim; y++) {
    for(int z = 0; z < zdim; z++) {
        for(int y = 0; y < ydim; y++) {
            for(int x = 0; x < xdim; x++) {     // mind it!! writing in row major order for C
                fprintf(fh, "%g %g %g \n", vectorfield[z*ydim*xdim + y*xdim + x].x, vectorfield[z*ydim*xdim + y*xdim + x].y, vectorfield[z*ydim*xdim + y*xdim + x].z);
            }
            fprintf(fh, "\n");
            // fprintf(fh, "# y done\n");
        }
        // fprintf(fh, "# yx done\n");
        fprintf(fh, "\n\n");
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

// Append 3D vector field to file
// ===============================
int append_vector3d_midrow(const Vector3* vectorfield, const int xdim, FILE* fh) {
    int status = 0;
    for(int x = 0; x < xdim; x++) {
        fprintf(fh, "%g %g %g \n", vectorfield[x].x, vectorfield[x].y, vectorfield[x].z);
    }
    fprintf(fh, "\n");
    return status ? EXIT_FAILURE : EXIT_SUCCESS;
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
int load_mask(const char *filename, BYTE **mask, int *xdim, int *ydim)
{
    assert( !strcmp(filename+strlen(filename)-3, "png") || !strcmp(filename+strlen(filename)-3, "PNG") );
    FIBITMAP *myimage = FreeImage_Load(FIF_PNG, filename, PNG_DEFAULT);
    // assert(FreeImage_GetWidth(myimage) == FreeImage_GetHeight(myimage));
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

/*
int load_mask_txt(const char *filename, byte **mask, int *xdim, int *ydim) {
    // mask = new byte[ydim*xdim](); // mask matrix
    FILE *fh = fopen(filename, "r");
    if(fh == NULL) {
        printf("FATAL ERROR: Error opening file %s\n", filename);
        return EXIT_FAILURE;
    }
    // obtain file size:
    fseek(fh, 0, SEEK_END);
    long lSize = ftell(fh);
    rewind(fh);
    char *buffer = new char[lSize]();
    char *bufferCopy = new char[lSize]();
    fread (buffer, 1, lSize, fh);
    strcpy(bufferCopy, buffer);
    // printf("lSize=%d buffer = \n%s\n", lSize, buffer);
    // printf("ASCII\n");
    // for(int i = 0; i < strlen(buffer); i++) {
        // if(buffer[i] == '\n') {
            // printf("%d [%3d '\\n']\n", i, buffer[i]);
            // break;
        // }
        // else
            // printf("%d [%3d '%c']\n", i, buffer[i], buffer[i]);
    // }
    int n = strchr(buffer, 10)-buffer;
    // printf("Newline first found at %d\n", n);
    // find xdim
    char *pch;
    pch = strtok (buffer, " ,\t");
    *xdim = 0;
    while (pch != NULL && pch < buffer+n)
    {
        (*xdim)++;
        // printf ("<%s>", pch);
        pch = strtok (NULL, " ,\t");
    }
    // printf("xdim = %d\n", *xdim);
    // find ydim
    strcpy(buffer, bufferCopy);
    pch = strtok (bufferCopy, "\n");
    *ydim = 0;
    while (pch != NULL)
    {
        (*ydim)++;
        // printf ("<%s>\n", pch);
        pch = strtok (NULL, "\n");
    }
    // printf("ydim = %d\n", *ydim);
    // load mask
    *mask = new byte[(*ydim) * (*xdim)]();

    pch = strtok (buffer, " ,\t\n");
    int i = 0;
    while (pch != NULL)
    {
        int a = atoi(pch);
        (*mask)[i++] = (a != 0);
        // printf (" <%d %d %s> ", i++, a, pch);
        pch = strtok (NULL, " ,\t\n");
    }
    return EXIT_SUCCESS;
}
*/
