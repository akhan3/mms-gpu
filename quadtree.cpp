#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <iostream>
#include <cmath>
#include "Box.hpp"
#define NEWLINE printf("\n");


using namespace std;
using std::cout;
using std::endl;


// prototype for FMM function
int fmm_calc(Box *root, unsigned int logN);


void create_tree_recurse(Box *thisBox, const unsigned int limit) {
    // Recursion
    if (thisBox->level < limit) {
        // function to perform on node
        thisBox->split(limit);
        for(int i=0; i<=3; i++)
            create_tree_recurse(thisBox->child[i], limit);
    }
}

void find_neighbors_recurse(Box *thisBox, Box *root, const unsigned int limit) {
    // function to perform on node
    thisBox->find_neighbors(root);
    // Recursion
    if (thisBox->level < limit) {
        for(int i=0; i<=3; i++)
            find_neighbors_recurse(thisBox->child[i], root, limit);
    }
}

int rand_atob(const int a, const int b) {
    double r = rand() / (double)RAND_MAX;
    r = a + (b-a) * r;
    return (int)r;
}

float frand_atob(const float a, const float b) {
    double r = rand() / (double)RAND_MAX;
    r = a + (b-a) * r;
    return (float)r;
}

//*************************************************************************//
//******************** Main function **************************************//
//*************************************************************************//
int main(int argc, char **argv) 
{
    unsigned int seed = time(NULL);
    unsigned int N = 256;
    if(argc >= 2) {
        sscanf(argv[1], "%u", &N);
        assert(N > 0);
    }
    if(argc >= 3) {
        sscanf(argv[2], "%u", &seed);
        assert(seed > 0);
    }
    srand(seed);
    const unsigned int logN = ceil(log2(N) / log2(4));
    N = (unsigned int)pow(4, logN);
    printf("N = %d, log4(N) = %d\n", N, logN);
    printf("sizeof(Box) = %ld\n", sizeof(Box));
    


    Box *root = new Box(0, logN);

// generate the tree
    create_tree_recurse(root, logN);
    find_neighbors_recurse(root, root, logN);
    
    
// // traverse to a random Box at the deepest level
// // ============================================================================
    // // int c[5] = {0,0,0,0,0};
    // // int cc = 0;
    // Box* n = root;
    // char idstring[100];
    // n->get_idstring(idstring);
    // printf("this Box is at L%d%s(%d,%d) = L%d(%.1f,%.1f)", n->level, idstring, n->x, n->y, logN, n->cx, n->cy);
    // NEWLINE;
    // while (n->child[0] != NULL)
    // {
        // n = n->child[rand_atob(0,4)];
        // // n = n->child[c[cc++]];
        // n->get_idstring(idstring);
        // printf("this Box is at L%d%s(%d,%d) = L%d(%.1f,%.1f)", n->level, idstring, n->x, n->y, logN, n->cx, n->cy);
        // NEWLINE;
        // // if(rand() / (double)RAND_MAX < 0.2) break;
    // }
    // // printf("its neighbors are ");
    // for(int i=0; i<8; i++)
        // if (n->neighbor[i] != NULL) {
            // n->neighbor[i]->get_idstring(idstring);
            // // printf("%s(%d,%d) ", idstring, n->neighbor[i]->x, n->neighbor[i]->y);
            // // printf("(%d,%d) ", n->neighbor[i]->x, n->neighbor[i]->y);
        // }
    // NEWLINE;
    // // printf("its interactions are ");
    // for(int i=0; i<27; i++)
        // if (n->interaction[i] != NULL) {
            // n->interaction[i]->get_idstring(idstring);
            // // printf("%s(%d,%d) ", idstring, n->interaction[i]->x, n->interaction[i]->y);
        // }
    // NEWLINE;
    

// // print out a matrix showing sample neighbor and interaction list
// // ============================================================================
    // unsigned int size = (unsigned int)pow(2, n->level);
    // byte *matrix = (byte*) malloc(size*size);
    // // printf("size of matrix = %dx%d\n", size, size);
    // for(unsigned int i=0; i<size; i++) {
        // for(unsigned int j=0; j<size; j++) {
            // matrix[i*size+j] = 'o';
        // }
    // }
    // matrix[n->y*size+n->x] = 'S';
    // for(int i=0; i<8; i++) {
        // if (n->neighbor[i] != NULL)
            // matrix[n->neighbor[i]->y*size+n->neighbor[i]->x] = 'N';
    // }
    // for (int i=0; i<27; i++) {
        // if (n->interaction[i] != NULL)
            // matrix[n->interaction[i]->y*size+n->interaction[i]->x] = '#';
    // }
    // // for(int i=size-1; i>=0; i--) {  // axis xy
    // for(int i=0; i<size; i++) {     // axis ij
        // for(unsigned int j=0; j<size; j++) {
            // printf("%c ", matrix[i*size+j]);
        // }
        // NEWLINE;
    // }
    // free(matrix);
    


// call the function for main FMM calculation 
    fmm_calc(root, logN);

    


    delete root;
    printf("SEED = %d\n", seed);
	return EXIT_SUCCESS;
}
