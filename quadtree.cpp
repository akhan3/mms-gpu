#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>

using namespace std;
using std::cout;
using std::endl;

typedef unsigned char byte;
#define NEWLINE printf("\n");


//*******************************************************************
//*********** Node class definition *********************************
//*******************************************************************
class node {
    public:
	node *parent;	    // pointer to parent
	node *child[4];     // pointers to 4 children
	node *neighbor[8];  // pointers to 8 neighbors
	node *interaction[27];  // pointers to 27 interaction nodes
	unsigned int level;
	unsigned int x;
	unsigned int y;
	byte *id;
    
// methods

// constructor - used only for root node
    node() {
        parent = NULL;
        for(int i=0; i<8; i++) 
            neighbor[i] = NULL;
        for(int i=0; i<27; i++) 
            interaction[i] = NULL;
        level = 0;
        x = 0;
        y = 0;
        id = (byte*) malloc((level+1) * sizeof(unsigned int));
        id[0] = 0;
    }

// calculate x and y spatial coords from id
    inline unsigned int calc_x(unsigned int level1, byte *id1) {
        unsigned int x1 = 0;
        for(int j=level1; j>=1; j--)
            x1 |= (id1[j]&1) << (level1-j);
        return x1;
    }
    inline unsigned int calc_y(unsigned int level1, byte *id1) {
        unsigned int y1 = 0;
        for(int j=level1; j>=1; j--)
            y1 |= ((id1[j]&2)>>1) << (level1-j);
        return y1;
    }
    // unsigned int calc_x() {
        // x = 0;
        // for(int j=level; j>=1; j--)
            // x |= (id[j]&1) << (level-j);
        // return x;
    // }

// calculate id from x and y spatial coords
    inline void calc_id(byte*id1, unsigned int level1, unsigned int x1, unsigned int y1) {
        for(unsigned int k=0; k<=level1; k++)
            id1[k] = 0;
        for(int k=level1; k>=0; k--) {
            id1[k] |= (x1 & (1<<(level1-k))) >> (level1-k);
            id1[k] |= ((y1 & (1<<(level1-k))) >> (level1-k)) << 1;
        }
    }

// Split up parent node into 4 children nodes
    void split() {
        for(int i=0; i<4; i++) {
            child[i] = (node*) malloc(sizeof(node));
            node *n = child[i];
            n->parent = this;
            for(int j=0; j<4; j++) 
                n->child[j] = NULL;
            n->level = level+1;
        // calculate node's id
            n->id = (byte*) malloc((n->level+1) * sizeof(byte));
            for(unsigned int j=0; j<=level; j++)
                n->id[j] = id[j];
            n->id[level+1] = i;
        // calculate spatial coordinates at level
            n->x = n->calc_x(n->level, n->id);
            n->y = n->calc_y(n->level, n->id);
        // initialize neighbor and interaction list to NULL
            for(int j=0; j<8; j++) 
                n->neighbor[j] = NULL;
            for(int j=0; j<27; j++) 
                n->interaction[j] = NULL;
        }
    }


// construct neighbors and intercation list
    void find_neighbors(node* root) {
        if (parent == NULL)
            return;
        // char *idstring = (char*) malloc(100);
        unsigned int Nx[8];
        unsigned int Ny[8];
        byte *Nid = (byte*) malloc(level+1);
        Nx[0] = x-1; Ny[0] = y-1;
        Nx[1] = x  ; Ny[1] = y-1;
        Nx[2] = x+1; Ny[2] = y-1;
        Nx[3] = x-1; Ny[3] = y  ;
        Nx[4] = x+1; Ny[4] = y  ;
        Nx[5] = x-1; Ny[5] = y+1;
        Nx[6] = x  ; Ny[6] = y+1;
        Nx[7] = x+1; Ny[7] = y+1;
        for(int j=0; j<=7; j++) {   // for each of 8 neighbors
            calc_id(Nid, level, Nx[j], Ny[j]);
            // put the pointer to neighbor in neighbors list
            for(unsigned int k=0; k<=level; k++) {
                if (Nid[0] == 0) {
                    node *nt = root;
                    for(unsigned int k=1; k<=level; k++)
                        nt = nt->child[Nid[k]];
                    neighbor[j] = nt;
                }
            }
        }
        
        // Now construct the interaction list
        // by taking difference of two sets
        if (level <= 1)
            return;
        node *fullList[32];
        int c = 0;
        for(int j=0; j<8; j++) {                // parent's neighbor index
            if (parent->neighbor[j]) {          // if parent has this neighbor?
                for(int k=0; k<4; k++) {        // parent's neighbor's child index
                    fullList[c++] = parent->neighbor[j]->child[k];
                }
            }
        }
        int c1 = 0;
        for (int i=0; i<c; i++) {
            int foundInNeighbor = 0; 
            for (int j=0; j<8; j++) { 
                if (fullList[i] == neighbor[j]) {
                    foundInNeighbor = 1; 
                    break; 
                }
            }
            if (!foundInNeighbor)
                interaction[c1++] = fullList[i];
        }
    } // function


// Gather node-id in a string
    void get_idstring(char *s) {
        int c = 0;
        c += sprintf(s, "[");
        for(unsigned int a=0; a<=level; a++)
            c += sprintf(s+c, "%d", id[a]);
        c += sprintf(s+c, "]");
    }
};  // struct definition finished



//******************************************************************

void create_tree_recurs(node *thisnode, const unsigned int limit) {
    if (thisnode->level >= limit)
        return;
    // printf("tree recursive call...\n");required by July 21. (Min. 4; max. 12)
    thisnode->split();
    for(int i=0; i<=3; i++)
        create_tree_recurs(thisnode->child[i], limit);
}

void find_neighbors_recurs(node *thisnode, node *root, const unsigned int limit) {
    // printf("find recursive call...\n");
    thisnode->find_neighbors(root);
    if (thisnode->level >= limit)
        return;
    for(int i=0; i<=3; i++)
        find_neighbors_recurs(thisnode->child[i], root, limit);
}

int rand_atob(const int a, const int b) {
    double r = rand() / (double)RAND_MAX;
    r = a + (b-a) * r;
    return (int)r;
}


//*************************************************************************//
//******************** Main function **************************************//
//*************************************************************************//
int main(void) {
    srand(time(NULL));
    const int N = 256;
    // const int N = 1000*1000;
    const int logN = ceil(log2(N) / log2(4));
    printf("N = %d, log4(N) = %d\n", N, logN);
    printf("sizeof(node) = %ld\n", sizeof(node));
    
    node *root = (node*) malloc(sizeof(node));
    *root = node();

// generate the tree
    create_tree_recurs(root, logN);
    find_neighbors_recurs(root, root, logN);
    
    
// traverse to a random node at the deepest level
    // int c[5] = {0,0,0,0,0};
    // int cc = 0;
    node* n = root;
    char *idstring = (char*) malloc(100);
    while (n->child[0] != NULL)
    {
        n = n->child[rand_atob(0,4)];
        // n = n->child[c[cc++]];
        n->get_idstring(idstring);
        printf("this node is at L%d%s(%d,%d)", n->level, idstring, n->x, n->y);
        NEWLINE;
        // if(rand() / (double)RAND_MAX < 0.2) break;
    }
    printf("its neighbors are ");
    for(int i=0; i<8; i++)
        if (n->neighbor[i] != NULL) {
            n->neighbor[i]->get_idstring(idstring);
            printf("%s(%d,%d) ", idstring, n->neighbor[i]->x, n->neighbor[i]->y);
            // printf("(%d,%d) ", n->neighbor[i]->x, n->neighbor[i]->y);
        }
    NEWLINE;
    printf("its interactions are ");
    for(int i=0; i<27; i++)
        if (n->interaction[i] != NULL) {
            n->interaction[i]->get_idstring(idstring);
            printf("%s(%d,%d) ", idstring, n->interaction[i]->x, n->interaction[i]->y);
        }
    NEWLINE;
    
    
    unsigned int size = (unsigned int)pow(2, n->level);
    byte *matrix = (byte*) malloc(size*size);
    // printf("size of matrix = %dx%d\n", size, size);
    for(unsigned int i=0; i<size; i++) {
        for(unsigned int j=0; j<size; j++) {
            matrix[i*size+j] = 'o';
        }
    }
    matrix[n->y*size+n->x] = 'S';
    for(int i=0; i<8; i++) {
        if (n->neighbor[i] != NULL)
            matrix[n->neighbor[i]->y*size+n->neighbor[i]->x] = 'N';
    }
    for (int i=0; i<27; i++) {
        if (n->interaction[i] != NULL)
            matrix[n->interaction[i]->y*size+n->interaction[i]->x] = '#';
    }
    for(int i=size-1; i>=0; i--) {
        for(unsigned int j=0; j<size; j++) {
            printf("%c ", matrix[i*size+j]);
        }
        NEWLINE;
    }
    
    // int a = n->child[1] == n->child[0];
    // int a[10];
    // printf("a=%ld\n", sizeof(n->child[1] == n->child[0]));
    // printf("%p %p %ld %d\n", n->parent, n->parent->parent, n->parent-n->parent->parent, n->parent>n->parent->parent);


    free(root);
	return EXIT_SUCCESS;
}
