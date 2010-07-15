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
struct node {
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
        for(int i=0; i<=7; i++) 
            neighbor[i] = NULL;
        level = 0;
        x = 0;
        y = 0;
        id = (byte*) malloc((level+1) * sizeof(unsigned int));
        id[0] = 0;
    }

// Split up parent node into 4 children nodes
    void split() {
        for(int i=0; i<=3; i++) {
            child[i] = (node*) malloc(sizeof(node));
            node *n = child[i];
            
            n->parent = this;
            for(int j=0; j<=3; j++) 
                n->child[j] = NULL;
            n->level = level+1;
            n->id = (byte*) malloc((n->level+1) * sizeof(byte));
            for(int j=0; j<=level; j++)
                n->id[j] = id[j];
            n->id[level+1] = i;
            // extract coordinates at level
            n->x = 0;
            n->y = 0;
            for(int j=n->level; j>=1; j--) {
                n->x |= (n->id[j]&1)      << (n->level-j);
                n->y |= ((n->id[j]&2)>>1) << (n->level-j);
            }            
            for(int j=0; j<=7; j++) 
                n->neighbor[j] = NULL;
        }
    }


// construct neighbors list
    void find_neighbors(node* n, node* root) {
        if (n->parent == NULL)
            return;
        unsigned int Nx[8];
        unsigned int Ny[8];
        byte *Nid = (byte*) malloc((n->level+1) * sizeof(byte));
        Nx[0] = n->x-1; Ny[0] = n->y-1;
        Nx[1] = n->x  ; Ny[1] = n->y-1;
        Nx[2] = n->x+1; Ny[2] = n->y-1;
        Nx[3] = n->x-1; Ny[3] = n->y  ;
        Nx[4] = n->x+1; Ny[4] = n->y  ;
        Nx[5] = n->x-1; Ny[5] = n->y+1;
        Nx[6] = n->x  ; Ny[6] = n->y+1;
        Nx[7] = n->x+1; Ny[7] = n->y+1;
        char *idstring = (char*) malloc(100);
        n->get_idstring(idstring);
        // printf("neighbors of L%d%s=(%d,%d) are ", n->level, idstring, n->x, n->y);
        for(int j=0; j<=7; j++) {
            // printf("(%d,%d)=", Nx[j], Ny[j]);
            for(int k=0; k<=n->level; k++)
                Nid[k] = 0;
            // printf("[");
            for(int k=n->level; k>=0; k--) {
                Nid[k] |= (Nx[j] & (1<<(n->level-k))) >> (n->level-k);
                Nid[k] |= ((Ny[j] & (1<<(n->level-k))) >> (n->level-k)) << 1;
                // printf("{%d:%d}", k, Nid[k]);
            }
            for(int k=0; k<=n->level; k++)
                // printf("%d", Nid[k]);
            // printf("]");
            // join the neighbors
            if (Nid[0] == 0) {
                node *nt = root;
                for(int k=1; k<=n->level; k++) {
                    nt = nt->child[Nid[k]];
                }
                n->neighbor[j] = nt;
                n->neighbor[j]->get_idstring(idstring);
                // printf("{%d:(%d,%d)=%s}", j, n->neighbor[j]->x, n->neighbor[j]->y, idstring);
            }
            // printf(" ");
        }
        // NEWLINE;

        n->get_idstring(idstring);
        // printf("Neighbors of L%d%s=(%d,%d) are ", n->level, idstring, n->x, n->y);
        for(int j=0; j<=7; j++) {
            // printf("%p ", n->neighbor[j]);
            if (n->neighbor[j] != NULL) {
                // printf("(%d,%d)", n->neighbor[j]->x, n->neighbor[j]->y);
                n->neighbor[j]->get_idstring(idstring);
                // printf("=");
                // printf("%s ", idstring);
            }
            else {
                // printf("( , )=");
                // printf("[] ");
            }
        }
        // NEWLINE;
    }

// Gather node-id in a string
    void get_idstring(char *s) {
        int c = 0;
        c += sprintf(s, "[");
        for(int a=0; a<=level; a++)
            c += sprintf(s+c, "%d", id[a]);
        c += sprintf(s+c, "]");
    }
};  // struct definition finished



//******************************************************************

void create_tree_recurs(node *thisnode, const int limit) {
    if (thisnode->level >= limit)
        return;
    // printf("tree recursive call...\n");
    thisnode->split();
    for(int i=0; i<=3; i++)
        create_tree_recurs(thisnode->child[i], limit);
}

void find_neighbors_recurs(node *thisnode, node *root, const int limit) {
    // printf("find recursive call...\n");
    thisnode->find_neighbors(thisnode, root);
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
    const int N = 1024;
    const int logN = ceil(log2(N) / log2(4));
    printf("N = %d, log4(N) = %d\n", N, logN);
    printf("sizeof(node) = %ld\n", sizeof(node));
    
    node *root = (node*) malloc(sizeof(node));
    *root = node();

// generate the tree
    create_tree_recurs(root, logN);
    find_neighbors_recurs(root, root, logN);
    
    



// traverse to a random deepest node
    node* n = root;
    char *idstring = (char*) malloc(100);
    int c[3] = {0,2,3};
    int cc = 0;
    while (n->child[0] != NULL)
    {
        n = n->child[rand_atob(0,4)];
        // n = n->child[c[cc++]];
        n->get_idstring(idstring);
        printf("this node is at L%d: %s ", n->level, idstring);
        printf("(%d,%d)", n->x, n->y);
        NEWLINE;
    }
    printf("and its neighbors are");
    for(int i=0; i<8; i++)
        if (n->neighbor[i] != NULL)
            printf("(%d,%d) ", n->neighbor[i]->x, n->neighbor[i]->y);
    NEWLINE;
    
    
    unsigned int size = (unsigned int)pow(2, n->level);
    byte *matrix = (byte*) malloc(size*size);
    printf("size of matrix = %dx%d\n", size, size);
    for(int i=0; i<size; i++) {
        for(int j=0; j<size; j++) {
            matrix[i*size+j] = '#';
        }
    }
    matrix[n->y*size+n->x] = '#';
    for(int i=0; i<8; i++) {
        if (n->neighbor[i] != NULL)
            matrix[n->neighbor[i]->y*size+n->neighbor[i]->x] = 219;
    }
    // for(int i=0; i<size; i++) {
    for(int i=size-1; i>=0; i--) {        
        for(int j=0; j<size; j++) {
            printf("%c ", matrix[i*size+j]);
        }
        NEWLINE;
    }
    
    
    
    // while (n->parent != NULL)
    // {
        // n = n->parent;
        // n->get_idstring(idstring);
        // printf("this node is at L%d: %s ", n->level, idstring);
        // printf("(%d,%d)", n->x, n->y);
        // NEWLINE;
    // }
    

    free(root);
	return EXIT_SUCCESS;
}
