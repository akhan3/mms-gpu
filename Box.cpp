#include <stdio.h>
#include <stdlib.h>
// #include <iostream>
// #include <string.h>
#include <cmath>
#include "Box.hpp"

// constructor
Box::Box(unsigned int level1, unsigned int index1, unsigned int limit) {
    level = level1;
    index = index1;
    pruned = 0;
    parent = NULL;
    for(int i=0; i<4; i++)
        child[i] = NULL;
    for(int i=0; i<8; i++)
        neighbor[i] = NULL;
    for(int i=0; i<27; i++)
        interaction[i] = NULL;

// calculate id from index
    // byte *id = new byte[level+1]();
    byte *id = new byte[limit+1]();
    if(id == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return;
    }
    idfromindex(id);

// calculate spatial coordinates at level
    x = calc_x(id);
    y = calc_y(id);
    delete []id;

// calculate spatial coordinates at deepest level
    cx = (unsigned int)pow(2, limit - level) * x + ((unsigned int)pow(2, limit - level) - 1) * 0.5;
    cy = (unsigned int)pow(2, limit - level) * y + ((unsigned int)pow(2, limit - level) - 1) * 0.5;
}

// destructor
Box::~Box() {
    for(int i=0; i<4; i++)
        if(child[0] != NULL)
            delete child[i];
}

// Split up parent Box into 4 children Boxs
void Box::split(unsigned int limit) {
    for(int i=0; i<4; i++) {
        child[i] = new Box(level+1, index*4 + i, limit);
        if(child[i] == NULL) {
            fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
            return;
        }
        child[i]->parent = this;
    }
}

// calculate id from index
void Box::idfromindex(byte *id) {
    int indext = index;
    id[level] = indext % 4; // indext % 4;
    int l = level - 1;
    while((indext /= 4) > 0) {
        id[l] = indext % 4;
        l--;
    }
    // for(; l >= 0; l--)
        // id[l] = 0;
}

// Prune the tree from this node downward
void Box::prune() {
    pruned = 1;
    if(child[0] != NULL)
        for(int i=0; i<=3; i++)
            child[i]->prune();
}

// Grow the tree from this node downward
void Box::grow() {
    pruned = 0;
    if(child[0] != NULL)
        for(int i=0; i<=3; i++)
            child[i]->grow();
}


// construct neighbors and intercation list
void Box::find_neighbors(Box* root) {
    unsigned int Nx[8];
    unsigned int Ny[8];
    byte *Nid = new byte[level+1]();
    if(Nid == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return;
    }
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
            if(Nid[0] == 0) {
                Box *nt = root;
                for(unsigned int k=1; k<=level; k++)
                    nt = nt->child[Nid[k]];
                neighbor[j] = nt;
            }
        }
    }
    delete []Nid;

    // Now construct the interaction list
    // by taking difference of two sets
    if(level <= 1)
        return;
    Box *fullList[32];
    int c = 0;
    for(int j=0; j<8; j++) {                // parent's neighbor index
        if(parent->neighbor[j]) {          // if parent has this neighbor?
            for(int k=0; k<4; k++) {        // parent's neighbor's child index
                fullList[c++] = parent->neighbor[j]->child[k];
            }
        }
    }
    int c1 = 0;
    for (int i=0; i<c; i++) {
        int foundInNeighbor = 0;
        for (int j=0; j<8; j++) {
            if(fullList[i] == neighbor[j]) {
                foundInNeighbor = 1;
                break;
            }
        }
        if(!foundInNeighbor)
            interaction[c1++] = fullList[i];
    }
} // find_neighbors function

// calculate x and y spatial coords from id
unsigned int Box::calc_x(const byte *id) {
    unsigned int x1 = 0;
    for(int j=level; j>=1; j--)
        x1 |= (id[j]&1) << (level-j);
    return x1;
}
unsigned int Box::calc_y(const byte *id) {
    unsigned int y1 = 0;
    for(int j=level; j>=1; j--)
        y1 |= ((id[j]&2)>>1) << (level-j);
    return y1;
}

// calculate id from x and y spatial coords
void Box::calc_id(byte *id, unsigned int level1, unsigned int x1, unsigned int y1) {
    for(unsigned int k=0; k<=level1; k++)
        id[k] = 0;
    for(int k=level1; k>=0; k--) {
        id[k] |= (x1 & (1<<(level1-k))) >> (level1-k);
        id[k] |= ((y1 & (1<<(level1-k))) >> (level1-k)) << 1;
    }
}

// Gather Box-id in a string
void Box::get_idstring(char *s) {
    byte *id = new byte[level+1]();
    if(id == NULL) {
        fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        return;
    }
    idfromindex(id);
    int c = 0;
    c += sprintf(s, "[");
    for(unsigned int a=0; a<=level; a++)
        c += sprintf(s+c, "%d", id[a]);
    c += sprintf(s+c, "]");
    delete []id;
}

void Box::create_tree_recurse(const unsigned int limit) {
    // Recursion
    if(level < limit) {
        // function to perform on node
        split(limit);
        for(int i=0; i<=3; i++)
            child[i]->create_tree_recurse(limit);
    }
}

void Box::find_neighbors_recurse(Box *root, const unsigned int limit) {
    // function to perform on node
    find_neighbors(root);
    // Recursion
    if(level < limit) {
        for(int i=0; i<=3; i++)
            child[i]->find_neighbors_recurse(root, limit);
    }
}

