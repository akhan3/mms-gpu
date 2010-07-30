#include <stdio.h>
// #include <iostream>
// #include <string.h>
#include <cmath>
#include "Box.hpp"

// constructor
Box::Box(unsigned int level1, unsigned int limit) {
    level = level1;
    parent = NULL;
    for(int i=0; i<4; i++)
        child[i] = NULL;
    for(int i=0; i<8; i++)
        neighbor[i] = NULL;
    for(int i=0; i<27; i++)
        interaction[i] = NULL;
    x = 0;
    y = 0;
    // cx = 0;
    // cy = 0;
    cx = (unsigned int)pow(2, limit-level) * x + ((unsigned int)pow(2, limit-level) - 1) * 0.5;
    cy = (unsigned int)pow(2, limit-level) * y + ((unsigned int)pow(2, limit-level) - 1) * 0.5;
    // potential = 0;
    pruned = 0;
    id = new byte[level+1]();
    id[0] = 0;
}

// destructor
Box::~Box() {
    delete []id;
    for(int i=0; i<4; i++)
        if (child[0] != NULL)
            delete child[i];
}

// Split up parent Box into 4 children Boxs
void Box::split(unsigned int limit) {
    for(int i=0; i<4; i++) {
        child[i] = new Box(level+1, limit);
        Box *n = child[i];
        n->parent = this;
        n->level = level+1;
    // calculate Box's id
        for(unsigned int j=0; j<=level; j++)
            n->id[j] = id[j];
        n->id[level+1] = i;
    // calculate spatial coordinates at level
        n->x = n->calc_x(n->level, n->id);
        n->y = n->calc_y(n->level, n->id);
    // calculate spatial coordinates at deepest level
        n->cx = (unsigned int)pow(2, limit-n->level) * n->x + ((unsigned int)pow(2, limit-n->level) - 1) * 0.5;
        n->cy = (unsigned int)pow(2, limit-n->level) * n->y + ((unsigned int)pow(2, limit-n->level) - 1) * 0.5;
    }
}


// Prune the tree from this node downward
void Box::prune() {
    pruned = 1;
    if (child[0] != NULL)
        for(int i=0; i<=3; i++) {
            child[i]->prune();
        }
}


// construct neighbors and intercation list
void Box::find_neighbors(Box* root) {
    unsigned int Nx[8];
    unsigned int Ny[8];
    byte *Nid = new byte[level+1]();
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
    if (level <= 1)
        return;
    Box *fullList[32];
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
} // find_neighbors function

// calculate x and y spatial coords from id
unsigned int Box::calc_x (unsigned int level1, byte *id1) {
    unsigned int x1 = 0;
    for(int j=level1; j>=1; j--)
        x1 |= (id1[j]&1) << (level1-j);
    return x1;
}
unsigned int Box::calc_y (unsigned int level1, byte *id1) {
    unsigned int y1 = 0;
    for(int j=level1; j>=1; j--)
        y1 |= ((id1[j]&2)>>1) << (level1-j);
    return y1;
}

// calculate id from x and y spatial coords
void Box::calc_id (byte *id1, unsigned int level1, unsigned int x1, unsigned int y1) {
    for(unsigned int k=0; k<=level1; k++)
        id1[k] = 0;
    for(int k=level1; k>=0; k--) {
        id1[k] |= (x1 & (1<<(level1-k))) >> (level1-k);
        id1[k] |= ((y1 & (1<<(level1-k))) >> (level1-k)) << 1;
    }
}

// Gather Box-id in a string
void Box::get_idstring(char *s) {
    int c = 0;
    c += sprintf(s, "[");
    for(unsigned int a=0; a<=level; a++)
        c += sprintf(s+c, "%d", id[a]);
    c += sprintf(s+c, "]");
}
