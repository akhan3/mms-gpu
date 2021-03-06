#ifndef _BOX_H_
#define _BOX_H_

#include "mydefs.hpp"

class Box {
// data members
public:
    Box *child[4];     // pointers to 4 children
    Box *neighbor[8];  // pointers to 8 neighbors
    Box *interaction[27];  // pointers to 27 interaction Boxs
    unsigned int level;
    fptype cx;
    fptype cy;
private:
    unsigned int index;
    Box *parent;        // pointer to parent
    unsigned int x;
    unsigned int y;
    int pruned;

public:
// constructors and destructor
    HOST Box(unsigned int level1, unsigned int index1, unsigned int limit);
    HOSTDEVICE Box();
    // ~Box();

// member functions
    HOST void            grow();
    HOST void            prune();
    HOST inline int      is_pruned() { return pruned; }
    HOST void            create_tree_bfs(const unsigned int limit, void **queue_mem);
    // void            create_tree_recurse(const unsigned int limit);
    HOST void            find_neighbors_recurse(Box *root, const unsigned int limit);
private:
    // void            split(Box *start_address, const unsigned int limit);
    HOST void            split(Box *firstchild_address, unsigned int limit);
    HOST unsigned int    calc_x(const byte *id);
    HOST unsigned int    calc_y(const byte *id);
    HOST void            idfromindex(byte *id);
    HOST void            calc_id(byte*id1, const unsigned int level1, const unsigned int x1, const unsigned int y1);
    HOST       void            get_idstring(char *s);
    HOST void            find_neighbors(Box* root);
};

#ifdef __CUDACC__
#include "Box.cu"
#endif

#endif // #ifndef  _BOX_H_
