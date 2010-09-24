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
    // byte *id;
    unsigned int index;
    Box *parent;        // pointer to parent
    unsigned int x;
    unsigned int y;
    int pruned;

public:
// constructors and destructor
    Box(unsigned int level1, unsigned int index1, unsigned int limit);
    ~Box();

// member functions
    void            grow();
    void            prune();
    inline int      is_pruned() { return pruned; }
    void            create_tree_recurse(const unsigned int limit);
    void            find_neighbors_recurse(Box *root, const unsigned int limit);
private:
    void            split(const unsigned int limit);
    unsigned int    calc_x(const byte *id);
    unsigned int    calc_y(const byte *id);
    void            idfromindex(byte *id);
    void            calc_id(byte*id1, const unsigned int level1, const unsigned int x1, const unsigned int y1);
    void            get_idstring(char *s);
    void            find_neighbors(Box* root);
};

// #ifdef __CUDACC__
// #include "Box.cpp"
// #endif

#endif // #ifndef  _BOX_H_
