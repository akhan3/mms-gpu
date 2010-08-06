#ifndef _BOX_H_
#define _BOX_H_

typedef unsigned char byte;

class Box {
// data members
public:
    byte *id;
    Box *parent;        // pointer to parent
    Box *child[4];     // pointers to 4 children
    Box *neighbor[8];  // pointers to 8 neighbors
    Box *interaction[27];  // pointers to 27 interaction Boxs
    unsigned int level;
    unsigned int x;
    unsigned int y;
    float cx;
    float cy;
private:
    int pruned;
    // float potential;

public:
// constructors and destructor
    Box(unsigned int level1, unsigned int limit);
    ~Box();

// member functions
    void            split(unsigned int limit);
    void            prune();
    void            grow();
    inline int      is_pruned() { return pruned; }
    void            find_neighbors(Box* root);
    unsigned int    calc_x(unsigned int level1, byte *id1);
    unsigned int    calc_y(unsigned int level1, byte *id1);
    void            calc_id(byte*id1, unsigned int level1, unsigned int x1, unsigned int y1);
    void            get_idstring(char *s);
};

#endif // #ifndef  _BOX_H_
