#ifndef _QUEUE_H_
#define _QUEUE_H_

#include "mydefs.hpp"

class Queue {
// data members
    void **contents;
public:
    int size;
    int count;
    int front;

// constructors and destructor
    Queue(int len);
    ~Queue();

// member functions
    void            enqueue(void* n);
    void*           dequeue();
    inline int      isEmpty() {
        return count == 0;
    }

    // void            printQueue();
};

#endif // #ifndef  _QUEUE_H_
