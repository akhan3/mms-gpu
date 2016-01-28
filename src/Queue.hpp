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
    // Queue(int len);
    HOST Queue(int len, void **queue_mem);
    // ~Queue();

// member functions
    HOST void            enqueue(void* n);
    HOST void*           dequeue();
    HOST inline int      isEmpty() {
        return count == 0;
    }

    // void            printQueue();
};

#ifdef __CUDACC__
#include "Queue.cu"
#endif

#endif // #ifndef  _QUEUE_H_
