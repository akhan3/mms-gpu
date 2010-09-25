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
    HOSTDEVICE Queue(int len, void **queue_mem);
    // ~Queue();

// member functions
    HOSTDEVICE void            enqueue(void* n);
    HOSTDEVICE void*           dequeue();
    HOSTDEVICE inline int      isEmpty() {
        return count == 0;
    }

    // void            printQueue();
};

#ifdef __CUDACC__
#include "Queue.cu"
#endif

#endif // #ifndef  _QUEUE_H_
