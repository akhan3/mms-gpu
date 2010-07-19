#ifndef _QUEUE_H_
#define _QUEUE_H_

#include "Box.hpp"

class Queue {
// data members
    Box **contents;
public:
    int size;
    int count;
    int front;

// constructors and destructor
    Queue(int len);
    ~Queue();

// member functions
    // void            enqueue(char n);
    // char            dequeue();
    void            enqueue(Box* n);
    Box*            dequeue();
   int             isEmpty();
    void            printQueue();
};

#endif // #ifndef  _QUEUE_H_
