#ifndef _QUEUE_H_
#define _QUEUE_H_

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
    int             isEmpty();
    // void            printQueue();
};

#endif // #ifndef  _QUEUE_H_
