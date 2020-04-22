#include <stdio.h>
#include <stdlib.h>
#include "Queue.hpp"


// constructors
// Queue::Queue(int len) {
    // size = len;
    // contents = new void*[len]();
    // if(contents == NULL) {
        // fprintf(stderr, "%s:%d Error allocating memory\n", __FILE__, __LINE__);
        // return;
    // }
    // front = 0;
    // count = 0;
// }

HOST
Queue::Queue(int len, void **queue_mem) {
    size = len;
    contents = queue_mem;
    front = 0;
    count = 0;
}

// destructor
// Queue::~Queue() {
    // delete []contents;
// }

HOST
void Queue::enqueue(void* n) {
    if(count >= size) {
        // printf("FATAL ERROR\n");
        // printf("FATAL ERROR: Cannot add more items. Queue is already full!\n");
        // printf("FATAL ERROR\n");
        return;
    }
    contents[(count+front) % size] = n;
    count++;
}

HOST
void* Queue::dequeue() {
    if(isEmpty()) {
        // printf("FATAL ERROR\n");
        // printf("FATAL ERROR: Cannot pick any item. Queue is already empty!\n");
        // printf("FATAL ERROR\n");
        return NULL;
    }
    void *n1 = contents[front];
    count--;
    front++;
    front = front % size;
    return n1;
}

// void Queue::printQueue() {
    // printf("Q = [");
    // for(int i=0; i<size; i++)
        // if(contents[i]->ch <=32)
            // printf("%c", 219);
        // else
            // printf("%c", contents[i]->ch);
    // printf("] count = %d\n", count);
    // printf("     ");
    // for(int i=0; i<size; i++) {
        // if(i == front)
            // printf("^");
        // else
            // printf(" ");
    // }
    // printf("\n");
// }
