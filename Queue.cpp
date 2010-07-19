#include "Queue.hpp"
#include <stdio.h>
#include <stdlib.h>


// constructor
Queue::Queue(int len) {
    size = len;
    contents = (Box**)malloc(len * sizeof(Box*));
    // contents = new char[size]();
    front = 0;
    count = 0;
}

// destructor
Queue::~Queue() {
    free(contents);
    // delete []contents;
}

// void Queue::enqueue(char n) {
    // if(count >= size) {
        // printf("FATAL ERROR\n");
        // printf("FATAL ERROR: Queue is already full!\n");
        // printf("FATAL ERROR\n");
        // return;
    // }
    // Box *b = new Box(n);
    // contents[(count+front) % size] = b;
    // count++;
// }

// char Queue::dequeue() {
    // if (isEmpty()) {
        // printf("FATAL ERROR\n");
        // printf("FATAL ERROR: Queue is already empty!\n");
        // printf("FATAL ERROR\n");
        // return NULL;
    // }
    // char n1 = contents[front]->ch;
    // count--;
    // front++; front = front % size;
    // return n1;    
// }

void Queue::enqueue(Box* n) {
    if(count >= size) {
        printf("FATAL ERROR\n");
        printf("FATAL ERROR: Queue is already full!\n");
        printf("FATAL ERROR\n");
        return;
    }
    contents[(count+front) % size] = n;
    count++;
}

Box* Queue::dequeue() {
    if (isEmpty()) {
        printf("FATAL ERROR\n");
        printf("FATAL ERROR: Queue is already empty!\n");
        printf("FATAL ERROR\n");
        return NULL;
    }
    Box *n1 = contents[front];
    count--;
    front++; front = front % size;
    return n1;    
}

int Queue::isEmpty() {
    return count == 0;
}

void Queue::printQueue() {
    printf("Q = [");
    for(int i=0; i<size; i++)
        if(contents[i]->ch <=32)
            printf("%c", 219);
        else
            printf("%c", contents[i]->ch);
    printf("] count = %d\n", count);
    printf("     ");
    for(int i=0; i<size; i++) {
        if (i == front) 
            printf("^");
        else
            printf(" ");
    }
    printf("\n");    
}
