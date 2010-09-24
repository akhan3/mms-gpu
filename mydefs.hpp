#ifndef _MYDEFS_H_
#define _MYDEFS_H_

#define NEWLINE printf("\n");

typedef float fptype;
typedef unsigned char byte;

// Cuda directives
#ifdef __CUDACC__
#define HOST __host__
#define DEVICE __device__
#define HOSTDEVICE __host__ __device__
#define ALIGN16 __align__(16)
#else
#define HOST
#define DEVICE
#define HOSTDEVICE
#define ALIGN16
#endif

#endif // #ifndef  _MYDEFS_H_
