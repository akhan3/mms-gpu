#ifndef _VECTOR3_H_
#define _VECTOR3_H_

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cmath>
#include "mydefs.hpp"

struct Vector3 {
// data members
    fptype x,y,z;

// constructors
    HOSTDEVICE Vector3();
    HOSTDEVICE Vector3(const fptype x1, const fptype y1, const fptype z1);
    // Vector3(const fptype *elem);

// member functions
// overloaded operators
    // Vector3& operator+=(const Vector3 &b);
    HOSTDEVICE void operator+=(const Vector3 &b);

    // inline Vector3 operator+(const Vector3 &b) const  { return Vector3(*this) += b; }
    HOSTDEVICE Vector3 operator+(const Vector3 &b) const;
    HOSTDEVICE Vector3 operator-(const Vector3 &b) const;
    HOSTDEVICE Vector3 operator/(const Vector3 &b) const;
    HOSTDEVICE Vector3 operator*(const fptype k) const;
    HOSTDEVICE friend Vector3 operator*(const fptype k, const Vector3 &vec);
    HOSTDEVICE Vector3 operator*(const Vector3 &b) const;  // element-by-element multiplication
    HOSTDEVICE fptype dot(const Vector3 &b) const;        // dot product
    HOSTDEVICE Vector3 cross(const Vector3 &b) const;      // cross product
    HOSTDEVICE fptype magnitude() const;
    HOSTDEVICE fptype colatitude() const;
    HOSTDEVICE fptype azimuth() const;
    HOST friend std::ostream& operator<<(std::ostream& output, const Vector3 &vec);
    HOST void print(FILE *fh) const;
};

// friend Vector3 operator*(const fptype k, const Vector3 &vec);
HOSTDEVICE inline Vector3 operator*(const fptype k, const Vector3 &vec) {
    return vec * k;
}

// friend ostream& operator<<(ostream& output, const Vector3 &vec);
HOST inline std::ostream& operator<<(std::ostream& output, const Vector3 &vec) {
    output << vec.x << ", " << vec.y << ", " << vec.z;
    return output;
}

#ifdef __CUDACC__
#include "Vector3.cu"
#endif

#endif // #ifndef  _VECTOR3_H_
