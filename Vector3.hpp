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
    Vector3();
    Vector3(const fptype x1, const fptype y1, const fptype z1);
    // Vector3(const fptype *elem);

// member functions
// overloaded operators
    // Vector3& operator+=(const Vector3 &b);
    void operator+=(const Vector3 &b);

    // inline Vector3 operator+(const Vector3 &b) const  { return Vector3(*this) += b; }
    Vector3 operator+(const Vector3 &b) const;
    Vector3 operator-(const Vector3 &b) const;
    Vector3 operator/(const Vector3 &b) const;
    Vector3 operator*(const fptype k) const;
    friend Vector3 operator*(const fptype k, const Vector3 &vec);
    Vector3 operator*(const Vector3 &b) const;  // element-by-element multiplication
    fptype dot(const Vector3 &b) const;        // dot product
    Vector3 cross(const Vector3 &b) const;      // cross product
    fptype magnitude() const;
    fptype colatitude() const;
    fptype azimuth() const;
    friend std::ostream& operator<<(std::ostream& output, const Vector3 &vec);
    void print(FILE *fh) const;
};

// friend Vector3 operator*(const fptype k, const Vector3 &vec);
inline Vector3 operator*(const fptype k, const Vector3 &vec) {
    return vec * k;
}

// friend ostream& operator<<(ostream& output, const Vector3 &vec);
inline std::ostream& operator<<(std::ostream& output, const Vector3 &vec) {
    output << vec.x << ", " << vec.y << ", " << vec.z;
    return output;
}

#endif // #ifndef  _VECTOR3_H_
