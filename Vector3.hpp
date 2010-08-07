#ifndef _VECTOR3_H_
#define _VECTOR3_H_

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <cmath>

struct Vector3 {
// data members
    float x,y,z;

// constructors
    Vector3();
    Vector3(const float x1, const float y1, const float z1);
    // Vector3(const float *elem);
    // Vector3(const float x1, const float y1, const float z1);

// member functions
// get functions

// set functions
    Vector3& init(const float x1, const float y1, const float z1);

// overloaded operators
    Vector3& operator+=(const Vector3 &b);

    Vector3 operator+(const Vector3 &b) const;
    Vector3 operator-(const Vector3 &b) const;
    Vector3 operator/(const Vector3 &b) const;
    Vector3 operator*(const float k) const;
    friend Vector3 operator*(const float k, const Vector3 &vec);
    Vector3 operator*(const Vector3 &b) const;  // element-by-element multiplication
    float dot(const Vector3 &b) const;        // dot product
    Vector3 cross(const Vector3 &b) const;      // cross product
    float magnitude() const;
    float colatitude() const;
    float azimuth() const;
    friend std::ostream& operator<<(std::ostream& output, const Vector3 &vec);
    void print(FILE *fh) const;
};

// friend Vector3 operator*(const float k, const Vector3 &vec);
inline Vector3 operator*(const float k, const Vector3 &vec) {
    return vec * k;
}

// friend ostream& operator<<(ostream& output, const Vector3 &vec);
inline std::ostream& operator<<(std::ostream& output, const Vector3 &vec) {
    output << vec.x << ", " << vec.y << ", " << vec.z;
    return output;
}

#endif // #ifndef  _VECTOR3_H_
