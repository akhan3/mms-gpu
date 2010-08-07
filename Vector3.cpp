#ifndef _VECTOR3_CPP_
#define _VECTOR3_CPP_

#include "Vector3.hpp"

Vector3::Vector3() {
    x = 0;
    y = 0;
    z = 0;
}

Vector3::Vector3(const float x1, const float y1, const float z1) {
    // std::cout << "constructor: " << x1 << " " << y1 << " " << z1 << std::endl;
    init(x1, y1, z1);
}

// Vector3::Vector3(const float* elem) {
    // init(elem[0], elem[1], elem[2]);
// }

Vector3& Vector3::init(const float x1, const float y1, const float z1) {
    // std::cout << "  init: " << x1 << " " << y1 << " " << z1 << std::endl;
    x = x1;
    y = y1;
    z = z1;
    // std::cout << "  init: " << x << " " << y << " " << z << std::endl;
    return *this;
}

Vector3& Vector3::operator+=(const Vector3 &b) {
    x += b.x;
    y += b.y;
    z += b.z;
    return *this;
}

Vector3 Vector3::operator+(const Vector3 &b) const {
    return Vector3( x + b.x,
                    y + b.y,
                    z + b.z );
}

Vector3 Vector3::operator-(const Vector3 &b) const {
    return Vector3( x - b.x,
                    y - b.y,
                    z - b.z );
}

Vector3 Vector3::operator/(const Vector3 &b) const {
    return Vector3( x / b.x,
                    y / b.y,
                    z / b.z );
}

Vector3 Vector3::operator*(const float k) const {
    return Vector3(k * x, k * y, k * z);
}

Vector3 Vector3::operator*(const Vector3 &b) const {
    return Vector3( x * b.x,
                    y * b.y,
                    z * b.z );
}

float Vector3::dot(const Vector3 &b) const {
    return x * b.x + y * b.y + z * b.z;
}

Vector3 Vector3::cross(const Vector3 &b) const {
    return Vector3( y * b.z - z * b.y,
                    z * b.x - x * b.z,
                    x * b.y - y * b.x   );
}

float Vector3::magnitude() const {
    return sqrt(x * x + y * y + z * z);
}

float Vector3::colatitude() const {
    return acos(z / sqrt(x * x + y * y + z * z));
}

float Vector3::azimuth() const {
    return atan2(y, x);
}

void Vector3::print(FILE *fh) const {
    fprintf(fh, "%g, %g, %g\n", x, y, z);
}

#endif // #ifndef  _VECTOR3_CPP_
