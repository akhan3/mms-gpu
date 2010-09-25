#include "Vector3.hpp"
// #ifdef _OPENMP
// #include <omp.h>
// #endif

HOSTDEVICE
Vector3::Vector3() {
    x = 0;
    y = 0;
    z = 0;
}

HOSTDEVICE
Vector3::Vector3(const fptype x1, const fptype y1, const fptype z1) {
    x = x1;
    y = y1;
    z = z1;
}

// Vector3::Vector3(const fptype* elem) {
    // init(elem[0], elem[1], elem[2]);
// }

// Vector3& Vector3::operator+=(const Vector3 &b) {
HOSTDEVICE
void Vector3::operator+=(const Vector3 &b) {
    x += b.x;
    y += b.y;
    z += b.z;
    // return *this;
}

HOSTDEVICE
Vector3 Vector3::operator+(const Vector3 &b) const {
    return Vector3( x + b.x,
                    y + b.y,
                    z + b.z );
}

HOSTDEVICE
Vector3 Vector3::operator-(const Vector3 &b) const {
    return Vector3( x - b.x,
                    y - b.y,
                    z - b.z );
}

HOSTDEVICE
Vector3 Vector3::operator/(const Vector3 &b) const {
    return Vector3( x / b.x,
                    y / b.y,
                    z / b.z );
}

HOSTDEVICE
Vector3 Vector3::operator*(const fptype k) const {
    return Vector3(k * x, k * y, k * z);
}

HOSTDEVICE
Vector3 Vector3::operator*(const Vector3 &b) const {
    return Vector3( x * b.x,
                    y * b.y,
                    z * b.z );
}

HOSTDEVICE
fptype Vector3::dot(const Vector3 &b) const {
    return x * b.x + y * b.y + z * b.z;
}

HOSTDEVICE
Vector3 Vector3::cross(const Vector3 &b) const {
    return Vector3( y * b.z - z * b.y,
                    z * b.x - x * b.z,
                    x * b.y - y * b.x   );
}

HOSTDEVICE
fptype Vector3::magnitude() const {
    return sqrt(x*x + y*y + z*z);
}

HOSTDEVICE
fptype Vector3::colatitude() const {
    return acos(z / sqrt(x * x + y * y + z * z));
}

HOSTDEVICE
fptype Vector3::azimuth() const {
    return atan2(y, x);
}

HOST
void Vector3::print(FILE *fh) const {
    fprintf(fh, "%g, %g, %g\n", x, y, z);
}
