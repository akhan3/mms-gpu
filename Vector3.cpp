#include "Vector3.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif

Vector3::Vector3() {
    x = 0;
    y = 0;
    z = 0;
}

Vector3::Vector3(const fptype x1, const fptype y1, const fptype z1) {
    // std::cout << "constructor: " << x1 << " " << y1 << " " << z1 << std::endl;
    init(x1, y1, z1);
}

// Vector3::Vector3(const fptype* elem) {
    // init(elem[0], elem[1], elem[2]);
// }

Vector3& Vector3::init(const fptype x1, const fptype y1, const fptype z1) {
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

Vector3 Vector3::operator*(const fptype k) const {
    return Vector3(k * x, k * y, k * z);
}

Vector3 Vector3::operator*(const Vector3 &b) const {
    return Vector3( x * b.x,
                    y * b.y,
                    z * b.z );
}

fptype Vector3::dot(const Vector3 &b) const {
    return x * b.x + y * b.y + z * b.z;
}

Vector3 Vector3::cross(const Vector3 &b) const {
    return Vector3( y * b.z - z * b.y,
                    z * b.x - x * b.z,
                    x * b.y - y * b.x   );
}

fptype Vector3::magnitude() const {
    fptype x2;
    fptype y2;
    fptype z2;
// #ifdef _OPENMP
    // #pragma omp parallel shared(x2, y2, z2)
    // {
        // if(omp_get_thread_num() == 0) x2 = x*x;
        // if(omp_get_thread_num() == 1) y2 = y*y;
        // if(omp_get_thread_num() == 0) z2 = z*z;
    // }
// #else
    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
// #endif
    return sqrt(x2 + y2 + z2);
}

fptype Vector3::colatitude() const {
    return acos(z / sqrt(x * x + y * y + z * z));
}

fptype Vector3::azimuth() const {
    return atan2(y, x);
}

void Vector3::print(FILE *fh) const {
    fprintf(fh, "%g, %g, %g\n", x, y, z);
}
