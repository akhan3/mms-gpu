#include <cmath>
#include "Cmpx.hpp"

// constructor
Cmpx::Cmpx() {
    x = 0;
    y = 0;
}

Cmpx::Cmpx(const float x, const float y) {
    init(x, y);
}

Cmpx& Cmpx::init(const float x1, const float y1) {
    x = x1;
    x = y1;
    return *this;
}

Cmpx Cmpx::conjugate() const {
    return Cmpx(x, -y);
}

Cmpx& Cmpx::operator+=(const Cmpx &z) {
    x += z.x;
    y += z.y;
    return *this;
}

Cmpx& Cmpx::operator-=(const Cmpx &z) {
    x -= z.x;
    y -= z.y;
    return *this;
}

Cmpx& Cmpx::operator*=(const Cmpx &z) {
    x = x*z.x - y*z.y;
    y = x*z.y + y*z.x;
    return *this;
}

// Cmpx& Cmpx::operator/=(const Cmpx &z) {
    // x = x*z.x - y*z.y;
    // y = x*z.y + y*z.x;
    // return *this;
// }

Cmpx& Cmpx::operator*=(const float k) {
    x *= k;
    y *= k;
    return *this;
}
