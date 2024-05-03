#ifndef MATH_HPP
#define MATH_HPP

#include "real.hpp"
#include <cassert>
#include <cmath>

class Vector3f {
public:
  real x;
  real y;
  real z;

  Vector3f() : x(0), y(0), z(0) {}
  Vector3f(real x, real y, real z) : x(x), y(y), z(z) {}
  Vector3f(const Vector3f &v) : x(v.x), y(v.y), z(v.z) {}

  Vector3f &operator=(const Vector3f &v) {
    x = v.x;
    y = v.y;
    z = v.z;
    return *this;
  };

  Vector3f &normalize() {
    real length_sqr = x*x + y*y + z*z;
    real length = std::sqrt(length_sqr);
    x /= length; y /= length; z /= length;
    return *this;
  }
  Vector3f operator-(const Vector3f &v) const {
    return Vector3f(x - v.x, y - v.y, z - v.z);
  }
  Vector3f operator+(const Vector3f &v) const {
    return Vector3f(x + v.x, y + v.y, z + v.z);
  }
  Vector3f operator/(real d) const {
    return Vector3f(x / d, y / d, z / d);
  }
};

inline real distance(const Vector3f &v1, const Vector3f &v2) {
  Vector3f &&d = v1 - v2;
  return std::sqrt(d.x*d.x + d.y*d.y + d.z*d.z);
}

inline Vector3f cross(const Vector3f &v1, const Vector3f &v2) {
  return Vector3f(v1.y*v2.z - v1.z*v2.y,
                  v1.z*v2.x - v1.x*v2.z,
                  v1.x*v2.y - v1.y*v2.x);
}

inline real dot(const Vector3f &v1, const Vector3f &v2) {
  return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

class Quadric4f {
public:
  real
  q11, q12, q13, q14,
       q22, q23, q24,
            q33, q34,
                 q44;

  Quadric4f() :
    q11(0), q12(0), q13(0), q14(0),
    q22(0), q23(0), q24(0),
    q33(0), q34(0),
    q44(0) {}
  Quadric4f(real q11, real q12, real q13, real q14,
            real q22, real q23, real q24,
            real q33, real q34,
            real q44) :
    q11(q11), q12(q12), q13(q13), q14(q14),
    q22(q22), q23(q23), q24(q24),
    q33(q33), q34(q34),
    q44(q44) {}
  Quadric4f(const Quadric4f &q) :
    q11(q.q11), q12(q.q12), q13(q.q13), q14(q.q14),
    q22(q.q22), q23(q.q23), q24(q.q24),
    q33(q.q33), q34(q.q34),
    q44(q.q44) {}

  Quadric4f &operator+=(const Quadric4f &q) {
    q11 += q.q11; q12 += q.q12; q13 += q.q13; q14 += q.q14;
    q22 += q.q22; q23 += q.q23; q24 += q.q24;
    q33 += q.q33; q34 += q.q34;
    q44 += q.q44;
    return *this;
  }

  Quadric4f operator+(const Quadric4f &q) const {
    return Quadric4f(q11 + q.q11, q12 + q.q12, q13 + q.q13, q14 + q.q14,
                     q22 + q.q22, q23 + q.q23, q24 + q.q24,
                     q33 + q.q33, q34 + q.q34,
                     q44 + q.q44);
  }

  real apply(const Vector3f &v) const {
    return q11*v.x*v.x + 2*q12*v.x*v.y + 2*q13*v.x*v.z + 2*q14*v.x + q22*v.y*v.y + 2*q23*v.y*v.z + 2*q24*v.y + q33*v.z*v.z + 2*q34*v.z + q44;
  }
};

#endif
