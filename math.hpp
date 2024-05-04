#ifndef MATH_HPP
#define MATH_HPP

#include "real.hpp"
#include <cassert>
#include <cmath>
#include <optional>

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

// taken from pbrt

static inline real diff_prod(real a, real b, real c, real d) {
  real cd = c * d;
  real res = std::fma(a, b, -cd);
  real error = std::fma(c, d, -cd);
  return res + error;
}

inline Vector3f cross(const Vector3f &v1, const Vector3f &v2) {
  return Vector3f(diff_prod(v1.y, v2.z, v1.z, v2.y),
                  diff_prod(v1.z, v2.x, v1.x, v2.z),
                  diff_prod(v1.x, v2.y, v1.y, v2.x));
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

class Matrix4f {
public:
  real m[4][4];

  Matrix4f(real m_[4][4]) {
    for (int i = 0; i < 4; i++) {
      for (int j = 0; j < 4; j++) {
        m[i][j] = m_[i][j];
      }
    }
  }
};

// taken from pbrt

static inline real scf(std::pair<real, real> cf) {
  return cf.first + cf.second;
}

static inline std::pair<real, real> TwoProd(real a, real b) {
    real ab = a * b;
    return {ab, std::fma(a, b, -ab)};
}

static inline std::pair<real, real> TwoSum(real a, real b) {
    real s = a + b, delta = s - a;
    return {s, (a - (s - delta)) + (b - delta)};
}

template <typename real>
static inline std::pair<real, real> InnerProduct(real a, real b) {
    return TwoProd(a, b);
}

// Accurate dot products with FMA: Graillat et al.,
// https://www-pequan.lip6.fr/~graillat/papers/posterRNC7.pdf
//
// Accurate summation, dot product and polynomial evaluation in complex
// floating point arithmetic, Graillat and Menissier-Morain.
template <typename real, typename... T>
static inline std::pair<real, real> InnerProduct(real a, real b, T... terms) {
    auto ab = TwoProd(a, b);
    auto tp = InnerProduct(terms...);
    auto sum = TwoSum(ab.first, tp.first);
    return {sum.first, ab.second + (tp.second + sum.second)};
}

inline std::optional<Matrix4f> inverse(const Matrix4f &m) {
  // Via: https://github.com/google/ion/blob/master/ion/math/matrixutils.cc,
  // (c) Google, Apache license.

  // For 4x4 do not compute the adjugate as the transpose of the cofactor
  // matrix, because this results in extra work. Several calculations can be
  // shared across the sub-determinants.
  //
  // This approach is explained in David Eberly's Geometric Tools book,
  // excerpted here:
  //   http://www.geometrictools.com/Documentation/LaplaceExpansionTheorem.pdf
  real s0 = diff_prod(m.m[0][0], m.m[1][1], m.m[1][0], m.m[0][1]);
  real s1 = diff_prod(m.m[0][0], m.m[1][2], m.m[1][0], m.m[0][2]);
  real s2 = diff_prod(m.m[0][0], m.m[1][3], m.m[1][0], m.m[0][3]);

  real s3 = diff_prod(m.m[0][1], m.m[1][2], m.m[1][1], m.m[0][2]);
  real s4 = diff_prod(m.m[0][1], m.m[1][3], m.m[1][1], m.m[0][3]);
  real s5 = diff_prod(m.m[0][2], m.m[1][3], m.m[1][2], m.m[0][3]);

  real c0 = diff_prod(m.m[2][0], m.m[3][1], m.m[3][0], m.m[2][1]);
  real c1 = diff_prod(m.m[2][0], m.m[3][2], m.m[3][0], m.m[2][2]);
  real c2 = diff_prod(m.m[2][0], m.m[3][3], m.m[3][0], m.m[2][3]);

  real c3 = diff_prod(m.m[2][1], m.m[3][2], m.m[3][1], m.m[2][2]);
  real c4 = diff_prod(m.m[2][1], m.m[3][3], m.m[3][1], m.m[2][3]);
  real c5 = diff_prod(m.m[2][2], m.m[3][3], m.m[3][2], m.m[2][3]);

  real determinant = scf(InnerProduct(s0, c5, -s1, c4, s2, c3, s3, c2, s5, c0, -s4, c1));
  if (determinant == 0)
    return {};
  real s = 1 / determinant;

  real inv[4][4] = {{s * scf(InnerProduct(m.m[1][1], c5, m.m[1][3], c3, -m.m[1][2], c4)),
                     s * scf(InnerProduct(-m.m[0][1], c5, m.m[0][2], c4, -m.m[0][3], c3)),
                     s * scf(InnerProduct(m.m[3][1], s5, m.m[3][3], s3, -m.m[3][2], s4)),
                     s * scf(InnerProduct(-m.m[2][1], s5, m.m[2][2], s4, -m.m[2][3], s3))},

                    {s * scf(InnerProduct(-m.m[1][0], c5, m.m[1][2], c2, -m.m[1][3], c1)),
                     s * scf(InnerProduct(m.m[0][0], c5, m.m[0][3], c1, -m.m[0][2], c2)),
                     s * scf(InnerProduct(-m.m[3][0], s5, m.m[3][2], s2, -m.m[3][3], s1)),
                     s * scf(InnerProduct(m.m[2][0], s5, m.m[2][3], s1, -m.m[2][2], s2))},

                    {s * scf(InnerProduct(m.m[1][0], c4, m.m[1][3], c0, -m.m[1][1], c2)),
                     s * scf(InnerProduct(-m.m[0][0], c4, m.m[0][1], c2, -m.m[0][3], c0)),
                     s * scf(InnerProduct(m.m[3][0], s4, m.m[3][3], s0, -m.m[3][1], s2)),
                     s * scf(InnerProduct(-m.m[2][0], s4, m.m[2][1], s2, -m.m[2][3], s0))},

                    {s * scf(InnerProduct(-m.m[1][0], c3, m.m[1][1], c1, -m.m[1][2], c0)),
                     s * scf(InnerProduct(m.m[0][0], c3, m.m[0][2], c0, -m.m[0][1], c1)),
                     s * scf(InnerProduct(-m.m[3][0], s3, m.m[3][1], s1, -m.m[3][2], s0)),
                     s * scf(InnerProduct(m.m[2][0], s3, m.m[2][2], s0, -m.m[2][1], s1))}};

  return Matrix4f(inv);
}

#endif
