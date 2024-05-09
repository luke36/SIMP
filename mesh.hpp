#ifndef MESH_HPP
#define MESH_HPP

#include "math.hpp"
#include <functional>
#include <vector>
#include <list>
#include <iostream>

class Face;
class Pair;
class Heap;

class Point : public Vector3f {
  friend class Face;
  friend class Pair;
private:
  Quadric4f Q;
  std::list<Pair *> ps;
  Point *fa;
public:
  Point(real x, real y, real z) : Vector3f(x, y, z), fa(nullptr) {}
  Point &merge(Point *p, const Vector3f &pos, Heap &pairs);
  bool useful() const { return fa == nullptr; }
  Point *repr() {
    if (fa == nullptr) {
      return this;
    } else {
      fa = fa->repr();
      return fa;
    }
  }
};

class Face {
  friend class Mesh;
private:
  Point *p1;
  Point *p2;
  Point *p3;
public:
  Face(Point *p1, Point *p2, Point *p3);
};

class Pair {
  friend class Point;
  friend class Mesh;
  friend class Heap;
private:
  size_t id; // index in the heap
  Point *p1;
  Point *p2;
  Vector3f opt;
  real error;
public:
  bool valid;
  Pair(Point *p1, Point *p2);
  void updateVertex(Point *x, Point *y, Heap &ps);
  bool degenerate() const { return p1 == p2; }
};

class Mesh {
private:
  std::vector<Point> points;
  std::list<Face> faces;
public:
  Mesh(std::istream &is);
  void dump(std::ostream &os, int precision = 8);
  Mesh &simplify(std::function<void (Mesh &, real ratio)> k,
                 std::vector<real> percentage, real epsilon = 0);
};

#endif
