#ifndef KD_HPP
#define KD_HPP

#include <memory>
#include "mesh.hpp"

using std::shared_ptr;
using std::make_shared;

enum class KDAxis {KDX, KDY, KDZ};

class KDTreeBase {
public:
  virtual void
  radiusSearch(const Vector3f &ref, real r, std::vector<Point *> &ret) = 0;
  virtual ~KDTreeBase() = default;
};

class KDTree {
  shared_ptr<KDTreeBase> t;
public:
  KDTree(shared_ptr<KDTreeBase> t) : t(t) {}
  void radiusSearch(const Vector3f &ref, real r, std::vector<Point *> &ret) {
    t->radiusSearch(ref, r, ret);
  }
};

class KDTreeLeaf : public KDTreeBase {
  Point *p;
public:
  KDTreeLeaf(Point *p) : p(p) {}
  virtual void
  radiusSearch(const Vector3f &ref, real r, std::vector<Point *> &ret) override {
    if (distance(ref, *p) <= r) {
      ret.push_back(p);
    }
  }
  ~KDTreeLeaf() = default;
};

class KDTreeNode : public KDTreeBase {
  // <= and >
  KDTree low;
  KDTree hig;
  KDAxis axis;
  real coord;
public:
  KDTreeNode(KDTree low, KDTree hig, KDAxis axis, real coord)
    : low(low), hig(hig), axis(axis), coord(coord) {}
  virtual void
  radiusSearch(const Vector3f &ref, real r, std::vector<Point *> &ret) override {
    real refc;
    switch (axis) {
    case KDAxis::KDX:
      refc = ref.x;
      break;
    case KDAxis::KDY:
      refc = ref.y;
      break;
    case KDAxis::KDZ:
      refc = ref.z;
      break;
    }
    if (refc <= coord) {
      low.radiusSearch(ref, r, ret);
      if (coord - refc < r) {
        hig.radiusSearch(ref, r, ret);
      }
    } else {
      hig.radiusSearch(ref, r, ret);
      if (refc - coord <= r) {
        low.radiusSearch(ref, r, ret);
      }
    }
  }
  ~KDTreeNode() = default;
};

// pts can be modified
KDTree buildKDTree(std::vector<Point *> &pts);

#endif
