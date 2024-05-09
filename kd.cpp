#include "kd.hpp"

static real variance(std::vector<Point *> &pts, size_t b, size_t e, KDAxis a) {
  real mean = 0.0;
  real var = 0.0;
  switch (a) {
  case KDAxis::KDX:
    for (size_t i = b; i != e; i++) {
      mean += pts[i]->x;
    }
    mean /= b - e + 1;
    for (size_t i = b; i != e; i++) {
      real diff = pts[i]->x - mean;
      var += diff * diff;
    }
    return var;
  case KDAxis::KDY:
    for (size_t i = b; i != e; i++) {
      mean += pts[i]->y;
    }
    mean /= b - e + 1;
    for (size_t i = b; i != e; i++) {
      real diff = pts[i]->y - mean;
      var += diff * diff;
    }
    return var;
  case KDAxis::KDZ:
    for (size_t i = b; i != e; i++) {
      mean += pts[i]->z;
    }
    mean /= b - e + 1;
    for (size_t i = b; i != e; i++) {
      real diff = pts[i]->z - mean;
      var += diff * diff;
    }
    return var;
  }
}

static size_t partite(std::vector<Point *> &pts, size_t b, size_t e,
                      KDAxis a) {
  real v;
  std::vector<Point *>::iterator h;
  size_t mid = (b + e) / 2;
  std::swap(pts[b], pts[mid]);
  switch (a) {
  case KDAxis::KDX:
    v = pts[b]->x;
    h = std::partition(pts.begin() + b + 1, pts.begin() + e + 1,
                       [v](Point *x){return x->x <= v;});
    break;
  case KDAxis::KDY:
    v = pts[b]->y;
    h = std::partition(pts.begin() + b + 1, pts.begin() + e + 1,
                       [v](Point *x){return x->y <= v;});
    break;
  case KDAxis::KDZ:
    v = pts[b]->z;
    h = std::partition(pts.begin() + b + 1, pts.begin() + e + 1,
                       [v](Point *x){return x->z <= v;});
    break;
  }
  std::swap(pts[b], *(h - 1));
  return h - 1 - pts.begin();
}

static KDAxis max(real &x, real &y, real &z) {
  if (x >= y) {
    if (x >= z) {
      return KDAxis::KDX;
    } else {
      return KDAxis::KDZ;
    }
  } else {
    if (y >= z) {
      return KDAxis::KDY;
    } else {
      return KDAxis::KDZ;
    }
  }
}

static real get_coord(Point *p, KDAxis a) {
  switch (a) {
  case KDAxis::KDX:
    return p->x;
  case KDAxis::KDY:
    return p->y;
  case KDAxis::KDZ:
    return p->z;
  }
}

static KDTree build_rec(std::vector<Point *> &pts, size_t b, size_t e) {
  if (b == e) {
    KDTreeBase *p = new KDTreeLeaf(pts[b]);
    return KDTree(shared_ptr<KDTreeBase>(p));
  } else {
    real
      var_x = variance(pts, b, e, KDAxis::KDX),
      var_y = variance(pts, b, e, KDAxis::KDY),
      var_z = variance(pts, b, e, KDAxis::KDZ);
    KDAxis a = max(var_x, var_y, var_z);
    a = KDAxis::KDX;
    size_t mid = partite(pts, b, e, a);
    if (mid == e) {
      mid -= 1;
    }
    KDTreeNode *p = new KDTreeNode(build_rec(pts, b, mid),
                                   build_rec(pts, mid + 1, e),
                                   a, get_coord(pts[mid], a));
    return KDTree(shared_ptr<KDTreeBase>(p));
  }
}

KDTree buildKDTree(std::vector<Point *> &pts) {
  return build_rec(pts, 0, pts.size() - 1);
}
