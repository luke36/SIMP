#include "mesh.hpp"
#include "kd.hpp"
#include "heap.hpp"
#include <cassert>
#include <string>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <set>
#include <map>

static void remember_pair(Point *a, Point *b,
                          std::set<std::pair<Point *, Point *>> &done) {
  if (a < b) {
    done.insert({a, b});
  } else {
    done.insert({b, a});
  }
}

static bool ask_pair(Point *a, Point *b,
                     std::set<std::pair<Point *, Point *>> &done) {
  if (a < b) {
    return done.count({a, b});
  } else {
    return done.count({b, a});
  }
}

Point &Point::merge(Point *p, const Vector3f &pos, Heap &pairs) {
  x = pos.x;
  y = pos.y;
  z = pos.z;
  Q += p->Q;

  p->fa = this;

  std::set<std::pair<Point *, Point *>> changed;
  ps.splice(ps.end(), p->ps);
  for (auto pr : ps) {
    if (pr->valid) {
      pr->updateVertex(p, this, pairs);
      if (pr->valid) {
        if (ask_pair(pr->p1, pr->p2, changed)) {
          pr->valid = false;
        } else {
          remember_pair(pr->p1, pr->p2, changed);
        }
      }
    }
  }
  return *this;
}

Face::Face(Point *p1, Point *p2, Point *p3)
  : p1(p1), p2(p2), p3(p3) {
  Vector3f &&norm = cross(*p2 - *p1, *p3 - *p1);
  norm.normalize();
  // check nan if the face is degenerate
  if (!std::isnormal(norm.x) || !std::isnormal(norm.y) || !std::isnormal(norm.z)) {
    return;
  }
  real
    a = norm.x,
    b = norm.y,
    c = norm.z,
    d = -dot(*p1, norm);
  Quadric4f Kp = Quadric4f(a*a, a*b, a*c, a*d,
                                b*b, b*c, b*d,
                                     c*c, c*d,
                                          d*d);
  p1->Q += Kp;
  p2->Q += Kp;
  p3->Q += Kp;
}

static void compute_optimal(const Vector3f &v1, const Vector3f &v2, const Quadric4f &Q,
                            Vector3f &opt, real &error) {
  real m[4][4] = {{Q.q11, Q.q12, Q.q13, Q.q14},
                  {Q.q12, Q.q22, Q.q23, Q.q24},
                  {Q.q13, Q.q23, Q.q33, Q.q34},
                  {0.0f,  0.0f,  0.0f,  1.0f}};
  Matrix4f d(m);
  std::optional<Matrix4f> &&inv = inverse(d);
  if (inv) {
    opt = Vector3f(inv.value().m[0][3],
                   inv.value().m[1][3],
                   inv.value().m[2][3]);
    error = Q.apply(opt);
    return;
  }

  Vector3f mid = (v1 + v2) / 2;
  real
    error_1 = Q.apply(v1),
    error_2 = Q.apply(v2),
    error_mid = Q.apply(mid);
  if (error_1 < error_2) {
    if (error_1 < error_mid) {
      opt = v1;
      error = error_1;
    } else {
      opt = mid;
      error = error_mid;
    }
  } else {
    if (error_2 < error_mid) {
      opt = v2;
      error = error_2;
    } else {
      opt = mid;
      error = error_mid;
    }
  }
}

Pair::Pair(Point *x, Point *y)
  : p1(x), p2(y), valid(true) {
  x->ps.emplace_back(this);
  y->ps.emplace_back(this);
  compute_optimal(*x, *y, x->Q + y->Q, opt, error);
}

void Pair::updateVertex(Point *x, Point *y, Heap &ps) {
  if (p1 == x) {
    p1 = y;
  }
  if (p2 == x) {
    p2 = y;
  }
  if (p1 != p2) {
    compute_optimal(*p1, *p2, p1->Q + p2->Q, opt, error);
    ps.update(this);
  } else {
    valid = false;
    ps.erase(this);
  }
}

static void add_pair(Point *a, Point *b,
                     std::set<std::pair<Point *, Point *>> &selected) {
  if (a > b) {
    std::swap(a, b);
  }
  if (selected.find({a, b}) == selected.end()) {
    selected.insert({a, b});
  }
}

Mesh &Mesh::simplify(std::function<void (Mesh &, real ratio)> k,
                     std::vector<real> percentage, real epsilon) {
  std::cerr << "initializing ..." << std::endl;

  // add edges
  std::set<std::pair<Point *, Point *>> selected;
  for (auto &f : faces) {
    add_pair(f.p1, f.p2, selected);
    add_pair(f.p2, f.p3, selected);
    add_pair(f.p3, f.p1, selected);
  }

  // add close vertices
  if (epsilon > 0) {
    std::vector<Point *> pts;
    for (auto &p : points) {
      pts.push_back(&p);
    }

    KDTree kdt = buildKDTree(pts);
    pts.clear();

    for (auto &p : points) {
      kdt.radiusSearch(p, epsilon, pts);
      for (auto &q : pts) {
        if (&p != q) {
          add_pair(&p, q, selected);
        }
      }
      pts.clear();
    }
  }

  std::vector<Pair *> pre_heap;
  for (auto &pp : selected) {
    pre_heap.push_back(new Pair(pp.first, pp.second));
  }
  Heap pairs(std::move(pre_heap));

  std::cerr << "initialization end." << std::endl;

  std::sort(percentage.begin(), percentage.end());
  std::vector<Pair *> removed;
  size_t n_points = points.size(), n = n_points;
  do {
    std::cerr << "next percentage: " << percentage.back() << std::endl;
    while (n > percentage.back() * n_points) {
      auto least = pairs.top();
      if (least->valid) {
        least->p1->merge(least->p2, least->opt, pairs);
        n -= 1;
      } else {
        pairs.erase(least);
      }
      removed.emplace_back(least);
    }
    k(*this, percentage.back());
    percentage.pop_back();
  } while (!percentage.empty());

  for  (auto &p : removed) {
    delete p;
  }

  return *this;
}

Mesh::Mesh(std::istream &is) {
  std::string op;
  std::string garbage;
  while (is >> op) {
    if (op == "v") {
      real x, y, z;
      is >> x >> y >> z;
      points.emplace_back(x, y, z);
    } else if (op == "f") {
      std::vector<size_t> vs;
      size_t v;
      is >> v;
      if (is.peek() == '/') {
        char slash;
        size_t vt;
        is >> slash >> vt;
        do {
          vs.emplace_back(v);
          is >> v >> slash >> vt;
        } while (is.good());
      } else {
        do {
          vs.emplace_back(v);
          is >> v;
        } while (is.good());
      }
      is.clear();
      for (size_t i = 1; i < vs.size() - 1; i++) {
        faces.emplace_front(&points[vs[0] - 1], &points[vs[i] - 1], &points[vs[i+1] - 1]);
      }
    } else {
      std::getline(is, garbage);
    }
  }
}

static std::tuple<Point *, Point *, Point *> sort3(Point *a, Point *b, Point *c) {
  if (c < b) {
    std::swap(c, b);
  }
  if (b < a) {
    std::swap(b, a);
  }
  if (c < b) {
    std::swap(c, b);
  }
  return {a, b, c};
}

void Mesh::dump(std::ostream &os, int precision) {
  os << std::setprecision(precision);

  std::map<const Point *, size_t> number;
  std::set<std::tuple<Point *, Point *, Point *>> face;
  size_t n = 0;
  for (auto &p : points) {
    if (p.useful()) {
      n += 1;
      number[&p] = n;
      os << "v" << ' '
         << p.x << ' '
         << p.y << ' '
         << p.z << '\n';
    }
  }
  for (auto &f : faces) {
    f.p1 = (f.p1)->repr();
    f.p2 = (f.p2)->repr();
    f.p3 = (f.p3)->repr();
    auto sorted = sort3(f.p1, f.p2, f.p3);
    if (f.p1 != f.p2 && f.p2 != f.p3 && f.p3 != f.p1 && !face.count(sorted)) {
      face.insert(sorted);
      os << "f" << ' '
         << number[f.p1] << ' '
         << number[f.p2] << ' '
         << number[f.p3] << '\n';
    }
  }
}
