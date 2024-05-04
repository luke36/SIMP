#include "mesh.hpp"
#include <cassert>
#include <string>
#include <iomanip>
#include <utility>
#include <algorithm>
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

Point &Point::merge(Point *p, const Vector3f &pos, SortedPairs &pairs) {
  x = pos.x;
  y = pos.y;
  z = pos.z;
  Q += p->Q;

  for (auto &f : p->fs) {
    f->updateVertex(p, this);
  }
  fs.splice(fs.end(), p->fs);

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
  p1->fs.emplace_back(this);
  p2->fs.emplace_back(this);
  p3->fs.emplace_back(this);
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

Face &Face::updateVertex(Point *x, Point *y) {
  if (p1 == x) {
    p1 = y;
  }
  if (p2 == x) {
    p2 = y;
  }
  if (p3 == x) {
    p3 = y;
  }
  return *this;
}

// todo
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

void Pair::updateVertex(Point *x, Point *y, SortedPairs &ps) {
  ps.erase(this);
  if (p1 == x) {
    p1 = y;
  }
  if (p2 == x) {
    p2 = y;
  }
  if (p1 != p2) {
    compute_optimal(*p1, *p2, p1->Q + p2->Q, opt, error);
    ps.insert(this);
  } else {
    valid = false;
  }
}

static void add_pair(Point *a, Point *b,
                     SortedPairs &ps,
                     std::set<std::pair<Point *, Point *>> &selected) {
  if (a > b) {
    auto t = a;
    a = b;
    b = t;
  }
  if (selected.find({a, b}) == selected.end()) {
    selected.insert({a, b});
    Pair *p = new Pair(a, b);
    ps.insert(p);
  }
}

Mesh &Mesh::simplify(real percentage, real epsilon) {
  SortedPairs pairs;
  // add edges
  std::set<std::pair<Point *, Point *>> selected;
  for (auto &f : faces) {
    add_pair(f.p1, f.p2, pairs, selected);
    add_pair(f.p2, f.p3, pairs, selected);
    add_pair(f.p3, f.p1, pairs, selected);
  }

  // add close vertices
  if (epsilon > 0) {
    std::vector<Point *> x_sort;
    for (auto &p : points) {
      x_sort.push_back(&p);
    }
    std::sort(x_sort.begin(), x_sort.end(),
              [](Point *a, Point *b){ return a->x < b->x; });
    for (auto i = x_sort.begin(); i != x_sort.end(); ++i) {
      for (auto j = i + 1; j != x_sort.end() && (*j)->x - (*i)->x < epsilon; ++j) {
        if (distance(**i, **j) < epsilon) {
          add_pair(*i, *j, pairs, selected);
          // std::cerr << distance(**i, **j) << ' ' << (*i)->x << std::endl;
        }
      }
    }
  }

  assert(percentage >= 0);
  std::vector<Pair *> removed;
  size_t n_points = points.size();
  size_t target = n_points * percentage;
  while (n_points > target) {
    auto least = pairs.begin();
    auto least_p = *least;
    if (least_p->valid) {
      least_p->p1->merge(least_p->p2, least_p->opt, pairs);
      n_points -= 1;
    } else {
      pairs.erase(least);
    }
    removed.emplace_back(least_p);
    // delete least_p;
  }

  for (auto &p : pairs) {
    delete p;
  }
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

void Mesh::dump(std::ostream &os, int precision) const {
  os << std::setprecision(precision);

  std::map<const Point *, size_t> number;
  size_t n = 0;
  for (const auto &p : points) {
    if (p.useful()) {
      n += 1;
      number[&p] = n;
      os << "v" << ' '
         << p.x << ' '
         << p.y << ' '
         << p.z << '\n';
    }
  }
  for (const auto &f : faces) {
    if (f.p1 != f.p2 && f.p2 != f.p3 && f.p3 != f.p1) {
      os << "f" << ' '
         << number[f.p1] << ' '
         << number[f.p2] << ' '
         << number[f.p3] << '\n';
    }
  }
}
