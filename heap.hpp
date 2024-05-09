#ifndef HEAP_H
#define HEAP_H

#include "mesh.hpp"

static inline size_t left(size_t x) {
  return x << 1;
}

static inline size_t right(size_t x) {
  return (x << 1) + 1;
}

static inline size_t mother(size_t x) {
  return x >> 1;
}

class Heap {
private:
  std::vector<Pair *> pts;

  bool le(size_t a, size_t b) {
    return pts[a]->error <= pts[b]->error;
  }

  void assign(size_t i, size_t j) {
    pts[i] = pts[j];
    pts[i]->id = i;
  }

  void up(size_t i) {
    Pair *v = pts[i];

    while (true) {
      if (mother(i) < 1) {
        break;
      } else if (pts[mother(i)]->error <= v->error) {
        break;
      } else {
        assign(i, mother(i));
        i = mother(i);
      }
    }
    pts[i] = v;
    v->id = i;
  }

  void down(size_t i) {
    Pair *v = pts[i];

    while (true) {
      if (left(i) >= pts.size()) {
        break;
      } else if (right(i) >= pts.size()) {
        if (v->error <= pts[left(i)]->error) {
          break;
        } else {
          assign(i, left(i));
          i = left(i);
        }
      } else if (v->error <= pts[left(i)]->error &&
                 v->error <= pts[right(i)]->error) {
        break;
      } else if (le(left(i), right(i))) {
        assign(i, left(i));
        i = left(i);
      } else {
        assign(i, right(i));
        i = right(i);
      }
    }
    pts[i] = v;
    v->id = i;
  }

public:
  Heap(std::vector<Pair *> &&pts_)
    : pts(std::move(pts_)) {
    pts.insert(pts.begin(), nullptr);
    for (size_t i = 1; i < pts.size(); i++) {
      pts[i]->id = i;
    }
    for (size_t i = pts.size(); i > 1; i--) {
      down(i - 1);
    }
  }

  void erase(Pair *p) {
    p->error = (-1.0) / 0.0; // - inf
    up(p->id);
    pop();
  }

  void insert(Pair *p) {
    p->id = pts.size();
    pts.push_back(p);
    up(p->id);
  }

  Pair *top() {
    return pts[1];
  }

  void pop() {
    assign(1, pts.size() - 1);
    pts.pop_back();
    down(1);
  }

  void update(Pair *p) {
    size_t id = p->id;
    if (id != 1 && !le(mother(id), id)) {
      up(id);
    } else {
      down(id);
    }
  }

  ~Heap() { // should I ?
    for (auto &p : pts) {
      delete p;
    }
  }
};

#endif
