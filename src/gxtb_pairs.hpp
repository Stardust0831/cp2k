/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef CP2K_GXTB_PAIRS_HPP
#define CP2K_GXTB_PAIRS_HPP
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace gxtb {
using Vec = std::array<double, 3>;
using Tensor = std::array<double, 9>;
constexpr double pi = 3.14159265358979323846;
inline double dot(const Vec &a, const Vec &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
inline Vec cross(const Vec &a, const Vec &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}
inline Vec difference(const Vec &a, const Vec &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

struct Pair {
  int i, j;
  Vec d;
  double r;
};

class Lattice {
  std::array<Vec, 3> lattice{}, reciprocal{};
  std::array<int, 3> periodic{};

public:
  Lattice(const double *cell, const int *pbc) {
    std::copy_n(pbc, 3, periodic.begin());
    if (periodic == std::array<int, 3>{})
      return;
    for (int k = 0; k < 3; ++k)
      std::copy_n(cell + 3 * k, 3, lattice[k].begin());
    const double volume = dot(lattice[0], cross(lattice[1], lattice[2]));
    if (!std::isfinite(volume) || std::abs(volume) < 1.e-12)
      throw std::invalid_argument("Singular g-xTB cell");
    for (int k = 0; k < 3; ++k) {
      reciprocal[k] = cross(lattice[(k + 1) % 3], lattice[(k + 2) % 3]);
      for (double &v : reciprocal[k])
        v /= volume;
    }
  }
  template <class Visitor>
  void images(const Vec &delta, double cutoff, Visitor visit) const {
    std::array<int, 3> lo{}, hi{};
    for (int k = 0; k < 3; ++k) {
      if (!periodic[k])
        continue;
      const double center = -dot(delta, reciprocal[k]);
      const double span =
          cutoff * std::sqrt(dot(reciprocal[k], reciprocal[k])) + 1.e-12;
      if (!std::isfinite(center) || !std::isfinite(span) ||
          std::abs(center) + span >= std::numeric_limits<int>::max() - 1)
        throw std::invalid_argument("g-xTB lattice image index out of range");
      lo[k] = static_cast<int>(std::ceil(center - span));
      hi[k] = static_cast<int>(std::floor(center + span));
    }
    for (int a = lo[0]; a <= hi[0]; ++a)
      for (int b = lo[1]; b <= hi[1]; ++b)
        for (int c = lo[2]; c <= hi[2]; ++c) {
          Vec d = delta;
          for (int k = 0; k < 3; ++k)
            d[k] += a * lattice[0][k] + b * lattice[1][k] + c * lattice[2][k];
          const double r = std::sqrt(dot(d, d));
          if (r <= cutoff)
            visit(std::array<int, 3>{a, b, c}, d, r);
        }
  }

  // An actual image gives an upper bound on the minimum. Reciprocal-space
  // bounds then include every competitor, also for skew or unwrapped cells.
  std::vector<Pair> nearest(int i, int j, const Vec &delta) const {
    constexpr double band = 0.3; // squared-distance window in bohr^2
    double upper = std::numeric_limits<double>::infinity();
    if (i == j) {
      for (int k = 0; k < 3; ++k)
        if (periodic[k])
          upper = std::min(upper, dot(lattice[k], lattice[k]));
      if (!std::isfinite(upper))
        return {};
    } else {
      Vec trial = delta;
      for (int k = 0; k < 3; ++k)
        if (periodic[k]) {
          const double shift = std::round(dot(delta, reciprocal[k]));
          for (int axis = 0; axis < 3; ++axis)
            trial[axis] -= shift * lattice[k][axis];
        }
      upper = dot(trial, trial);
    }
    std::vector<Pair> found;
    double best = upper;
    images(delta, std::sqrt(upper + band) + 1.e-12,
           [&](const auto &shift, const Vec &d, double r) {
             if (i == j && shift == std::array<int, 3>{})
               return;
             if (r < 1.e-6)
               throw std::invalid_argument("Coincident g-xTB exchange centers");
             best = std::min(best, r * r);
             found.push_back({i, j, d, r});
           });
    found.erase(
        std::remove_if(found.begin(), found.end(),
                       [&](const Pair &p) { return p.r * p.r > best + band; }),
        found.end());
    if (found.empty())
      throw std::runtime_error("Missing nearest g-xTB exchange image");
    return found;
  }
};

inline Vec position(const double *xyz, int i) {
  return {xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2]};
}

template <class Visitor>
void for_each_pair(int n, const double *xyz, const double *cell, const int *pbc,
                   double cutoff, Visitor visit) {
  Lattice lattice(cell, pbc);
  for (size_t k = 0; k < 3 * static_cast<size_t>(n); ++k)
    if (!std::isfinite(xyz[k]))
      throw std::invalid_argument("Nonfinite g-xTB coordinate");
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      lattice.images(difference(position(xyz, j), position(xyz, i)), cutoff,
                     [&](const auto &shift, const Vec &d, double r) {
                       // Count the undirected edge once, including for self
                       // images.
                       if (shift < std::array<int, 3>{} ||
                           (shift == std::array<int, 3>{} && j <= i))
                         return;
                       if (r < 1.e-6)
                         throw std::invalid_argument(
                             "Coincident atoms or images in g-xTB");
                       visit(Pair{i, j, d, r});
                     });
}

struct Gradient {
  std::vector<Vec> atoms;
  Tensor strain{};
  explicit Gradient(size_t n) : atoms(n) {}
  void add(const Pair &edge, double slope) {
    for (int a = 0; a < 3; ++a) {
      const double g = slope * edge.d[a] / edge.r;
      if (edge.i != edge.j) {
        atoms[edge.i][a] -= g;
        atoms[edge.j][a] += g;
      }
      for (int b = 0; b < 3; ++b)
        strain[3 * b + a] += edge.d[b] * g;
    }
  }
  void output(double *g, double *s) const {
    for (size_t i = 0; i < atoms.size(); ++i)
      std::copy_n(atoms[i].data(), 3, g + 3 * i);
    std::copy(strain.begin(), strain.end(), s);
  }
};

} // namespace gxtb
#endif
