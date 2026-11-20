/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "gxtb_revd4.h"
#include "gxtb_pairs.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <new>
#include <stdexcept>
#include <vector>

namespace {
using gxtb::difference;
using gxtb::dot;
using gxtb::Gradient;
using gxtb::pi;
using gxtb::Vec;
using Spectrum = std::array<double, 23>;
constexpr double cn_radius = 30.0, pair_radius = 50.0, atm_radius = 25.0;
struct Element {
  int first, count;
  double radius, r4r2;
  double charge[4];
};
struct Reference {
  double cn;
  int gaussians;
  Spectrum spectrum;
};
#include "gxtb_revd4_data.inc"

std::array<double, 2> charge_factor(int z, double q) {
  const auto &p = elements[z].charge;
  const double t = std::tanh(p[2] * q + p[3]);
  return {p[0] + p[1] * t, p[1] * p[2] * (1.0 - t * t)};
}

struct Response {
  Spectrum alpha{}, derivative{};
};
Response response(int z, int partner, double cn) {
  const auto &element = elements[z];
  const int hi = std::max(z, partner), lo = std::min(z, partner);
  const double width = pair_widths[hi * (hi + 1) / 2 + lo];
  Response result;
  double norm = 0.0, dnorm = 0.0;
  const Reference *last = references + element.first;
  for (int j = 0; j < element.count; ++j) {
    const auto &ref = references[element.first + j];
    if (ref.cn > last->cn)
      last = &ref;
    const double delta = cn - ref.cn;
    double weight = 0.0, slope = 0.0;
    for (int m = 1; m <= ref.gaussians; ++m) {
      const double term = std::exp(-m * width * delta * delta);
      weight += term;
      slope -= 2 * m * width * delta * term;
    }
    norm += weight;
    dnorm += slope;
    for (size_t k = 0; k < ref.spectrum.size(); ++k) {
      result.alpha[k] += weight * ref.spectrum[k];
      result.derivative[k] += slope * ref.spectrum[k];
    }
  }
  if (norm < std::numeric_limits<double>::min()) {
    result = {};
    for (int j = 0; j < element.count; ++j) {
      const auto &ref = references[element.first + j];
      if (std::abs(ref.cn - last->cn) < 1.e-12)
        for (size_t k = 0; k < ref.spectrum.size(); ++k)
          result.alpha[k] += ref.spectrum[k];
    }
    return result;
  }
  for (size_t k = 0; k < result.alpha.size(); ++k) {
    result.alpha[k] /= norm;
    result.derivative[k] =
        (result.derivative[k] - result.alpha[k] * dnorm) / norm;
  }
  return result;
}

struct Coefficient {
  double k = 0, di = 0, dj = 0;
  double neutral = 0, rho = 0, sum_damping = 0;
  int i, j;
};
struct Edge : gxtb::Pair {
  double dcn = 0, inverse = 0, log_slope = 0;
  Coefficient *coefficient = nullptr;
};

double cube(double x) { return x * x * x; }
double damping6(const Edge &edge) {
  const double t = cube(edge.inverse);
  return t * t;
}

class Geometry {
public:
  std::vector<int> z;
  std::vector<double> atm_energy;
  Gradient atm_gradient;
  Geometry(const Geometry &) = delete;
  Geometry &operator=(const Geometry &) = delete;
  Geometry(int n, const int *numbers, const double *xyz, const double *cell,
           const int *pbc)
      : z(n), atm_energy(n), atm_gradient(n) {
    std::vector<double> cn(n);
    for (int i = 0; i < n; ++i) {
      if (numbers[i] < 1 || numbers[i] > 103)
        throw std::invalid_argument("D4Srev requires Z=1..103");
      z[i] = numbers[i] - 1;
    }
    gxtb::for_each_pair(
        n, xyz, cell, pbc, pair_radius, [&](const gxtb::Pair &pair) {
          Edge edge{pair};
          if (edge.r <= cn_radius) {
            const double beta =
                7.5 / (elements[z[edge.i]].radius + elements[z[edge.j]].radius);
            const double u = beta * edge.r - 7.5;
            const double count = 0.5 * std::erfc(u);
            cn[edge.i] += count;
            cn[edge.j] += count;
            edge.dcn = -beta / std::sqrt(pi) * std::exp(-u * u);
          }
          edges.push_back(edge);
        });
    std::vector<std::vector<const Edge *>> forward(n);
    std::map<std::pair<int, int>, Response> spectra;
    auto spectrum = [&](int atom, int partner) -> const Response & {
      auto entry = spectra.try_emplace({atom, z[partner]});
      if (entry.second)
        entry.first->second = response(z[atom], z[partner], cn[atom]);
      return entry.first->second;
    };
    for (Edge &edge : edges) {
      const int i = std::min(edge.i, edge.j), j = std::max(edge.i, edge.j);
      auto entry = coefficients.try_emplace({i, j});
      Coefficient &coeff = entry.first->second;
      if (entry.second) {
        const auto &x = spectrum(i, j), &y = spectrum(j, i);
        for (size_t k = 0; k < x.alpha.size(); ++k) {
          coeff.k += x.alpha[k] * y.alpha[k];
          coeff.di += x.derivative[k] * y.alpha[k];
          coeff.dj += x.alpha[k] * y.derivative[k];
        }
        coeff.i = i;
        coeff.j = j;
        coeff.rho = std::sqrt(3 * elements[z[i]].r4r2 * elements[z[j]].r4r2);
        coeff.neutral = std::sqrt(std::abs(coeff.k * charge_factor(z[i], 0)[0] *
                                           charge_factor(z[j], 0)[0]));
      }
      damp(edge, coeff);
      coeff.sum_damping += damping6(edge);
      if (edge.r <= atm_radius)
        forward[edge.i].push_back(&edge);
    }
    // Each periodic triple has exactly one lowest site. Both of its
    // forward edges share that site, even when basis indices repeat.
    std::vector<double> dcn(n);
    for (const auto &neighbors : forward) {
      for (size_t j = 1; j < neighbors.size(); ++j)
        for (size_t k = 0; k < j; ++k) {
          const Edge &first = *neighbors[j], &second = *neighbors[k];
          const Vec d = difference(second.d, first.d);
          Edge third{{first.j, second.j, d, std::sqrt(dot(d, d))}};
          if (third.r > atm_radius)
            continue;
          damp(third, coefficients.at(std::minmax(third.i, third.j)));
          add_triple({&first, &second, &third}, dcn);
        }
    }
    add_cn_gradient(atm_gradient, dcn);
  }

  void electronic(const double *q, double *e, double *v) const {
    std::vector<std::array<double, 2>> f(z.size());
    for (size_t i = 0; i < z.size(); ++i)
      f[i] = charge_factor(z[i], q[i]);
    std::fill_n(e, z.size(), 0.0);
    std::fill_n(v, z.size(), 0.0);
    for (const auto &entry : coefficients) {
      const auto &c = entry.second;
      const double w = -c.k * c.sum_damping;
      const double half = 0.5 * w * f[c.i][0] * f[c.j][0];
      e[c.i] += half;
      e[c.j] += half;
      v[c.i] += w * f[c.i][1] * f[c.j][0];
      v[c.j] += w * f[c.j][1] * f[c.i][0];
    }
  }
  void gradient(const double *q, double *g, double *s) const {
    Gradient result(z.size());
    std::vector<double> f(z.size()), dcn(z.size());
    for (size_t i = 0; i < z.size(); ++i)
      f[i] = charge_factor(z[i], q[i])[0];
    for (const Edge &edge : edges) {
      const auto &c = *edge.coefficient;
      const double factor = -f[c.i] * f[c.j] * damping6(edge);
      result.add(edge, 6 * factor * c.k * edge.log_slope);
      dcn[c.i] += factor * c.di;
      dcn[c.j] += factor * c.dj;
    }
    add_cn_gradient(result, dcn);
    result.output(g, s);
  }

private:
  std::vector<Edge> edges;
  std::map<std::pair<int, int>, Coefficient> coefficients;
  static void damp(Edge &edge, Coefficient &c) {
    constexpr double a2 = 8.5302496868, a3 = 0.5395219242, a4 = 0.654;
    const double x = a3 * (edge.r - a4 * c.rho);
    edge.inverse = 1 / (edge.r + 0.5 * a2 * std::erfc(x));
    edge.log_slope =
        -edge.inverse * (1 - a2 * a3 / std::sqrt(pi) * std::exp(-x * x));
    edge.coefficient = &c;
  }
  void add_cn_gradient(Gradient &g, const std::vector<double> &dcn) const {
    for (const auto &edge : edges)
      if (edge.dcn != 0)
        g.add(edge, edge.dcn * (dcn[edge.i] + dcn[edge.j]));
  }
  void add_triple(const std::array<const Edge *, 3> &sides,
                  std::vector<double> &dcn) {
    double scale = 1.826410691;
    Vec r{}, cosine{}, radial{};
    for (int a = 0; a < 3; ++a) {
      const auto &edge = *sides[a];
      scale *= edge.coefficient->neutral * cube(edge.inverse);
      r[a] = edge.r;
    }
    if (scale == 0)
      return;
    for (int a = 0; a < 3; ++a) {
      const int b = (a + 1) % 3, c = (a + 2) % 3;
      cosine[a] = (r[b] * r[b] + r[c] * r[c] - r[a] * r[a]) / (2 * r[b] * r[c]);
    }
    const double energy = scale * (1 + 3 * cosine[0] * cosine[1] * cosine[2]);
    // Differentiate each cosine in the law of cosines with respect to
    // its opposite side and its two adjacent sides.
    for (int a = 0; a < 3; ++a) {
      const int b = (a + 1) % 3, c = (a + 2) % 3;
      const double angular = 3 * scale * cosine[b] * cosine[c];
      radial[a] -= angular * r[a] / (r[b] * r[c]);
      radial[b] += angular * (1 / r[c] - cosine[a] / r[b]);
      radial[c] += angular * (1 / r[b] - cosine[a] / r[c]);
    }
    for (int a = 0; a < 3; ++a) {
      const auto &edge = *sides[a];
      const auto &c = *edge.coefficient;
      atm_gradient.add(edge, radial[a] + 3 * energy * edge.log_slope);
      const double factor = 0.5 * energy / c.k;
      dcn[c.i] += factor * c.di;
      dcn[c.j] += factor * c.dj;
    }
    atm_energy[sides[0]->i] += energy / 3;
    atm_energy[sides[0]->j] += energy / 3;
    atm_energy[sides[1]->j] += energy / 3;
  }
};

void *create(int n, const int *z, const double *xyz, const double *cell,
             const int *pbc, char *msg, size_t len) {
  try {
    if (n < 1)
      throw std::invalid_argument("Empty D4Srev structure");
    return new Geometry(n, z, xyz, cell, pbc);
  } catch (const std::exception &failure) {
    if (len)
      std::snprintf(msg, len, "%s", failure.what());
    return nullptr;
  }
}
void destroy(void *p) { delete static_cast<Geometry *>(p); }
void scf(const void *p, const double *q, double *e, double *v) {
  const auto &geom = *static_cast<const Geometry *>(p);
  try {
    geom.electronic(q, e, v);
  } catch (const std::bad_alloc &) {
    std::fill_n(e, geom.z.size(), std::numeric_limits<double>::quiet_NaN());
    std::fill_n(v, geom.z.size(), std::numeric_limits<double>::quiet_NaN());
  }
}
void derivatives(const void *p, const double *q, double *g, double *s) {
  const auto &geom = *static_cast<const Geometry *>(p);
  try {
    geom.gradient(q, g, s);
  } catch (const std::bad_alloc &) {
    std::fill_n(g, 3 * geom.z.size(), std::numeric_limits<double>::quiet_NaN());
    std::fill_n(s, 9, std::numeric_limits<double>::quiet_NaN());
  }
}
void atm(const void *p, double *e, double *g, double *s) {
  const auto &geom = *static_cast<const Geometry *>(p);
  std::copy(geom.atm_energy.begin(), geom.atm_energy.end(), e);
  geom.atm_gradient.output(g, s);
}
} // namespace

extern "C" const cp2k_gxtb_d4_v1 *cp2k_gxtb_native_d4(void) {
  static const cp2k_gxtb_d4_v1 methods{create, destroy, scf, derivatives, atm};
  return &methods;
}
