/* SPDX-License-Identifier: GPL-2.0-or-later */
#if defined(__TBLITE)
#include "gxtb_integrals.hpp"
#include "gxtb_onecenter.hpp"
#include "gxtb_parameters.hpp"

namespace gxtb {
namespace {
constexpr int pbc[] = {1, 1, 1};
constexpr int indices[][2] = {{0, 0}, {0, 1}, {1, 1}, {0, 2}, {1, 2}, {2, 2}};
constexpr double damp[] = {0.3405910191, 0.1691310614, 0.0740343390, -0.02};
extern "C" void gxtb_ewald_alpha(const double *, double *);

struct Tensors {
  Vec sd{};
  Tensor dd{};
  std::array<double, 6> sq{};
  std::array<double, 18> dq{};
  std::array<double, 36> qq{};
  // Cartesian derivatives of 1/r through fourth order. Two equal
  // off-diagonal quadrupole components occupy one packed matrix entry.
  void real(const Vec &x, const double b[4][5]) {
    for (int a = 0; a < 3; ++a) {
      sd[a] += x[a] * b[0][1];
      for (int c = 0; c < 3; ++c)
        dd[3 * c + a] += (a == c) * b[1][1] - 3 * x[a] * x[c] * b[1][2];
    }
    for (int u = 0; u < 6; ++u) {
      const auto [a, bx] = indices[u];
      const int m = a == bx ? 1 : 2;
      sq[u] += m * (x[a] * x[bx] * b[1][2] - (a == bx) * b[1][1] / 3);
      for (int c = 0; c < 3; ++c)
        dq[3 * u + c] +=
            m * (5 * x[a] * x[bx] * x[c] * b[2][3] -
                 ((a == bx) * x[c] + (a == c) * x[bx] + (bx == c) * x[a]) *
                     b[2][2]);
      for (int v = 0; v < 6; ++v) {
        const auto [c, d] = indices[v];
        const double rr = (a == bx) * x[c] * x[d] + (a == c) * x[bx] * x[d] +
                          (a == d) * x[bx] * x[c] + (bx == c) * x[a] * x[d] +
                          (bx == d) * x[a] * x[c] + (c == d) * x[a] * x[bx];
        const int dd0 =
            (a == bx && c == d) + (a == c && bx == d) + (a == d && bx == c);
        qq[u + 6 * v] += m * (c == d ? 1 : 2) *
                         (35 * x[a] * x[bx] * x[c] * x[d] * b[3][4] -
                          5 * rr * b[3][3] + dd0 * b[3][2]);
      }
    }
  }
  void reciprocal(const Vec &k, double c, double s) {
    std::array<double, 6> kk;
    for (int u = 0; u < 6; ++u) {
      const auto [a, b] = indices[u];
      kk[u] = (a == b ? 1 : 2) * k[a] * k[b];
      sq[u] -= damp[1] * c * kk[u] / 3;
    }
    for (int a = 0; a < 3; ++a) {
      sd[a] += damp[0] * s * k[a];
      for (int b = 0; b < 3; ++b)
        dd[a + 3 * b] += damp[1] * c * k[a] * k[b];
      for (int u = 0; u < 6; ++u)
        dq[a + 3 * u] -= damp[2] * s * k[a] * kk[u] / 3;
    }
    for (int v = 0; v < 6; ++v)
      for (int u = 0; u < 6; ++u)
        qq[u + 6 * v] += damp[3] * c * kk[u] * kk[v] / 3;
  }
};
struct KPoint {
  Vec k;
  double charge, multi;
};
} // namespace

void coulomb_integrals(int n, const int *z, const double *xyz,
                       const double *cell, const Orbitals &orb, Matrices &m) {
  const size_t ns = orb.shell_atom.size(), n2 = size_t(n) * n;
  m.gamma.resize(ns * ns);
  m.tau.resize(ns * ns);
  m.sd.resize(3 * n2);
  m.dd.resize(9 * n2);
  m.sq.resize(6 * n2);
  m.dq.resize(18 * n2);
  m.qq.resize(36 * n2);
  const Lattice lattice(cell, pbc);
  const Vec a = position(cell, 0), b = position(cell, 1), c = position(cell, 2);
  const double det = dot(a, cross(b, c)), volume = std::abs(det);
  Tensor reciprocal{};
  for (int h = 0; h < 3; ++h) {
    const Vec normal =
        cross(position(cell, (h + 1) % 3), position(cell, (h + 2) % 3));
    for (int k = 0; k < 3; ++k)
      reciprocal[3 * h + k] = 2 * pi * normal[k] / det;
  }
  // The charge matrix omits the uniform G=0 mode in the tblite convention.
  // Its additive constant therefore depends on the charge splitting parameter.
  double alpha;
  gxtb_ewald_alpha(cell, &alpha);
  const double beta = std::sqrt(pi) / std::cbrt(volume);
  std::vector<KPoint> kpoints;
  Lattice(reciprocal.data(), pbc)
      .images({}, 14 * std::max(alpha, beta),
              [&](const auto &, const Vec &k, double g) {
                if (g == 0)
                  return;
                const double fac = 4 * pi / (volume * g * g);
                kpoints.push_back(
                    {k, fac * std::exp(-std::pow(g / (2 * alpha), 2)),
                     fac * std::exp(-std::pow(g / (2 * beta), 2))});
              });
  std::vector<size_t> offset(n + 1);
  Array hard(ns), dipscale(n);
  for (int i = 0; i < n; ++i) {
    const auto &e = basis_data(z[i]);
    offset[i + 1] = offset[i] + e.count;
    const double root = std::sqrt(m.cn[i] + 1.e-12) - 1.e-6;
    dipscale[i] = 1. + e.dipole_cn * root;
    for (int sh = 0; sh < e.count; ++sh)
      hard[offset[i] + sh] = e.hardness[sh] * (1. + e.hardness_cn * root);
  }
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j) {
      const Vec delta = difference(position(xyz, i), position(xyz, j));
      Tensors t;
      double ewald = 0.;
      const double radius = pair_radius(z[i], z[j]);
      const double rcut = std::max({7 / alpha, 7 / beta, radius + 16.});
      lattice.images(delta, rcut, [&](const auto &, const Vec &x, double r) {
        if (r < 1.e-8)
          return;
        ewald += std::erfc(alpha * r) / r;
        // Recurrence of the screened Coulomb radial derivatives. Keeping
        // erfc(beta*r) and the short-range damping correction separate avoids
        // cancellation against the bare inverse powers at long range.
        const double inv2 = 1. / (r * r), br = beta * r;
        double screened = std::erfc(br) / r;
        double gauss = 2 * beta / std::sqrt(pi) * std::exp(-br * br);
        double inverse = 1. / r, radial[4][5]{};
        double short_range[4];
        for (int l = 0; l < 4; ++l)
          short_range[l] = 0.5 * std::erfc((l == 0 ? 0.5 : 1.) * (r - radius));
        for (int l = 1; l <= 4; ++l) {
          inverse *= inv2;
          screened = (screened + gauss) * inv2;
          gauss *= 2 * beta * beta / (2 * l + 1);
          for (int k = 0; k < 4; ++k)
            radial[k][l] = damp[k] * (screened - inverse * short_range[k]);
        }
        t.real(x, radial);
      });
      for (const auto &kp : kpoints) {
        const double theta = dot(kp.k, delta);
        ewald += kp.charge * std::cos(theta);
        t.reciprocal(kp.k, kp.multi * std::cos(theta),
                     kp.multi * std::sin(theta));
      }
      if (i == j) {
        ewald -= 2 * alpha / std::sqrt(pi);
        const double s3 = 4 * std::pow(beta, 3) / (3 * std::sqrt(pi)) * damp[1];
        const double s5 =
            -8 * std::pow(beta, 5) / (15 * std::sqrt(pi)) * damp[3];
        for (int k = 0; k < 3; ++k)
          t.dd[4 * k] -= s3;
        for (int u = 0; u < 6; ++u) {
          const auto [a0, b0] = indices[u];
          if (a0 == b0)
            t.sq[u] += s3 / 3;
          for (int v = 0; v < 6; ++v) {
            const auto [c0, d0] = indices[v];
            t.qq[u + 6 * v] +=
                (a0 == b0 ? 1 : 2) * (c0 == d0 ? 1 : 2) * s5 *
                ((a0 == b0 && c0 == d0) + (a0 == c0 && b0 == d0) +
                 (a0 == d0 && b0 == c0));
          }
        }
      }
      for (int a0 = 0; a0 < 3; ++a0) {
        m.sd[3 * (j + size_t(n) * i) + a0] = dipscale[j] * t.sd[a0];
        for (int b0 = 0; b0 < 3; ++b0)
          m.dd[3 * j + a0 + 3 * size_t(n) * (3 * i + b0)] =
              dipscale[j] * dipscale[i] * t.dd[a0 + 3 * b0];
        for (int u = 0; u < 6; ++u)
          m.dq[3 * j + a0 + 3 * size_t(n) * (6 * i + u)] =
              dipscale[j] * t.dq[a0 + 3 * u];
      }
      for (int u = 0; u < 6; ++u) {
        m.sq[6 * (j + size_t(n) * i) + u] = t.sq[u];
        for (int v = 0; v < 6; ++v)
          m.qq[6 * j + u + 6 * size_t(n) * (6 * i + v)] = t.qq[u + 6 * v];
      }
      // The short-range correction to gamma is defined on the closest images.
      const auto near = lattice.nearest(i, j, delta);
      const double closest = std::min_element(near.begin(), near.end(),
                                              [](const Pair &x, const Pair &y) {
                                                return x.r < y.r;
                                              })
                                 ->r;
      Array weight;
      double sum = 0.;
      for (const auto &image : near) {
        const double x =
            std::clamp((image.r * image.r - closest * closest) / 0.3, 0., 1.);
        const double w = std::pow(1. - x, 3) * (1. + 3 * x + 6 * x * x);
        sum += w;
        weight.push_back(w);
      }
      // tau uses a fixed direct-space translation box around each Wigner-Seitz
      // representative (10 bohr). Keep the same image convention at the
      // boundary.
      int extent[3];
      for (int k = 0; k < 3; ++k)
        extent[k] =
            int(std::ceil(10 *
                          std::sqrt(dot(position(reciprocal.data(), k),
                                        position(reciprocal.data(), k))) /
                          (2 * pi)));
      Array distances;
      for (const auto &im : near)
        for (int x = -extent[0]; x <= extent[0]; ++x)
          for (int y = -extent[1]; y <= extent[1]; ++y)
            for (int w = -extent[2]; w <= extent[2]; ++w) {
              Vec d = im.d;
              for (int k = 0; k < 3; ++k)
                d[k] += x * a[k] + y * b[k] + w * c[k];
              const double r = std::sqrt(dot(d, d));
              if (r > 1.e-8)
                distances.push_back(r);
            }
      for (size_t u = offset[i]; u < offset[i + 1]; ++u)
        for (size_t v = offset[j]; v < offset[j + 1]; ++v) {
          const double h = 2 * hard[u] * hard[v] / (hard[u] + hard[v]);
          double g = ewald + (i == j ? h : 0.),
                 tau = i == j ? -0.2093327496 * h * h : 0.;
          for (size_t k = 0; k < near.size(); ++k) {
            const double r = near[k].r, d = std::exp(-0.2946211550 * r) / h;
            g -= weight[k] / sum * d / (r * (r + d));
          }
          for (double r : distances) {
            const double x = r / h;
            tau += 2.3 / near.size() * x * (1. - 0.65 * x) * std::exp(-1.3 * x);
          }
          m.gamma[u + ns * v] = g;
          m.tau[u + ns * v] = tau;
        }
    }
}
} // namespace gxtb
#endif
