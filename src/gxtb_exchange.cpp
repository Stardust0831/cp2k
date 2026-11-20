/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "gxtb_exchange.hpp"
#include "gxtb_parameters.hpp"

extern "C" void gxtb_product(int n, int transpose_a, int transpose_b,
                             double beta, const double *a, const double *b,
                             double *c);

namespace gxtb {
namespace {
constexpr double fraction = 0.15, omega = 0.2;
constexpr double shell_scale[] = {0.05, 2.5, 2.5, 2.5};
constexpr double charge_slope[] = {1.1, 0.55, 0.275, 0.1375};
} // namespace

Exchange::Exchange(const cp2k_gxtb_scf_problem &p)
    : ns(p.nsh), shell(p.nao), atom(p.nsh), kq(ns), kernels(size_t(ns) * ns) {
  std::vector<int> size(ns), local(ns);
  std::vector<std::vector<int>> on_atom(p.nat);
  for (int i = 0; i < p.nao; ++i) {
    shell[i] = p.ao_shell[i] - 1;
    ++size.at(shell[i]);
  }
  for (int u = 0; u < ns; ++u) {
    atom[u] = p.shell_atom[u] - 1;
    auto &list = on_atom.at(atom[u]);
    local[u] = int(list.size());
    const int l = (size[u] - 1) / 2;
    if (l < 0 || l > 3 || size[u] != 2 * l + 1 || local[u] >= 4)
      throw std::invalid_argument("Unsupported g-xTB shell");
    kq[u] = charge_slope[l];
    list.push_back(u);
  }
  Lattice lattice(p.cell, p.pbc);
  for (int i = 0; i < p.nat; ++i)
    for (int j = 0; j <= i; ++j) {
      const auto &a = element(p.z[i]), &b = element(p.z[j]);
      const auto images = lattice.nearest(
          i, j, difference(position(p.xyz, j), position(p.xyz, i)));
      const double radius = pair_radius(p.z[i], p.z[j]);
      const double kappa = -0.2960502355 + 0.0349817377 * radius;
      const double bond_scale = (a.bond_scale + b.bond_scale) * 0.5;
      const double bond_radius = std::sqrt(a.bond_radius * b.bond_radius);
      for (int u : on_atom[i])
        for (int v : on_atom[j]) {
          if (i == j && v > u)
            continue;
          const int l = (size[u] - 1) / 2, m = (size[v] - 1) / 2;
          const double h1 = a.hardness[local[u]], h2 = b.hardness[local[v]];
          const double geometric = std::sqrt(h1 * h2);
          const double power =
              std::max(a.average_power[local[u]], b.average_power[local[v]]);
          const double hardness =
              geometric * std::pow(2 * geometric / (h1 + h2), power - 1);
          const double scale = std::sqrt(shell_scale[l] * shell_scale[m]);
          const size_t uv = u + size_t(ns) * v;
          auto &k = kernels[uv];
          if (i == j) {
            k.mulliken = fraction * hardness * (u == v ? 1.3826597204 : scale);
            const int lo = std::min(l, m), hi = std::max(l, m);
            k.onsite = fraction * a.onsite[lo * (7 - lo) / 2 + hi];
            if (u == v && l > 0)
              k.rotation = 0.5 / size[u];
          }
          for (const auto &image : images) {
            const double r = image.r, w = 1. / images.size();
            // Evaluate the inverse denominator without forming an overflowing
            // exponential in a large vacuum cell.
            const double logd = -kappa * r - std::log(hardness * scale);
            double inv = 0, damp_fraction = 1;
            if (logd < std::log(std::numeric_limits<double>::max()) - 1) {
              const double d = std::exp(logd);
              inv = 1. / (r + d);
              damp_fraction = d * inv;
            }
            const double numerator =
                fraction + (1 - fraction) * std::erf(omega * r);
            const double dn = (1 - fraction) * 2 * omega / std::sqrt(pi) *
                              std::exp(-omega * omega * r * r);
            const double g = numerator * inv;
            const double dg = dn * inv - g * (inv - kappa * damp_fraction);
            const double arg = (r - bond_radius) / radius;
            const double bond = 0.5 * bond_scale * std::erfc(arg);
            const double db =
                -bond_scale * std::exp(-arg * arg) / (radius * std::sqrt(pi));
            k.mulliken += w * g;
            k.bond += w * bond;
            radial.push_back({image, uv, w * dg, w * db});
          }
          kernels[v + size_t(ns) * u] = k;
        }
    }
}

void Exchange::evaluate(cp2k_gxtb_scf_problem &p, double *energy,
                        bool derivatives) const {
  using Matrix = std::vector<double>;
  const size_t n = p.nao, nn = n * n;
  Matrix dG(derivatives ? kernels.size() : 0), dK(dG.size());
  for (int spin = 0; spin < p.nspin; ++spin) {
    const double *density = p.p + spin * nn;
    const double factor = p.nspin == 1 ? 1. : 2.;
    Matrix ps(nn), sps(nn), dp(nn), dx(nn), dy(nn);
    gxtb_product(n, 0, 0, 0., density, p.s, ps.data());
    gxtb_product(n, 0, 0, 0., p.s, ps.data(), sps.data());

    // X=P*S, Y=S*P*S. At fixed charges the energy is quadratic in P.
    // Accumulate partial derivatives wrt P, X and Y before differentiating
    // the two matrix products. Both AO matrices vary symmetrically.
    for (size_t j = 0; j < n; ++j)
      for (size_t i = 0; i < n; ++i) {
        const size_t ij = i + n * j, ji = j + n * i, ii = i + n * i,
                     jj = j + n * j;
        const int u = shell[i], v = shell[j];
        const size_t uv = u + size_t(ns) * v;
        const auto &k = kernels[uv];
        const double lambda =
            k.onsite * (1. - 0.5 * (kq[u] * p.qsh[u] + kq[v] * p.qsh[v]));
        const double t = ps[ij] * ps[ji], h = density[ij] * sps[ij];
        const double diag = density[ii] * sps[jj] + ps[ii] * ps[jj];
        const double d_m = -factor * (t + h) / 8, d_b = factor * t / 2;
        const double d_lambda = factor * (k.rotation * (t + h) / 4 -
                                          (h + ps[ij] * ps[ij] + diag) / 16);
        if (energy)
          energy[atom[u]] +=
              k.mulliken * d_m + k.bond * d_b + lambda * d_lambda;
        if (derivatives) {
          dG[uv] += d_m;
          dK[uv] += d_b;
        } else {
          p.vsh[u] -= 0.5 * k.onsite * kq[u] * d_lambda;
          p.vsh[v] -= 0.5 * k.onsite * kq[v] * d_lambda;
        }
        const double common =
            factor * (-k.mulliken / 8 + lambda * k.rotation / 4);
        const double ct = common + factor * k.bond / 2,
                     co = -factor * lambda / 16, ch = common + co;
        dx[ij] += ct * ps[ji] + 2 * co * ps[ij];
        dx[ji] += ct * ps[ij];
        dp[ij] += ch * sps[ij];
        dy[ij] += ch * density[ij];
        dp[ii] += co * sps[jj];
        dy[jj] += co * density[ii];
        dx[ii] += co * ps[jj];
        dx[jj] += co * ps[ii];
      }
    gxtb_product(n, 0, 0, 1., p.s, dy.data(), dx.data());
    if (derivatives) {
      Matrix ds(nn);
      gxtb_product(n, 0, 1, 0., dy.data(), ps.data(), ds.data());
      gxtb_product(n, 0, 0, 1., density, dx.data(), ds.data());
      for (size_t j = 0; j < n; ++j)
        for (size_t i = 0; i < n; ++i)
          p.weighted_p[spin * nn + i + n * j] -=
              0.5 * (ds[i + n * j] + ds[j + n * i]);
    } else {
      gxtb_product(n, 0, 0, 1., dx.data(), p.s, dp.data());
      for (size_t j = 0; j < n; ++j)
        for (size_t i = 0; i < n; ++i)
          p.exchange[spin * nn + i + n * j] +=
              0.5 * (dp[i + n * j] + dp[j + n * i]);
    }
  }
  if (derivatives) {
    Gradient result(p.nat);
    for (const auto &entry : radial) {
      const size_t uv = entry.shells, vu = uv / ns + (uv % ns) * ns;
      double slope = dG[uv] * entry.dmulliken + dK[uv] * entry.dbond;
      if (uv != vu)
        slope += dG[vu] * entry.dmulliken + dK[vu] * entry.dbond;
      result.add(entry.image, slope);
    }
    for (size_t i = 0; i < result.atoms.size(); ++i)
      for (int k = 0; k < 3; ++k)
        p.geometry_gradient[3 * i + k] += result.atoms[i][k];
    for (int k = 0; k < 9; ++k)
      p.geometry_strain[k] += result.strain[k];
  }
}
} // namespace gxtb
