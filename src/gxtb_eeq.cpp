/* SPDX-License-Identifier: GPL-2.0-or-later */
#if defined(__TBLITE)
#include "gxtb_integrals.hpp"
#include "gxtb_onecenter.hpp"
#include "gxtb_parameters.hpp"

extern "C" int gxtb_linear_solve(int, double *, double *);
namespace gxtb {
Array auxiliary_charges(int n, const int *z, const double *xyz,
                        const double *cell) {
  constexpr int periodic[] = {1, 1, 1};
  Array cn(n), local(n), width(n), hardness(n), chi(n), cap(n);
  for_each_pair(n, xyz, cell, periodic, 25., [&](const Pair &p) {
    const auto &a = basis_data(z[p.i]).eeq, &b = basis_data(z[p.j]).eeq;
    const double radius = a.cov + b.cov;
    const double count =
        0.5 * std::erfc(2 * (p.r - radius) / std::pow(radius, 0.75));
    cn[p.i] += count;
    cn[p.j] += count;
    local[p.i] += count * (b.en - a.en);
    local[p.j] -= count * (b.en - a.en);
  });
  for (int i = 0; i < n; ++i) {
    const auto &a = basis_data(z[i]).eeq;
    width[i] = a.width * (1. - 0.30 * cn[i] / std::pow(a.avg_cn, 0.75));
    hardness[i] = a.eta + std::sqrt(2 / pi) / width[i];
    chi[i] = -a.chi + a.kcn * cn[i] + a.kq * local[i];
  }
  const Lattice lattice(cell, periodic);
  const size_t dim = size_t(n) + 1;
  Array matrix(dim * dim), rhs(dim);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j <= i; ++j) {
      const auto &a = basis_data(z[i]).eeq, &b = basis_data(z[j]).eeq;
      const double radius =
          d3_radius_angstrom(z[i], z[j]) * 0.5 * (a.vdw + b.vdw);
      const double prefactor = 0.5 * std::sqrt(a.cap * b.cap);
      const double gamma = 1. / std::hypot(width[i], width[j]);
      double capacitance = 0, potential = 0;
      // Sum the periodic bond capacitance once per image. At this radius
      // erfc(8) is already below 1.2e-29. Self images contribute to the
      // diagonal hardness and cancel from the electronegativity difference.
      lattice.images(
          difference(position(xyz, i), position(xyz, j)),
          radius * (1. + 8. / 0.75), [&](const auto &, const Vec &, double r) {
            if (r < 1.e-8)
              return;
            const double c = prefactor * std::erfc(0.75 * (r / radius - 1.));
            capacitance += c;
            potential += c * std::erf(r * gamma) / r;
          });
      cap[i] += capacitance;
      matrix[i + dim * j] -= potential;
      if (i != j) {
        cap[j] += capacitance;
        matrix[j + dim * i] -= potential;
        const double delta = capacitance * (chi[i] - chi[j]);
        rhs[i] += delta;
        rhs[j] -= delta;
      }
    }
  for (int i = 0; i < n; ++i) {
    matrix[i * (dim + 1)] += 1. + hardness[i] * cap[i];
    matrix[i + dim * n] = matrix[n + dim * i] = 1.;
  }
  if (gxtb_linear_solve(dim, matrix.data(), rhs.data()))
    throw std::runtime_error("Singular g-xTB EEQBC charge system");
  rhs.resize(n);
  return rhs;
}
} // namespace gxtb
#endif
