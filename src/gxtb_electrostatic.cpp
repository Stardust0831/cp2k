/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "gxtb_electrostatic.hpp"
#include "gxtb_pairs.hpp"
#include "gxtb_parameters.hpp"

extern "C" void gxtb_matvec(int rows, int cols, int transpose, double scale,
                            const double *matrix, const double *input,
                            double *output);
namespace gxtb {
Electrostatic::Electrostatic(const cp2k_gxtb_scf_problem &p)
    : ipea(p.nsh), third(p.nsh), fourth(p.nsh), reference_field(p.nsh),
      spin_matrix(p.nspin == 2 ? size_t(p.nsh) * p.nsh : 0) {
  std::vector<int> size(p.nsh), ordinal(p.nat), angular(p.nsh);
  std::vector<double> reference(p.nsh);
  for (int i = 0; i < p.nao; ++i)
    ++size.at(p.ao_shell[i] - 1);
  for (int u = 0; u < p.nsh; ++u) {
    const int i = p.shell_atom[u] - 1, local = ordinal.at(i)++;
    const int l = (size[u] - 1) / 2;
    if (local >= 4 || l < 0 || l > 3)
      throw std::invalid_argument("Invalid g-xTB electrostatic shell");
    angular[u] = l;
    const Element &e = element(p.z[i]);
    const double root = std::sqrt(p.cn[i] + 1.e-12) - 1.e-6;
    ipea[u] = e.ipea[local] * (1 + e.ipea_cn * root);
    third[u] = e.hubbard_derivative * (1 + 0.15 * l);
    fourth[u] = 0.036 * (1 + 0.15 * l);
    reference[u] = e.reference_charge[local];
  }
  for (int v = 0; v < p.nsh; ++v)
    for (int u = 0; u < p.nsh; ++u) {
      const size_t uv = u + size_t(p.nsh) * v;
      if (p.shell_atom[u] != p.shell_atom[v]) {
        reference_field[u] += p.gamma[uv] * reference[v];
      } else if (p.nspin == 2) {
        const auto [lo, hi] = std::minmax(angular[u], angular[v]);
        spin_matrix[uv] =
            element(p.z[p.shell_atom[u] - 1]).spin[lo * (7 - lo) / 2 + hi];
      }
    }
}

void Electrostatic::evaluate(cp2k_gxtb_scf_problem &p, double *energy) const {
  const int n = p.nat, ns = p.nsh;
  // Symmetric quadratic form of the atomic charges, dipoles and traceless
  // quadrupoles. Off-diagonal blocks contribute to both moment potentials.
  gxtb_matvec(3 * n, n, 0, 1., p.sd, p.qat, p.vdip);
  gxtb_matvec(3 * n, n, 1, 1., p.sd, p.dip, p.vat);
  gxtb_matvec(3 * n, 3 * n, 0, 1., p.dd, p.dip, p.vdip);
  gxtb_matvec(6 * n, n, 0, 1., p.sq, p.qat, p.vquad);
  gxtb_matvec(6 * n, n, 1, 1., p.sq, p.quad, p.vat);
  gxtb_matvec(3 * n, 6 * n, 0, 1., p.dq, p.quad, p.vdip);
  gxtb_matvec(3 * n, 6 * n, 1, 1., p.dq, p.dip, p.vquad);
  gxtb_matvec(6 * n, 6 * n, 0, 1. / 3, p.qq, p.quad, p.vquad);
  std::vector<double> split(n), ds(n);
  for (int i = 0; i < n; ++i) {
    if (energy) {
      double e = p.qat[i] * p.vat[i];
      for (int k = 3 * i; k < 3 * i + 3; ++k)
        e += p.dip[k] * p.vdip[k];
      for (int k = 6 * i; k < 6 * i + 6; ++k)
        e += p.quad[k] * p.vquad[k];
      energy[i] += 0.5 * e;
    }
    const double low = p.qat[i] - 0.6666666660, high = p.qat[i] + 0.6666666660;
    split[i] = 1. + 0.0125 * (std::erf(low) + std::erf(high));
    ds[i] =
        0.025 / std::sqrt(pi) * (std::exp(-low * low) + std::exp(-high * high));
  }
  std::vector<double> gq(ns), tq(ns), q2(ns);
  gxtb_matvec(ns, ns, 0, 1., p.gamma, p.qsh, gq.data());
  gxtb_matvec(ns, ns, 0, 1., p.tau, p.qsh, tq.data());
  if (p.nspin == 2)
    gxtb_matvec(ns, ns, 0, 1., spin_matrix.data(), p.qsh + ns, p.vsh + ns);
  for (int u = 0; u < ns; ++u) {
    const int i = p.shell_atom[u] - 1;
    const double q = p.qsh[u], qsq = q * q,
                 linear = ipea[u] * split[i] + reference_field[u];
    q2[u] = third[u] * qsq;
    p.vat[i] += ipea[u] * q * ds[i];
    p.vsh[u] += linear + gq[u] + 2. / 3 * third[u] * q * tq[u] +
                fourth[u] * qsq * q / 6;
    if (energy) {
      energy[i] += q * (linear + 0.5 * gq[u]) + q2[u] * tq[u] / 3 +
                   fourth[u] * qsq * qsq / 24;
      if (p.nspin == 2)
        energy[i] += 0.5 * p.qsh[u + ns] * p.vsh[u + ns];
    }
  }
  gxtb_matvec(ns, ns, 0, 1. / 3, p.tau, q2.data(), p.vsh);
}
} // namespace gxtb
