/* SPDX-License-Identifier: GPL-2.0-or-later */
#if defined(__TBLITE)
#include "gxtb_electrostatic.hpp"
#include "gxtb_integrals.hpp"
#include "gxtb_scf.h"
#include <cstdio>
#include <numeric>

namespace gxtb {
namespace {
using Problem = cp2k_gxtb_scf_problem;
void operators(Problem &p, const Matrices &m) {
  p.s = m.s.data();
  p.h0 = m.h.data();
  p.dip_int = m.dip.data();
  p.quad_int = m.quad.data();
  p.cn = m.cn.data();
  p.gamma = m.gamma.data();
  p.tau = m.tau.data();
  p.sd = m.sd.data();
  p.dd = m.dd.data();
  p.sq = m.sq.data();
  p.dq = m.dq.data();
  p.qq = m.qq.data();
}
struct Potentials {
  Array at, sh, ao, dp, qp;
  explicit Potentials(Problem &p)
      : at(size_t(p.nat) * p.nspin), sh(size_t(p.nsh) * p.nspin),
        ao(size_t(p.nao) * p.nspin), dp(at.size() * 3), qp(at.size() * 6) {
    p.vat = at.data();
    p.vsh = sh.data();
    p.vao = ao.data();
    p.vdip = dp.data();
    p.vquad = qp.data();
  }
};

double electrostatic_energy(Problem p, const Matrices &m) {
  operators(p, m);
  Potentials field(p);
  Array energy(p.nat);
  Electrostatic(p).evaluate(p, energy.data());
  return std::accumulate(energy.begin(), energy.end(), 0.);
}

// Explicit geometry response of the converged electronic Lagrangian. The
// orbitals stay fixed; the overlap response includes the Pulay and exchange
// terms. Repulsion, dispersion and the radial exchange response are analytic.
double response(const Problem &p, const Matrices &plus, const Matrices &minus,
                double step) {
  double derivative =
      electrostatic_energy(p, plus) - electrostatic_energy(p, minus);
  const size_t n = p.nao;
  for (int spin = 0; spin < p.nspin; ++spin) {
    const double sign = spin == 0 ? 1. : -1.;
    for (size_t i = 0; i < n; ++i) {
      const int atom = p.shell_atom[p.ao_shell[i] - 1] - 1;
      for (size_t j = 0; j < n; ++j) {
        const size_t ij = i * n + j, k = spin * n * n + ij;
        const double ds = plus.s[ij] - minus.s[ij];
        double h = plus.h[ij] - minus.h[ij];
        for (int channel = 0; channel < p.nspin; ++channel) {
          const size_t site = atom + channel * p.nat;
          double shift = p.vao[i + channel * n] * ds;
          for (int u = 0; u < 3; ++u)
            shift += p.vdip[3 * site + u] *
                     (plus.dip[3 * ij + u] - minus.dip[3 * ij + u]);
          for (int u = 0; u < 6; ++u)
            shift += p.vquad[6 * site + u] *
                     (plus.quad[6 * ij + u] - minus.quad[6 * ij + u]);
          h -= (channel ? sign : 1.) * shift;
        }
        derivative += p.p[k] * h - p.weighted_p[k] * ds;
      }
    }
  }
  return derivative / (2 * step);
}
} // namespace
} // namespace gxtb

extern "C" int cp2k_gxtb_internal(int n, const int *z, const double *xyz,
                                  const double *cell, const double *radii,
                                  int forces, int virial, double *total,
                                  double *grad, double *strain, char *msg,
                                  size_t length) {
  using namespace gxtb;
  try {
    if (radii[0] != 0 || radii[1] != 0)
      throw std::invalid_argument(
          "Omit EXCHANGE_RADII for the CP2K g-xTB solver");
    Orbitals orb(n, z);
    Matrices m(n, z, xyz, cell, orb);
    Problem p{};
    p.nat = n;
    p.nsh = orb.neutral.size();
    p.nao = orb.ao_shell.size();
    const int electrons = std::lround(
        std::accumulate(orb.neutral.begin(), orb.neutral.end(), 0.));
    p.nspin = 1 + (electrons % 2);
    p.nel[0] = (electrons + 1) / 2;
    p.nel[1] = electrons / 2;
    p.max_iter = 250;
    p.energy_tol = 1.e-11;
    p.residual_tol = 1.e-9;
    constexpr int periodic[] = {1, 1, 1};
    p.z = z;
    p.pbc = periodic;
    p.xyz = xyz;
    p.cell = cell;
    p.ao_shell = orb.ao_shell.data();
    p.shell_atom = orb.shell_atom.data();
    p.neutral_shell = orb.neutral.data();
    p.geometry_gradient = grad;
    p.geometry_strain = strain;
    operators(p, m);
    const size_t nao = p.nao, nn = nao * nao * p.nspin,
                 sites = size_t(n) * p.nspin;
    Array density(nn), coeff(nn), exchange(nn), weighted(nn),
        eps(nao * p.nspin), occ(2 * nao), qat(sites),
        qsh(size_t(p.nsh) * p.nspin), dip(3 * sites), quad(6 * sites),
        energy(n);
    p.p = density.data();
    p.c = coeff.data();
    p.eps = eps.data();
    p.occ = occ.data();
    p.exchange = exchange.data();
    p.weighted_p = weighted.data();
    p.qat = qat.data();
    p.qsh = qsh.data();
    p.dip = dip.data();
    p.quad = quad.data();
    p.energy = energy.data();
    Potentials field(p);
    std::copy(m.aux.begin(), m.aux.end(), qat.begin());
    Array neutral_atom(n);
    std::vector<int> shell_size(p.nsh);
    for (int sh : orb.ao_shell)
      ++shell_size[sh - 1];
    for (int sh = 0; sh < p.nsh; ++sh)
      neutral_atom[orb.shell_atom[sh] - 1] += orb.neutral[sh];
    for (int sh = 0; sh < p.nsh; ++sh) {
      const int a = orb.shell_atom[sh] - 1;
      qsh[sh] = qat[a] * orb.neutral[sh] / neutral_atom[a];
    }
    for (size_t i = 0; i < nao; ++i) {
      const int sh = orb.ao_shell[i] - 1;
      for (int spin = 0; spin < p.nspin; ++spin)
        density[spin * nao * nao + i * (nao + 1)] =
            (0.5 * orb.neutral[sh] - qsh[sh] / p.nspin) / shell_size[sh];
    }
    if (!cp2k_gxtb_scf(&p, sizeof(p), msg, length))
      return 1;
    *total = std::accumulate(energy.begin(), energy.end(), 0.);
    constexpr double dr = 1.e-5, de = 5.e-7;
    if (forces)
      for (int k = 0; k < 3 * n; ++k) {
        Array pos(xyz, xyz + 3 * n);
        pos[k] += dr;
        Matrices plus(n, z, pos.data(), cell, orb);
        pos[k] -= 2 * dr;
        Matrices minus(n, z, pos.data(), cell, orb);
        grad[k] += response(p, plus, minus, dr);
      }
    if (virial)
      for (int a = 0; a < 3; ++a)
        for (int b = 0; b <= a; ++b) {
          // Fixed fractional coordinates; the model is invariant under
          // rotation.
          auto stretched = [&](double e) {
            Array pos(xyz, xyz + 3 * n), box(cell, cell + 9);
            for (int i = 0; i < n; ++i)
              pos[3 * i + a] += e * xyz[3 * i + b];
            for (int i = 0; i < 3; ++i)
              box[3 * i + a] += e * cell[3 * i + b];
            return Matrices(n, z, pos.data(), box.data(), orb);
          };
          const double v = response(p, stretched(de), stretched(-de), de);
          strain[3 * b + a] += v;
          if (a != b)
            strain[3 * a + b] += v;
        }
    return 0;
  } catch (const std::exception &e) {
    std::snprintf(msg, length, "%s", e.what());
    return 1;
  }
}
#endif
