/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "gxtb_scf.h"
#include "gxtb_electrostatic.hpp"
#include "gxtb_exchange.hpp"
#include "gxtb_repulsion.hpp"
#include "gxtb_revd4.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <vector>

extern "C" int gxtb_cholesky(int, double *);
extern "C" int gxtb_linear_solve(int, double *, double *);
extern "C" int gxtb_orbitals(int, const double *, double *, double *);
extern "C" void gxtb_density(int, const double *, const double *, double *);

namespace {
using Problem = cp2k_gxtb_scf_problem;
using Vector = std::vector<double>;

// Constrained least-squares extrapolation of successive AO Fock residuals.
// The normalization keeps the small DIIS system well scaled near convergence.
class Mixer {
  std::vector<Vector> inputs, residuals;

public:
  void next(Vector &f, const Vector &raw) {
    if (inputs.size() == 7) {
      inputs.erase(inputs.begin());
      residuals.erase(residuals.begin());
    }
    inputs.push_back(f);
    residuals.emplace_back(f.size());
    for (size_t k = 0; k < f.size(); ++k)
      residuals.back()[k] = raw[k] - f[k];
    const int m = static_cast<int>(inputs.size());
    const int dim = m + 1;
    Vector system(dim * dim), weights(dim);
    double scale = 0;
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j <= i; ++j) {
        const double overlap = std::inner_product(
            residuals[i].begin(), residuals[i].end(), residuals[j].begin(), 0.);
        system[i + dim * j] = system[j + dim * i] = overlap;
      }
      scale = std::max(scale, system[i + dim * i]);
    }
    if (scale == 0) {
      f = raw;
      return;
    }
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j)
        system[i + dim * j] /= scale;
      system[i + dim * i] += 1.e-13;
      system[i + dim * m] = system[m + dim * i] = 1.;
    }
    weights[m] = 1.;
    if (gxtb_linear_solve(dim, system.data(), weights.data()) != 0) {
      for (size_t k = 0; k < f.size(); ++k)
        f[k] += 0.2 * residuals.back()[k];
      return;
    }
    std::fill(f.begin(), f.end(), 0.);
    for (int i = 0; i < m; ++i)
      for (size_t k = 0; k < f.size(); ++k)
        f[k] += weights[i] * (inputs[i][k] + 0.2 * residuals[i][k]);
  }
};

void populations(Problem &a) {
  const size_t n = a.nao, ns = a.nsh, na = a.nat, nn = n * n;
  std::fill_n(a.qsh, ns * a.nspin, 0.);
  std::copy_n(a.neutral_shell, ns, a.qsh);
  std::fill_n(a.qat, na * a.nspin, 0.);
  std::fill_n(a.dip, 3 * na * a.nspin, 0.);
  std::fill_n(a.quad, 6 * na * a.nspin, 0.);
  std::fill_n(a.energy, na, 0.);
  for (int spin = 0; spin < a.nspin; ++spin) {
    for (size_t i = 0; i < n; ++i) {
      const int shell = a.ao_shell[i] - 1, atom = a.shell_atom[shell] - 1;
      for (size_t j = 0; j < n; ++j) {
        const size_t ij = i * n + j;
        const double p = a.p[ij + spin * nn];
        a.energy[atom] += p * a.h0[ij];
        for (int channel = 0; channel < a.nspin; ++channel) {
          const double value = spin && channel ? p : -p;
          const size_t site = atom + channel * na;
          a.qsh[shell + channel * ns] += value * a.s[ij];
          for (int k = 0; k < 3; ++k)
            a.dip[3 * site + k] += value * a.dip_int[3 * ij + k];
          for (int k = 0; k < 6; ++k)
            a.quad[6 * site + k] += value * a.quad_int[6 * ij + k];
        }
      }
    }
  }
  for (int channel = 0; channel < a.nspin; ++channel)
    for (size_t sh = 0; sh < ns; ++sh)
      a.qat[channel * na + a.shell_atom[sh] - 1] += a.qsh[channel * ns + sh];
}

void assemble(Problem &a, Vector &f) {
  const size_t n = a.nao, nn = n * n;
  for (int channel = 0; channel < a.nspin; ++channel) {
    for (int sh = 0; sh < a.nsh; ++sh)
      a.vsh[sh + channel * a.nsh] +=
          a.vat[a.shell_atom[sh] - 1 + channel * a.nat];
    for (size_t i = 0; i < n; ++i)
      a.vao[i + channel * n] += a.vsh[a.ao_shell[i] - 1 + channel * a.nsh];
  }
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j <= i; ++j) {
      const size_t ij = i * n + j, ji = j * n + i;
      double v[2] = {a.h0[ij], 0.};
      for (int c = 0; c < a.nspin; ++c) {
        const size_t x = a.shell_atom[a.ao_shell[i] - 1] - 1 + c * a.nat;
        const size_t y = a.shell_atom[a.ao_shell[j] - 1] - 1 + c * a.nat;
        double shift = (a.vao[c * n + i] + a.vao[c * n + j]) * a.s[ij];
        for (int k = 0; k < 3; ++k)
          shift += a.vdip[3 * x + k] * a.dip_int[3 * ij + k] +
                   a.vdip[3 * y + k] * a.dip_int[3 * ji + k];
        for (int k = 0; k < 6; ++k)
          shift += a.vquad[6 * x + k] * a.quad_int[6 * ij + k] +
                   a.vquad[6 * y + k] * a.quad_int[6 * ji + k];
        v[c] -= 0.5 * shift;
      }
      for (int spin = 0; spin < a.nspin; ++spin)
        f[spin * nn + ij] = f[spin * nn + ji] =
            v[0] + (spin ? -v[1] : v[1]) + a.exchange[spin * nn + ij];
    }
}

void density(Problem &a, bool weighted) {
  const size_t n = a.nao, nn = n * n;
  Vector w(n);
  for (int spin = 0; spin < a.nspin; ++spin) {
    for (size_t j = 0; j < n; ++j) {
      w[j] = a.nspin == 1 ? a.occ[j] + a.occ[n + j] : a.occ[spin * n + j];
      if (weighted)
        w[j] *= a.eps[spin * n + j];
    }
    gxtb_density(a.nao, a.c + spin * nn, w.data(),
                 (weighted ? a.weighted_p : a.p) + spin * nn);
  }
}

class Interactions {
  const cp2k_gxtb_d4_v1 *d4 = cp2k_gxtb_native_d4();
  std::unique_ptr<void, void (*)(void *)> dispersion{nullptr, d4->destroy};
  gxtb::Repulsion repulsion;
  gxtb::Exchange exchange;
  gxtb::Electrostatic electrostatic;
  Vector ed4, vd4;

public:
  explicit Interactions(const Problem &a)
      : repulsion(a.nat, a.z, a.xyz, a.cell, a.pbc), exchange(a),
        electrostatic(a), ed4(a.nat), vd4(a.nat) {
    char message[256]{};
    dispersion.reset(
        d4->create(a.nat, a.z, a.xyz, a.cell, a.pbc, message, sizeof(message)));
    if (!dispersion)
      throw std::runtime_error(message);
  }
  void evaluate(Problem &a, double *energy) {
    const size_t n = a.nao, sites = size_t(a.nat) * a.nspin;
    std::fill_n(a.vat, sites, 0.);
    std::fill_n(a.vdip, 3 * sites, 0.);
    std::fill_n(a.vquad, 6 * sites, 0.);
    std::fill_n(a.vsh, size_t(a.nsh) * a.nspin, 0.);
    std::fill_n(a.vao, n * a.nspin, 0.);
    std::fill_n(a.exchange, n * n * a.nspin, 0.);
    electrostatic.evaluate(a, energy);
    repulsion.add_scf(a.qat, energy, a.vat);
    exchange.evaluate(a, energy);
    d4->scf(dispersion.get(), a.qat, ed4.data(), vd4.data());
    for (int i = 0; i < a.nat; ++i) {
      if (energy)
        energy[i] += ed4[i];
      a.vat[i] += vd4[i];
    }
  }
  void finish(Problem &a) {
    repulsion.add_core(a.energy);
    repulsion.gradient(a.qat, a.geometry_gradient, a.geometry_strain);
    exchange.evaluate(a, nullptr, true);
    Vector grad(3 * size_t(a.nat)), atm_grad(grad.size());
    gxtb::Tensor strain{}, atm_strain{};
    d4->derivatives(dispersion.get(), a.qat, grad.data(), strain.data());
    d4->atm(dispersion.get(), ed4.data(), atm_grad.data(), atm_strain.data());
    for (int i = 0; i < a.nat; ++i)
      a.energy[i] += ed4[i];
    for (size_t k = 0; k < grad.size(); ++k)
      a.geometry_gradient[k] += grad[k] + atm_grad[k];
    for (size_t k = 0; k < strain.size(); ++k)
      a.geometry_strain[k] += strain[k] + atm_strain[k];
  }
};

int solve(Problem &a, char *message, size_t length) {
  if (a.nao < 1 || a.nat < 1 || a.nsh < 1 || (a.nspin != 1 && a.nspin != 2))
    throw std::invalid_argument("Invalid g-xTB SCF dimensions");
  const size_t n = a.nao, nn = n * n;
  Vector chol(a.s, a.s + nn), raw(nn * a.nspin), f(raw.size());
  Interactions interactions(a);
  if (gxtb_cholesky(a.nao, chol.data()) != 0)
    throw std::runtime_error("Non-positive-definite g-xTB overlap");
  for (int spin = 0; spin < 2; ++spin) {
    if (!std::isfinite(a.nel[spin]) || a.nel[spin] < -1.e-8 ||
        a.nel[spin] > n + 1.e-8)
      throw std::invalid_argument("Invalid g-xTB electron count");
    for (size_t j = 0; j < n; ++j)
      a.occ[spin * n + j] = std::clamp(a.nel[spin] - double(j), 0., 1.);
  }
  Mixer mixer;
  // The first operator uses the supplied EEQBC guess. Every subsequent
  // population and multipole is a contraction of the same AO density.
  interactions.evaluate(a, nullptr);
  assemble(a, f);
  double previous = std::numeric_limits<double>::infinity(), de = previous,
         residual = previous;
  for (int iteration = 1; iteration <= a.max_iter; ++iteration) {
    if (!std::all_of(f.begin(), f.end(),
                     [](double x) { return std::isfinite(x); }))
      throw std::runtime_error("Nonfinite g-xTB Fock matrix");
    std::copy(f.begin(), f.end(), a.c);
    for (int spin = 0; spin < a.nspin; ++spin)
      if (gxtb_orbitals(a.nao, chol.data(), a.c + spin * nn,
                        a.eps + spin * n) != 0)
        throw std::runtime_error("g-xTB orbital diagonalization failed");
    density(a, false);
    populations(a);
    interactions.evaluate(a, a.energy);
    const double energy = std::accumulate(a.energy, a.energy + a.nat, 0.);
    assemble(a, raw);
    if (!std::isfinite(energy) ||
        !std::all_of(raw.begin(), raw.end(),
                     [](double x) { return std::isfinite(x); }))
      throw std::runtime_error(
          "Nonfinite g-xTB electronic energy or potential");
    residual = 0.;
    for (size_t k = 0; k < f.size(); ++k)
      residual = std::max(residual, std::abs(raw[k] - f[k]));
    de = std::abs(energy - previous);
    if (de < a.energy_tol && residual < a.residual_tol) {
      density(a, true);
      interactions.finish(a);
      return iteration;
    }
    previous = energy;
    mixer.next(f, raw);
  }
  std::snprintf(message, length,
                "CP2K g-xTB SCF exceeded %d cycles (dE %.3g Ha, max Fock "
                "residual %.3g Ha)",
                a.max_iter, de, residual);
  return 0;
}
} // namespace

extern "C" int cp2k_gxtb_scf(Problem *p, size_t size, char *msg, size_t len) {
  try {
    if (!p || size != sizeof(*p))
      throw std::invalid_argument(
          "Incompatible g-xTB electronic problem descriptor");
    return solve(*p, msg, len);
  } catch (const std::exception &e) {
    std::snprintf(msg, len, "%s", e.what());
    return 0;
  }
}
