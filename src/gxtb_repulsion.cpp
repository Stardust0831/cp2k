/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "gxtb_repulsion.hpp"
#include "gxtb_parameters.hpp"
#include <map>

namespace gxtb {
namespace {
constexpr double cutoff = 25.0, cn_steepness = 2.068, cn_epsilon = 1.e-6;
constexpr double short_weight = 0.0046511298, short_exponent = 0.73;
constexpr double radius_unit = 3.75372519958925, k3 = 0.5868618150,
                 k4 = 2.3046216529;
} // namespace

Repulsion::Repulsion(int n, const int *z, const double *xyz, const double *cell,
                     const int *pbc)
    : species(z, z + n) {
  std::vector<double> cn(n), alpha(n), dalpha(n);
  for_each_pair(n, xyz, cell, pbc, cutoff, [&](const Pair &p) {
    const double radius =
        element(species[p.i]).radius + element(species[p.j]).radius;
    const double slope = cn_steepness / radius;
    const double u = slope * (p.r - radius);
    const double count = 0.5 * std::erfc(u);
    cn[p.i] += count;
    cn[p.j] += count;
    edges.push_back(
        {p, -slope * std::exp(-u * u) / std::sqrt(pi), 0., 0., 0., 0.});
  });
  for (int i = 0; i < n; ++i) {
    const Element &el = element(species[i]);
    const double root = std::sqrt(cn[i] + cn_epsilon * cn_epsilon);
    alpha[i] = el.alpha * (1.0 + el.kcn * (root - cn_epsilon));
    dalpha[i] = 0.5 * el.alpha * el.kcn / root;
  }
  // Distinct images of the same basis pair have a common SCF charge factor.
  std::map<std::pair<int, int>, double> summed;
  for (auto &edge : edges) {
    const int i = edge.i, j = edge.j;
    const Element &a = element(species[i]), &b = element(species[j]);
    const double radius = pair_radius(species[i], species[j]);
    const double inv = 1. / edge.r, x = radius * inv / radius_unit;
    const double k2 = 0.5 * ((species[i] <= 2 ? 0.1704681501 : 0.1203923899) +
                             (species[j] <= 2 ? 0.1704681501 : 0.1203923899));
    const double t1 = 0.5 * (a.k1 + b.k1) * inv, t2 = k2 * x * x,
                 t3 = k3 * x * x * x, t4 = k4 * x * x * x * x;
    const double poly = 1.0 + t1 + t2 + t3 + t4;
    const double dpoly = -inv * (t1 + 2 * t2 + 3 * t3 + 4 * t4);
    const double sum = alpha[i] + alpha[j], width = alpha[i] * alpha[j] / sum;
    const double r = edge.r + std::sqrt(a.offset * b.offset),
                 root = std::sqrt(r);
    const double u = r * root, v = short_exponent * r * r;
    const double long_term = std::exp(-width * u),
                 short_term = short_weight * std::exp(-width * v);
    const double damp = long_term + short_term;
    const double dk_dwidth = -poly * (u * long_term + v * short_term);
    edge.kernel = poly * damp;
    edge.dk_dr = dpoly * damp - poly * width *
                                    (1.5 * root * long_term +
                                     2 * short_exponent * r * short_term);
    edge.dk_dcn_i = dk_dwidth * std::pow(alpha[j] / sum, 2) * dalpha[i];
    edge.dk_dcn_j = dk_dwidth * std::pow(alpha[i] / sum, 2) * dalpha[j];
    summed[std::minmax(i, j)] += edge.kernel;
  }
  for (const auto &entry : summed)
    coefficients.push_back(
        {entry.first.first, entry.first.second, entry.second});
}

double Repulsion::effective_charge(size_t i, const double *charge) const {
  const Element &e = element(species[i]);
  return e.zeff * (1. - e.kq * charge[i]);
}

void Repulsion::add_core(double *energy) const {
  for (size_t i = 0; i < species.size(); ++i)
    energy[i] += element(species[i]).core;
}

void Repulsion::add_scf(const double *charge, double *energy,
                        double *potential) const {
  for (const auto &pair : coefficients) {
    const int i = pair.i, j = pair.j;
    const Element &a = element(species[i]), &b = element(species[j]);
    const double zi = effective_charge(i, charge),
                 zj = effective_charge(j, charge);
    if (energy) {
      const double half = 0.5 * zi * zj * pair.kernel;
      energy[i] += half;
      energy[j] += half;
    }
    potential[i] -= a.zeff * a.kq * zj * pair.kernel;
    potential[j] -= b.zeff * b.kq * zi * pair.kernel;
  }
}

void Repulsion::gradient(const double *charge, double *grad,
                         double *strain) const {
  Gradient result(species.size());
  std::vector<double> dE_dcn(species.size());
  for (const auto &edge : edges) {
    const double zz =
        effective_charge(edge.i, charge) * effective_charge(edge.j, charge);
    result.add(edge, zz * edge.dk_dr);
    dE_dcn[edge.i] += zz * edge.dk_dcn_i;
    dE_dcn[edge.j] += zz * edge.dk_dcn_j;
  }
  for (const auto &edge : edges)
    result.add(edge, edge.dcn_dr * (dE_dcn[edge.i] + dE_dcn[edge.j]));
  result.output(grad, strain);
}
} // namespace gxtb
