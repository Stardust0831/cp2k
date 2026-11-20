/* SPDX-License-Identifier: GPL-2.0-or-later */
#if defined(__TBLITE)
#include "gxtb_integrals.hpp"
#include "gxtb_onecenter.hpp"
#include "gxtb_parameters.hpp"

namespace gxtb {
namespace {
#include "gxtb_basis_data.inc"
struct Gaussian {
  int l, n;
  double exponent[12]{}, coefficient[12]{};
};
extern "C" void gxtb_moments(const Gaussian *, const Gaussian *, const double *,
                             double *, double *, double *);
extern "C" void gxtb_overlap(const Gaussian *, const Gaussian *, const double *,
                             double *);
extern "C" void gxtb_diatomic(int, int, const double *, const double *,
                              double *);
constexpr int pbc[3] = {1, 1, 1};

double primitive_norm(int l, double a) {
  constexpr double factorial[] = {1., 1., 3., 15.};
  return std::pow(2 * a / pi, 0.75) *
         std::sqrt(std::pow(4 * a, l) / factorial[l]);
}
Gaussian contracted(const ShellData &data, double scale, double charge) {
  Gaussian g{data.l, data.n};
  double norm = 0;
  for (int i = 0; i < g.n; ++i) {
    g.exponent[i] = scale * data.alpha[i];
    g.coefficient[i] = data.c0[i] + charge * data.c1[i];
    for (int j = 0; j <= i; ++j) {
      const double a = data.alpha[i], b = data.alpha[j];
      const double overlap =
          std::pow(2 * std::sqrt(a * b) / (a + b), g.l + 1.5);
      norm += (i == j ? 1 : 2) * g.coefficient[i] * g.coefficient[j] * overlap;
    }
  }
  for (int i = 0; i < g.n; ++i)
    g.coefficient[i] *= primitive_norm(g.l, g.exponent[i]) / std::sqrt(norm);
  return g;
}
struct AtomBasis {
  std::vector<Gaussian> scf, h0, proj;
  Array level;
  int first, projector_start;
  int size() const {
    const int shells = scf.size();
    return shells * shells;
  }
};
} // namespace

const BasisData &basis_data(int z) {
  if (z < 1 || z > 103)
    throw std::invalid_argument("g-xTB requires Z=1..103");
  return basis_table[z - 1];
}

Orbitals::Orbitals(int nat, const int *z) {
  for (int i = 0; i < nat; ++i) {
    start.push_back(ao_shell.size());
    const auto &b = basis_data(z[i]);
    for (int sh = 0; sh < b.count; ++sh) {
      shell_atom.push_back(i + 1);
      neutral.push_back(b.occ[sh]);
      ao_shell.insert(ao_shell.end(), 2 * shell_table[b.start + sh].l + 1,
                      shell_atom.size());
    }
  }
}

Matrices::Matrices(int nat, const int *z, const double *xyz, const double *cell,
                   const Orbitals &orb)
    : cn(nat), aux(auxiliary_charges(nat, z, xyz, cell)) {
  const size_t n = orb.ao_shell.size();
  s.resize(n * n);
  h.resize(n * n);
  dip.resize(3 * n * n);
  quad.resize(6 * n * n);
  Array basis_cn(nat);
  std::vector<Vec> anisotropy(nat);
  // Both coordination numbers use the corresponding error-function counting
  // rule and its 25-bohr sphere. The anisotropy is a vector sum over neighbors.
  for_each_pair(nat, xyz, cell, pbc, 25., [&](const Pair &e) {
    const auto &a = basis_data(z[e.i]), &b = basis_data(z[e.j]);
    const double cn_radius = element(z[e.i]).radius + element(z[e.j]).radius;
    const double c = 0.5 * std::erfc(2.068 * (e.r / cn_radius - 1.));
    const double bc =
        0.5 * std::erfc(3.75 * (e.r / (a.cn_radius + b.cn_radius) - 1.));
    cn[e.i] += c;
    cn[e.j] += c;
    basis_cn[e.i] += bc;
    basis_cn[e.j] += bc;
    const double v = 0.25 * (a.anisotropy + b.anisotropy) *
                     std::erfc(1.5 * (e.r - pair_radius(z[e.i], z[e.j])));
    for (int k = 0; k < 3; ++k) {
      anisotropy[e.i][k] += e.d[k] * v;
      anisotropy[e.j][k] -= e.d[k] * v;
    }
  });
  std::vector<AtomBasis> atoms(nat);
  int np = 0;
  double min_basis = 1.e30, min_proj = 1.e30;
  for (int i = 0; i < nat; ++i) {
    auto &at = atoms[i];
    const auto &b = basis_data(z[i]);
    at.first = orb.start[i];
    at.projector_start = np;
    const double q = aux[i], bc = basis_cn[i];
    const double charge = q * (1 - b.charge[1] * q);
    const double root = b.charge[0] * (std::sqrt(bc + 1.e-12) - 1.e-6);
    const double cross = b.charge[2] * bc * q;
    for (int u = 0; u < b.count; ++u) {
      const auto &sh = shell_table[b.start + u];
      at.scf.push_back(contracted(sh, 1., charge + root + cross));
      at.h0.push_back(contracted(sh, b.h_exp[u],
                                 b.scale[0] * charge + b.scale[1] * root +
                                     b.scale[2] * cross));
      for (int p = 0; p < sh.n; ++p)
        min_basis =
            std::min({min_basis, sh.alpha[p], sh.alpha[p] * b.h_exp[u]});
      at.level.push_back(b.level[u] - cn[i] * b.cn_level[u]);
    }
    for (int u = 0; u < b.projectors; ++u) {
      Gaussian g{b.proj_l[u], 1};
      g.exponent[0] = b.proj_exp[u];
      g.coefficient[0] = primitive_norm(g.l, g.exponent[0]);
      at.proj.push_back(g);
      np += 2 * g.l + 1;
      min_proj = std::min(min_proj, g.exponent[0]);
    }
  }
  Array projector(size_t(np) * n), proj_level(np);
  const double cutoff = std::min(
      40.,
      std::sqrt(25. * (1 / min_basis + 1 / std::min(min_proj, min_basis))));
  Lattice lattice(cell, pbc);
  constexpr double wolfsberg[] = {0.7129720394, 0.7981384369, 1.425, 2.85};
  for (int i = 0; i < nat; ++i) {
    const auto &a = atoms[i];
    const auto &ai = basis_data(z[i]);
    for (int j = 0; j < nat; ++j) {
      const auto &b = atoms[j];
      const auto &bj = basis_data(z[j]);
      const Vec delta = difference(position(xyz, i), position(xyz, j));
      lattice.images(delta, cutoff, [&](const auto &, const Vec &d, double r) {
        const bool onsite = r < 1.e-10;
        const int na = a.size(), nb = b.size();
        Array block(na * nb);
        std::array<double, 49> ov{}, scaled{};
        std::array<double, 3 * 49> dp{};
        std::array<double, 6 * 49> qp{};
        for (const auto &ga : a.scf) {
          for (const auto &gb : b.scf) {
            gxtb_moments(&gb, &ga, d.data(), ov.data(), dp.data(), qp.data());
            if (!onsite)
              gxtb_overlap(&b.h0[gb.l], &a.h0[ga.l], d.data(), scaled.data());
            for (int x = 0; x < 2 * ga.l + 1; ++x)
              for (int y = 0; y < 2 * gb.l + 1; ++y) {
                const int ij = y + (2 * gb.l + 1) * x, ia = ga.l * ga.l + x,
                          jb = gb.l * gb.l + y;
                const size_t idx = b.first + jb + n * (a.first + ia);
                s[idx] += ov[ij];
                for (int k = 0; k < 6; ++k)
                  quad[6 * idx + k] += qp[6 * ij + k];
                for (int k = 0; k < 3; ++k) {
                  dip[3 * idx + k] += dp[3 * ij + k];
                  h[idx] -= 0.5 * (anisotropy[i][k] * dp[3 * ij + k] +
                                   anisotropy[j][k] *
                                       (dp[3 * ij + k] + d[k] * ov[ij]));
                }
                block[jb + nb * ia] = onsite ? ov[ij] : scaled[ij];
              }
          }
          int offset = b.projector_start;
          for (int p = 0; p < bj.projectors; ++p) {
            const auto &proj = b.proj[p];
            gxtb_overlap(&proj, &ga, d.data(), ov.data());
            const int rows = 2 * proj.l + 1;
            for (int y = 0; y < rows; ++y)
              proj_level[offset + y] = bj.proj_level[p];
            for (int x = 0; x < 2 * ga.l + 1; ++x)
              for (int y = 0; y < rows; ++y)
                projector[offset + y +
                          size_t(np) * (a.first + ga.l * ga.l + x)] +=
                    ov[y + rows * x];
            offset += rows;
          }
        }
        if (!onsite) {
          Vec scale{};
          for (int k = 0; k < 3; ++k) {
            const double x = ai.diatomic[k], y = bj.diatomic[k];
            if (x * y != 0)
              scale[k] = 2 * x * y / (x + y);
          }
          gxtb_diatomic(bj.count - 1, ai.count - 1, d.data(), scale.data(),
                        block.data());
        }
        const double rr =
            r / (ai.h_radius + bj.h_radius) * (1. + 0.0402348406 * r);
        for (const auto &ga : a.scf)
          for (const auto &gb : b.scf) {
            double t = 0.5 * (a.level[ga.l] + b.level[gb.l]);
            if (!onsite) {
              t *= 0.5 * (wolfsberg[ga.l] + wolfsberg[gb.l]) *
                   (1. + (2. + ga.l) * ai.poly * rr) *
                   (1. + (2. + gb.l) * bj.poly * rr);
              if (ga.l + gb.l == 5)
                t *= 1.5;
            }
            for (int x = ga.l * ga.l; x < (ga.l + 1) * (ga.l + 1); ++x)
              for (int y = gb.l * gb.l; y < (gb.l + 1) * (gb.l + 1); ++y)
                h[b.first + y + n * (a.first + x)] += t * block[y + nb * x];
          }
      });
    }
  }
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j <= i; ++j) {
      double correction = 0;
      for (int k = 0; k < np; ++k)
        correction +=
            proj_level[k] * projector[k + np * i] * projector[k + np * j];
      h[i * n + j] += correction;
      if (i != j)
        h[j * n + i] += correction;
    }
  coulomb_integrals(nat, z, xyz, cell, orb, *this);
}
} // namespace gxtb
#endif
