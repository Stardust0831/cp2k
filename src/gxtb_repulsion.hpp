/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef CP2K_GXTB_REPULSION_HPP
#define CP2K_GXTB_REPULSION_HPP
#include "gxtb_pairs.hpp"

namespace gxtb {
class Repulsion {
public:
  Repulsion(int n, const int *z, const double *xyz, const double *cell,
            const int *pbc);
  void add_core(double *energy) const;
  void add_scf(const double *charge, double *energy, double *potential) const;
  void gradient(const double *charge, double *grad, double *strain) const;

private:
  struct Edge : Pair {
    double dcn_dr, kernel, dk_dr, dk_dcn_i, dk_dcn_j;
  };
  struct Coefficient {
    int i, j;
    double kernel;
  };
  std::vector<int> species;
  std::vector<Edge> edges;
  std::vector<Coefficient> coefficients;
  double effective_charge(size_t i, const double *charge) const;
};
} // namespace gxtb
#endif
