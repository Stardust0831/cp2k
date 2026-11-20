/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef CP2K_GXTB_INTEGRALS_HPP
#define CP2K_GXTB_INTEGRALS_HPP
#include "gxtb_pairs.hpp"
namespace gxtb {
using Array = std::vector<double>;
struct Orbitals {
  std::vector<int> ao_shell, shell_atom, start;
  Array neutral;
  Orbitals(int nat, const int *z);
};
struct Matrices {
  Array s, h, dip, quad, cn, aux, gamma, tau, sd, dd, sq, dq, qq;
  Matrices(int nat, const int *z, const double *xyz, const double *cell,
           const Orbitals &orb);
};
void coulomb_integrals(int nat, const int *z, const double *xyz,
                       const double *cell, const Orbitals &orb, Matrices &m);
Array auxiliary_charges(int n, const int *z, const double *xyz,
                        const double *cell);
} // namespace gxtb
#endif
