/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef CP2K_GXTB_EXCHANGE_HPP
#define CP2K_GXTB_EXCHANGE_HPP
#include "gxtb_pairs.hpp"
#include "gxtb_scf.h"

namespace gxtb {
class Exchange {
  struct ShellPair {
    double mulliken = 0, bond = 0, onsite = 0, rotation = 0;
  };
  struct Radial {
    Pair image;
    size_t shells;
    double dmulliken, dbond;
  };
  int ns;
  std::vector<int> shell, atom;
  std::vector<double> kq;
  std::vector<ShellPair> kernels;
  std::vector<Radial> radial;

public:
  explicit Exchange(const cp2k_gxtb_scf_problem &p);
  // At the final density also differentiate the kernels and AO overlap.
  void evaluate(cp2k_gxtb_scf_problem &p, double *energy,
                bool derivatives = false) const;
};
} // namespace gxtb
#endif
