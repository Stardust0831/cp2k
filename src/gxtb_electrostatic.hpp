/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef CP2K_GXTB_ELECTROSTATIC_HPP
#define CP2K_GXTB_ELECTROSTATIC_HPP
#include "gxtb_scf.h"
#include <vector>

namespace gxtb {
class Electrostatic {
  std::vector<double> ipea, third, fourth, reference_field, spin_matrix;

public:
  explicit Electrostatic(const cp2k_gxtb_scf_problem &p);
  void evaluate(cp2k_gxtb_scf_problem &p, double *energy) const;
};
} // namespace gxtb
#endif
