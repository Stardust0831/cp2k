/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "gxtb_parameters.hpp"
#include <algorithm>
#include <stdexcept>

namespace gxtb {
namespace {
#include "gxtb_electronic_data.inc"
}
const Element &element(int z) {
  if (z < 1 || z > 103)
    throw std::invalid_argument("g-xTB requires Z=1..103");
  return elements[z - 1];
}
double d3_radius_angstrom(int z, int partner) {
  const auto [lo, hi] = std::minmax(z, partner);
  return d3_pair_radius[hi * (hi - 1) / 2 + lo - 1];
}
double pair_radius(int z, int partner) {
  return d3_radius_angstrom(z, partner) * 0.5 *
         (element(z).vdw_scale + element(partner).vdw_scale);
}
} // namespace gxtb
