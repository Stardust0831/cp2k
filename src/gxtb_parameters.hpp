/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef CP2K_GXTB_PARAMETERS_HPP
#define CP2K_GXTB_PARAMETERS_HPP
namespace gxtb {
struct Element {
  double core, radius, alpha, kcn, zeff, kq, offset, k1, vdw_scale;
  double bond_scale, bond_radius, hardness[4], average_power[4], onsite[10];
  double ipea_cn, hubbard_derivative, ipea[4], reference_charge[4], spin[10];
};
const Element &element(int z);
double pair_radius(int z, int partner);
double d3_radius_angstrom(int z, int partner);
} // namespace gxtb
#endif
