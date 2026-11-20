/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef CP2K_GXTB_ONECENTER_HPP
#define CP2K_GXTB_ONECENTER_HPP

namespace gxtb {
struct ShellData {
  int l, n;
  double alpha[12], c0[12], c1[12];
};
struct BasisData {
  int start, count, projectors;
  double charge[3], scale[3], cn_radius, h_radius, poly, anisotropy;
  double level[4], cn_level[4], h_exp[4], diatomic[3], occ[4];
  int proj_l[4];
  double proj_level[4], proj_exp[4];
  double hardness[4], hardness_cn, dipole_cn;
  struct {
    double chi, eta, width, kcn, kq, cap, cov, avg_cn, vdw, en;
  } eeq;
};
const BasisData &basis_data(int z);
} // namespace gxtb
#endif
