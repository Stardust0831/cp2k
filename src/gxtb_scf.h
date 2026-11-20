/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef CP2K_GXTB_SCF_H
#define CP2K_GXTB_SCF_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/* Column-major arrays owned by one electronic force evaluation. ao_shell
 * and shell_atom use one-based indices. Dipoles/quadrupoles have their
 * component index first. Spin moments are in charge/magnetization order;
 * P, C, eps and the nonlocal exchange operator are in alpha/beta order. */
struct cp2k_gxtb_scf_problem {
  int nao, nsh, nat, nspin, max_iter;
  double nel[2], energy_tol, residual_tol;
  const int *ao_shell, *shell_atom;
  const double *s, *h0, *dip_int, *quad_int, *neutral_shell;
  double *p, *c, *eps, *occ, *qsh, *qat, *dip, *quad;
  double *vat, *vsh, *vao, *vdip, *vquad, *exchange, *energy, *weighted_p;
  const int *z, *pbc;
  const double *xyz, *cell;
  double *geometry_gradient, *geometry_strain;
  const double *cn, *gamma, *tau, *sd, *dd, *sq, *dq, *qq;
};
typedef int cp2k_gxtb_scf_fn(struct cp2k_gxtb_scf_problem *,
                             size_t problem_size, char *message,
                             size_t capacity);
cp2k_gxtb_scf_fn cp2k_gxtb_scf;
#ifdef __cplusplus
}
#endif
#endif
