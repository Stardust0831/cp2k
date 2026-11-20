/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef CP2K_GXTB_BACKEND_H
#define CP2K_GXTB_BACKEND_H
#include "gxtb_revd4.h"
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/* Neutral 3D Gamma supercell. Each consecutive triple is one position or
 * lattice vector in bohr. Radii are in bohr; zeros select provider defaults.
 * Return energy in Ha, dE/dr in Ha/bohr and dE/d(strain) in Ha (column-major).
 * A nonzero return code aborts the calculation with the supplied message. */
typedef int cp2k_gxtb_fn(int n, const int *z, const double *positions,
                         const double *lattice, const double *radii,
                         double *energy, double *gradient, double *strain,
                         char *message, size_t capacity);
cp2k_gxtb_fn cp2k_gxtb_eval_v1;
typedef int cp2k_gxtb_fn_v2(const struct cp2k_gxtb_d4_v1 *dispersion, int n,
                            const int *z, const double *positions,
                            const double *lattice, const double *radii,
                            double *energy, double *gradient, double *strain,
                            char *message, size_t capacity);
cp2k_gxtb_fn_v2 cp2k_gxtb_eval_v2;
const char *cp2k_gxtb_name_v1(void);
#ifdef __cplusplus
}
#endif
#endif
