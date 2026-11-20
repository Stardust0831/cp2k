/* SPDX-License-Identifier: GPL-2.0-or-later */
#ifndef CP2K_GXTB_REVD4_H
#define CP2K_GXTB_REVD4_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif

/* Atomic units. Positions and lattice vectors are consecutive triples.
 * A geometry belongs to its caller until destroy. The two SCF outputs
 * are atom energies and derivatives with respect to Mulliken charges. */
struct cp2k_gxtb_d4_v1 {
  void *(*create)(int, const int *, const double *, const double *, const int *,
                  char *, size_t);
  void (*destroy)(void *);
  void (*scf)(const void *, const double *, double *, double *);
  void (*derivatives)(const void *, const double *, double *, double *);
  void (*atm)(const void *, double *, double *, double *);
};
const struct cp2k_gxtb_d4_v1 *cp2k_gxtb_native_d4(void);

#ifdef __cplusplus
}
#endif
#endif
