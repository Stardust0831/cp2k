/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "gxtb_backend.h"
#include <dlfcn.h>
#include <stdio.h>

extern cp2k_gxtb_fn atomli_evaluate;

const char *cp2k_gxtb_name_v1(void) { return "Atomli " ATOMLI_VERSION " CPU g-xTB"; }

int cp2k_gxtb_eval_v1(int n, const int *z, const double *xyz,
                      const double *cell, const double *radii, double *energy,
                      double *gradient, double *strain, char *message,
                      size_t capacity) {
  /* Python extension modules resolve their interpreter symbols globally. */
  void *python = dlopen(LIBPYTHON_NAME, RTLD_NOW | RTLD_GLOBAL);
  if (!python) {
    snprintf(message, capacity, "%s", dlerror());
    return 1;
  }
  int status = atomli_evaluate(n, z, xyz, cell, radii, energy, gradient, strain,
                              message, capacity);
  dlclose(python);
  /* CFFI returns zero when the interpreter or callback could not start. */
  if (status == 0)
    snprintf(message, capacity, "Atomli Python initialization failed; see stderr");
  return status != 1;
}
