/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "gxtb_backend.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifndef _WIN32
#include <dlfcn.h>
#endif

struct gxtb_handle {
  void *library;
  cp2k_gxtb_fn *evaluate;
};

void *cp_gxtb_open(const char *path, char *message, size_t capacity) {
#ifdef _WIN32
  snprintf(message, capacity, "The g-xTB loader requires POSIX shared libraries");
  return NULL;
#else
  void *library = dlopen(path, RTLD_NOW | RTLD_LOCAL);
  if (!library) {
    snprintf(message, capacity, "%s", dlerror());
    return NULL;
  }
  const char *(*name)(void) = dlsym(library, "cp2k_gxtb_name_v1");
  cp2k_gxtb_fn *evaluate = dlsym(library, "cp2k_gxtb_eval_v1");
  if (!name || !evaluate) {
    snprintf(message, capacity, "Library is missing its g-xTB v1 entry points");
    dlclose(library);
    return NULL;
  }
  struct gxtb_handle *handle = malloc(sizeof(*handle));
  if (!handle) {
    snprintf(message, capacity, "Could not allocate the g-xTB handle");
    dlclose(library);
    return NULL;
  }
  handle->library = library;
  handle->evaluate = evaluate;
  snprintf(message, capacity, "%s", name());
  return handle;
#endif
}

int cp_gxtb_eval(const struct gxtb_handle *handle, int n, const int *z,
                 const double *pos, const double *cell, const double *radii,
                 double *energy, double *grad, double *strain, char *message,
                 size_t capacity) {
  /* A successful call must set all outputs, including zero derivatives. */
  *energy = NAN;
  for (size_t i = 0; i < 3 * (size_t)n; i++) grad[i] = NAN;
  for (size_t i = 0; i < 9; i++) strain[i] = NAN;
  message[0] = '\0';
  int status = handle->evaluate(n, z, pos, cell, radii, energy, grad, strain,
                                message, capacity);
  message[capacity - 1] = '\0';
  return status;
}

void cp_gxtb_close(struct gxtb_handle *handle) {
  if (!handle) return;
#ifndef _WIN32
  dlclose(handle->library);
#endif
  free(handle);
}
