/* SPDX-License-Identifier: LGPL-3.0-or-later */
#include "../../../src/gxtb_revd4.h"
#include <stdio.h>
#include <stdlib.h>

struct cache {
  const struct cp2k_gxtb_d4_v1 *ops;
  void *geometry;
};

#ifdef GXTB_STANDALONE
#define cp2k_save_d4_host cp2k_gxtb_native_d4
#else
const struct cp2k_gxtb_d4_v1 *cp2k_save_d4_host(void);
#endif

void *cp2k_save_d4_new(int n, const int *z, const double *xyz,
                       const double *cell, const int *pbc, char *msg,
                       size_t len) {
  const struct cp2k_gxtb_d4_v1 *ops = cp2k_save_d4_host();
  if (!ops) {
    snprintf(msg, len, "Missing CP2K D4Srev host");
    return NULL;
  }
  struct cache *c = malloc(sizeof(*c));
  if (!c) {
    snprintf(msg, len, "Could not allocate the D4Srev cache");
    return NULL;
  }
  c->ops = ops;
  c->geometry = ops->create(n, z, xyz, cell, pbc, msg, len);
  if (!c->geometry) {
    free(c);
    return NULL;
  }
  return c;
}
void cp2k_save_d4_free(struct cache *c) {
  c->ops->destroy(c->geometry);
  free(c);
}
void cp2k_save_d4_scf(const struct cache *c, const double *q, double *e,
                      double *v) {
  c->ops->scf(c->geometry, q, e, v);
}
void cp2k_save_d4_derivatives(const struct cache *c, const double *q, double *g,
                              double *s) {
  c->ops->derivatives(c->geometry, q, g, s);
}
void cp2k_save_d4_atm(const struct cache *c, double *e, double *g, double *s) {
  c->ops->atm(c->geometry, e, g, s);
}
