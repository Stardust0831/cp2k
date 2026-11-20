/* SPDX-License-Identifier: GPL-2.0-or-later */
#include "../../src/gxtb_backend.h"
#include <stdbool.h>
#include <stdio.h>
#include <tblite.h>

const char *cp2k_gxtb_name_v1(void) {
  return "save_tblite electronics + CP2K native D4Srev";
}

/* tblite constructs its dispersion cache synchronously on the calling
 * thread. Each cache keeps the host callbacks for its entire lifetime. */
static _Thread_local const struct cp2k_gxtb_d4_v1 *host;
const struct cp2k_gxtb_d4_v1 *cp2k_save_d4_host(void) { return host; }
static _Thread_local cp2k_gxtb_scf_fn *electronic_scf;
cp2k_gxtb_scf_fn *cp2k_save_scf_solver(void) { return electronic_scf; }

static int evaluate(const struct cp2k_gxtb_d4_v1 *dispersion,
                    cp2k_gxtb_scf_fn *scf, int n, const int *z,
                    const double *xyz, const double *cell, const double *radii,
                    double *energy, double *gradient, double *strain,
                    char *message, size_t capacity) {
  if (radii[0] != 0 || radii[1] != 0) {
    snprintf(message, capacity,
             "save_tblite uses its own exchange kernel; omit EXCHANGE_RADII");
    return 1;
  }
  const struct cp2k_gxtb_d4_v1 *previous = host;
  cp2k_gxtb_scf_fn *previous_scf = electronic_scf;
  host = dispersion;
  electronic_scf = scf;
  const bool periodic[3] = {true, true, true};
  const int length = (int)capacity;
  tblite_error error = tblite_new_error();
  tblite_context ctx = tblite_new_context();
  tblite_result result = tblite_new_result();
  tblite_calculator calc = NULL;
  tblite_structure mol =
      tblite_new_structure(error, n, z, xyz, NULL, NULL, cell, periodic);
  int status = 1;
  if (!tblite_check_error(error)) {
    calc = tblite_new_gxtb_calculator(ctx, mol, NULL);
    if (calc && !tblite_check_context(ctx)) {
      tblite_set_context_verbosity(ctx, 0);
      tblite_set_calculator_accuracy(ctx, calc, 0.01);
      tblite_set_calculator_temperature(ctx, calc, 0.0);
      tblite_set_calculator_guess(ctx, calc, TBLITE_GUESS_EEQBC);
      tblite_get_singlepoint(ctx, mol, calc, result);
      if (!tblite_check_context(ctx)) {
        tblite_get_result_energy(error, result, energy);
        tblite_get_result_gradient(error, result, gradient);
        tblite_get_result_virial(error, result, strain);
        status = tblite_check_error(error);
      }
    }
  }
  if (status) {
    if (tblite_check_error(error))
      tblite_get_error(error, message, &length);
    else if (tblite_check_context(ctx))
      tblite_get_context_error(ctx, message, &length);
    else
      snprintf(message, capacity, "Could not construct the g-xTB calculator");
  }
  tblite_delete_calculator(&calc);
  tblite_delete_structure(&mol);
  tblite_delete_result(&result);
  tblite_delete_context(&ctx);
  tblite_delete_error(&error);
  host = previous;
  electronic_scf = previous_scf;
  return status;
}

int cp2k_gxtb_eval_v2(const struct cp2k_gxtb_d4_v1 *d4, int n, const int *z,
                      const double *xyz, const double *cell,
                      const double *radii, double *energy, double *gradient,
                      double *strain, char *message, size_t capacity) {
  return evaluate(d4, NULL, n, z, xyz, cell, radii, energy, gradient, strain,
                  message, capacity);
}

int cp2k_gxtb_eval_v3(cp2k_gxtb_scf_fn *scf, int n, const int *z,
                      const double *xyz, const double *cell,
                      const double *radii, double *energy, double *gradient,
                      double *strain, char *message, size_t capacity) {
  if (!scf) {
    snprintf(message, capacity, "Missing CP2K g-xTB electronic solver");
    return 1;
  }
  return evaluate(NULL, scf, n, z, xyz, cell, radii, energy, gradient, strain,
                  message, capacity);
}
