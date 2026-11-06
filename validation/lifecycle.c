/* End-to-end libcp2k test: repeated real/complex solves and BLACS lifetimes. */
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include "libcp2k.h"

static void exercise(const char *name, MPI_Comm comm, int group,
                     const char *phase) {
  int rank, size, natom;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  char input[256], output[256], report[256];
  snprintf(input, sizeof(input), "%s.inp", name);
  snprintf(output, sizeof(output), "%s-%s-g%d.out", phase, name, group);
  snprintf(report, sizeof(report), "observables-g%d.jsonl", group);
  force_env_t env;
  cp2k_create_force_env_comm(&env, input, output, (int)MPI_Comm_c2f(comm));
  cp2k_get_natom(env, &natom);
  if (natom != 8)
    MPI_Abort(MPI_COMM_WORLD, 10);
  double positions[24], forces[24], energy;
  cp2k_get_positions(env, positions, 24);
  for (int step = 0; step < 2; ++step) {
    /* Changing S and H between calls catches accidental retention of matrices. */
    if (step == 1) {
      positions[0] += 0.01;
      cp2k_set_positions(env, positions, 24);
    }
    cp2k_calc_energy_force(env);
    cp2k_get_potential_energy(env, &energy);
    cp2k_get_forces(env, forces, 24);
    if (!isfinite(energy))
      MPI_Abort(MPI_COMM_WORLD, 11);
    for (int i = 0; i < 24; ++i)
      if (!isfinite(forces[i]))
        MPI_Abort(MPI_COMM_WORLD, 12);
    if (rank == 0) {
      FILE *f = fopen(report, "a");
      if (!f)
        MPI_Abort(MPI_COMM_WORLD, 13);
      fprintf(f, "{\"phase\":\"%s\",\"case\":\"%s\",\"group\":%d,"
                 "\"ranks\":%d,\"step\":%d,\"energy_ha\":%.17g,\"forces\":[",
              phase, name, group, size, step, energy);
      for (int i = 0; i < 24; ++i)
        fprintf(f, "%s%.17g", i ? "," : "", forces[i]);
      fprintf(f, "]}\n");
      fclose(f);
    }
  }
  cp2k_destroy_force_env(env);
}

int main(int argc, char **argv) {
  int rank, size, provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size != 4 || provided < MPI_THREAD_MULTIPLE)
    MPI_Abort(MPI_COMM_WORLD, 1);
  cp2k_init_without_mpi();
  const char *cases[] = {"real-standard", "real-subset", "real-partial",
                         "real-full", "complex-subset", "complex-full"};
  for (int i = 0; i < 6; ++i)
    exercise(cases[i], MPI_COMM_WORLD, 0, "world");

  /* Concurrent independent grids, then destruction and world-grid recreation. */
  MPI_Comm subgroup;
  int group = rank / 2;
  MPI_Comm_split(MPI_COMM_WORLD, group, rank, &subgroup);
  exercise(group == 0 ? "real-subset" : "complex-subset", subgroup,
           group, "split");
  MPI_Comm_free(&subgroup);
  MPI_Barrier(MPI_COMM_WORLD);
  exercise("real-subset", MPI_COMM_WORLD, 0, "recreated");
  cp2k_finalize_without_mpi();
  if (rank == 0)
    printf("CP2K_CUSOLVER_LIFECYCLE_COMPLETE\n");
  MPI_Finalize();
  return 0;
}
