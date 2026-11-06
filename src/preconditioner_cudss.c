/*----------------------------------------------------------------------------*/
/*  CP2K: A general program to perform molecular dynamics simulations         */
/*  Copyright 2000-2026 CP2K developers group <https://cp2k.org>              */
/*                                                                            */
/*  SPDX-License-Identifier: GPL-2.0-or-later                                 */
/*----------------------------------------------------------------------------*/

#if defined(__CUDSS)
#include "offload/offload_library.h"
#include <cuda_runtime_api.h>
#include <cudss.h>
#include <limits.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

_Static_assert(CUDSS_VERSION_MAJOR >= 1,
               "The cuDSS preconditioner requires cuDSS 1.0 or newer");

#define CHECK(call)                                                            \
  do {                                                                         \
    int status = (int)(call);                                                  \
    if (status != 0) {                                                         \
      fprintf(stderr, "cuDSS preconditioner: %s failed (%d) at %s:%d\n",       \
              #call, status, __FILE__, __LINE__);                              \
      MPI_Abort(ctx->comm, status);                                            \
      abort();                                                                 \
    }                                                                          \
  } while (0)

typedef struct {
  int row, col;
  double value;
} entry_t;

typedef struct {
  MPI_Comm comm;
  int rank, size, n, device, batch_size;
  cudaStream_t stream;
  cudssHandle_t handle;
  cudssConfig_t config;
  cudssData_t data;
  cudssMatrix_t matrix;
  int *offsets, *columns;
  double *values, *rhs, *solution;
  double analysis_seconds, factor_seconds, solve_seconds;
  long solve_calls;
} preconditioner_t;

static void *allocate(size_t size, MPI_Comm comm) {
  void *p = malloc(size ? size : 1);
  if (!p) {
    fprintf(stderr, "cuDSS preconditioner: host allocation failed\n");
    MPI_Abort(comm, 1);
    abort();
  }
  return p;
}

static int compare_entries(const void *a, const void *b) {
  const entry_t *x = a, *y = b;
  if (x->row != y->row)
    return (x->row > y->row) - (x->row < y->row);
  return (x->col > y->col) - (x->col < y->col);
}

static void check_numerical_info(preconditioner_t *ctx) {
  int info = 0;
  size_t written = 0;
  CHECK(cudssDataGet(ctx->handle, ctx->data, CUDSS_DATA_INFO, &info,
                     sizeof(info), &written));
  if (info) {
    fprintf(stderr, "cuDSS preconditioner: SPD factorization/solve info=%d\n",
            info);
    MPI_Abort(ctx->comm, 1);
    abort();
  }
}

static void dense_descriptors(preconditioner_t *ctx, int nrhs, cudssMatrix_t *x,
                              cudssMatrix_t *b) {
  CHECK(cudssMatrixCreateDn(b, ctx->n, nrhs, ctx->n, ctx->rhs, CUDA_R_64F,
                            CUDSS_LAYOUT_COL_MAJOR));
  CHECK(cudssMatrixCreateDn(x, ctx->n, nrhs, ctx->n, ctx->solution, CUDA_R_64F,
                            CUDSS_LAYOUT_COL_MAJOR));
}

/* The prototype gathers sparse entries on rank zero; numerical factorization
 * and solves use cuDSS MGMN across the original CP2K communicator. */
void cp_cudss_create(void **context, int fortran_comm, int n, int local_nnz,
                     int batch_size, const int *local_rows,
                     const int *local_cols, const double *local_values) {
  MPI_Comm comm = MPI_Comm_f2c(fortran_comm);
  preconditioner_t *ctx = allocate(sizeof(*ctx), comm);
  *ctx = (preconditioner_t){0};
  ctx->comm = comm;
  CHECK(MPI_Comm_dup(comm, &ctx->comm));
  CHECK(MPI_Comm_rank(ctx->comm, &ctx->rank));
  CHECK(MPI_Comm_size(ctx->comm, &ctx->size));
  ctx->n = n;
  ctx->batch_size = batch_size;
  if (n < 1 || batch_size < 1 || n > INT_MAX / batch_size)
    MPI_Abort(ctx->comm, 1);
  ctx->device = offload_get_chosen_device();
  CHECK(cudaSetDevice(ctx->device));
  int *counts = allocate(ctx->size * sizeof(int), ctx->comm);
  int *displs = allocate(ctx->size * sizeof(int), ctx->comm);
  CHECK(MPI_Allgather(&local_nnz, 1, MPI_INT, counts, 1, MPI_INT, ctx->comm));
  int nnz = 0;
  for (int p = 0; p < ctx->size; p++) {
    if (counts[p] < 0 || counts[p] > INT_MAX - nnz) {
      fprintf(stderr,
              "cuDSS prototype: sparse entry count exceeds 32-bit indexing\n");
      MPI_Abort(ctx->comm, 1);
    }
    displs[p] = nnz;
    nnz += counts[p];
  }
  int *rows = NULL, *cols = NULL, *offsets = NULL;
  double *values = NULL;
  if (ctx->rank == 0) {
    rows = allocate((size_t)nnz * sizeof(int), ctx->comm);
    cols = allocate((size_t)nnz * sizeof(int), ctx->comm);
    values = allocate((size_t)nnz * sizeof(double), ctx->comm);
  }
  CHECK(MPI_Gatherv(local_rows, local_nnz, MPI_INT, rows, counts, displs,
                    MPI_INT, 0, ctx->comm));
  CHECK(MPI_Gatherv(local_cols, local_nnz, MPI_INT, cols, counts, displs,
                    MPI_INT, 0, ctx->comm));
  CHECK(MPI_Gatherv(local_values, local_nnz, MPI_DOUBLE, values, counts, displs,
                    MPI_DOUBLE, 0, ctx->comm));
  free(counts);
  free(displs);
  if (ctx->rank == 0) {
    entry_t *entries = allocate((size_t)nnz * sizeof(entry_t), ctx->comm);
    for (int i = 0; i < nnz; i++)
      entries[i] = (entry_t){rows[i], cols[i], values[i]};
    qsort(entries, nnz, sizeof(entry_t), compare_entries);
    offsets = allocate(((size_t)n + 1) * sizeof(int), ctx->comm);
    for (int i = 0; i <= n; i++)
      offsets[i] = 0;
    int count = 0;
    for (int i = 0; i < nnz; i++) {
      if (entries[i].row < 0 || entries[i].row >= n ||
          entries[i].col < entries[i].row || entries[i].col >= n) {
        MPI_Abort(ctx->comm, 1);
      }
      if (count && entries[i].row == rows[count - 1] &&
          entries[i].col == cols[count - 1]) {
        values[count - 1] += entries[i].value;
      } else {
        rows[count] = entries[i].row;
        cols[count] = entries[i].col;
        values[count++] = entries[i].value;
        offsets[entries[i].row + 1]++;
      }
    }
    nnz = count;
    for (int i = 0; i < n; i++)
      offsets[i + 1] += offsets[i];
    free(entries);
    free(rows);
  }
  CHECK(MPI_Bcast(&nnz, 1, MPI_INT, 0, ctx->comm));
  CHECK(cudaStreamCreate(&ctx->stream));
  CHECK(cudssCreate(&ctx->handle));
  CHECK(cudssSetStream(ctx->handle, ctx->stream));
  CHECK(cudssConfigCreate(&ctx->config));
  CHECK(cudssDataCreate(ctx->handle, &ctx->data));
  if (ctx->size > 1) {
    const char *layer = getenv("CUDSS_COMM_LIB");
    if (!layer)
      layer = "libcudss_commlayer_openmpi.so.1";
    CHECK(cudssSetCommLayer(ctx->handle, layer));
    CHECK(cudssDataSet(ctx->handle, ctx->data, CUDSS_DATA_COMM_DEVICE,
                       &ctx->comm, sizeof(ctx->comm)));
    CHECK(cudssDataSet(ctx->handle, ctx->data, CUDSS_DATA_COMM_HOST, &ctx->comm,
                       sizeof(ctx->comm)));
  }
  if (ctx->rank == 0) {
    CHECK(cudaMalloc((void **)&ctx->offsets, ((size_t)n + 1) * sizeof(int)));
    CHECK(cudaMalloc((void **)&ctx->columns, (size_t)nnz * sizeof(int)));
    CHECK(cudaMalloc((void **)&ctx->values, (size_t)nnz * sizeof(double)));
    CHECK(cudaMalloc((void **)&ctx->rhs,
                     (size_t)n * ctx->batch_size * sizeof(double)));
    CHECK(cudaMalloc((void **)&ctx->solution,
                     (size_t)n * ctx->batch_size * sizeof(double)));
    CHECK(cudaMemcpy(ctx->offsets, offsets, ((size_t)n + 1) * sizeof(int),
                     cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(ctx->columns, cols, (size_t)nnz * sizeof(int),
                     cudaMemcpyHostToDevice));
    CHECK(cudaMemcpy(ctx->values, values, (size_t)nnz * sizeof(double),
                     cudaMemcpyHostToDevice));
    free(offsets);
    free(cols);
    free(values);
  }
  CHECK(cudssMatrixCreateCsr(&ctx->matrix, n, n, nnz, ctx->offsets, NULL,
                             ctx->columns, ctx->values, CUDA_R_32I, CUDA_R_32I,
                             CUDA_R_64F, CUDSS_MTYPE_SPD, CUDSS_MVIEW_UPPER,
                             CUDSS_BASE_ZERO));
  cudssMatrix_t x, b;
  dense_descriptors(ctx, ctx->batch_size, &x, &b);
  double t = MPI_Wtime();
  CHECK(cudssExecute(ctx->handle, CUDSS_PHASE_ANALYSIS, ctx->config, ctx->data,
                     ctx->matrix, x, b));
  CHECK(cudaStreamSynchronize(ctx->stream));
  ctx->analysis_seconds = MPI_Wtime() - t;
  check_numerical_info(ctx);
  t = MPI_Wtime();
  CHECK(cudssExecute(ctx->handle, CUDSS_PHASE_FACTORIZATION, ctx->config,
                     ctx->data, ctx->matrix, x, b));
  CHECK(cudaStreamSynchronize(ctx->stream));
  ctx->factor_seconds = MPI_Wtime() - t;
  check_numerical_info(ctx);
  CHECK(cudssMatrixDestroy(x));
  CHECK(cudssMatrixDestroy(b));
  if (ctx->rank == 0) {
    printf(" CUDSS| FULL_KINETIC n=%d upper_nnz=%d ranks=%d rhs_batch=%d\n", n,
           nnz, ctx->size, ctx->batch_size);
    printf(" CUDSS| Analysis %.6f s factorization %.6f s\n",
           ctx->analysis_seconds, ctx->factor_seconds);
    fflush(stdout);
  }
  *context = ctx;
}

void cp_cudss_solve(void *context, int nrhs, double *values) {
  preconditioner_t *ctx = context;
  if (nrhs < 1 || nrhs > ctx->batch_size || ctx->n > INT_MAX / nrhs)
    MPI_Abort(ctx->comm, 1);
  CHECK(cudaSetDevice(ctx->device));
  const size_t bytes = (size_t)ctx->n * nrhs * sizeof(double);
  double t = MPI_Wtime();
  if (ctx->rank == 0)
    CHECK(cudaMemcpyAsync(ctx->rhs, values, bytes, cudaMemcpyHostToDevice,
                          ctx->stream));
  cudssMatrix_t x, b;
  dense_descriptors(ctx, nrhs, &x, &b);
  CHECK(cudssExecute(ctx->handle, CUDSS_PHASE_SOLVE, ctx->config, ctx->data,
                     ctx->matrix, x, b));
  if (ctx->rank == 0)
    CHECK(cudaMemcpyAsync(values, ctx->solution, bytes, cudaMemcpyDeviceToHost,
                          ctx->stream));
  CHECK(cudaStreamSynchronize(ctx->stream));
  check_numerical_info(ctx);
  CHECK(cudssMatrixDestroy(x));
  CHECK(cudssMatrixDestroy(b));
  CHECK(MPI_Bcast(values, ctx->n * nrhs, MPI_DOUBLE, 0, ctx->comm));
  ctx->solve_seconds += MPI_Wtime() - t;
  ctx->solve_calls++;
}

void cp_cudss_destroy(void **context) {
  preconditioner_t *ctx = *context;
  if (!ctx)
    return;
  int device;
  CHECK(cudaGetDevice(&device));
  CHECK(cudaSetDevice(ctx->device));
  CHECK(cudssMatrixDestroy(ctx->matrix));
  CHECK(cudssDataDestroy(ctx->handle, ctx->data));
  CHECK(cudssConfigDestroy(ctx->config));
  CHECK(cudssDestroy(ctx->handle));
  if (ctx->rank == 0) {
    CHECK(cudaFree(ctx->offsets));
    CHECK(cudaFree(ctx->columns));
    CHECK(cudaFree(ctx->values));
    CHECK(cudaFree(ctx->rhs));
    CHECK(cudaFree(ctx->solution));
    printf(" CUDSS| %ld solve batches, total solve/transfer %.6f s\n",
           ctx->solve_calls, ctx->solve_seconds);
    fflush(stdout);
  }
  CHECK(cudaStreamDestroy(ctx->stream));
  CHECK(cudaSetDevice(device));
  CHECK(MPI_Comm_free(&ctx->comm));
  free(ctx);
  *context = NULL;
}
#endif
