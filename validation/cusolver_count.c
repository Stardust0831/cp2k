/* External lifecycle-only interposer; the timed benchmarks do not load it. */
#define _GNU_SOURCE
#include <cusolverMp.h>
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>

static unsigned long created, destroyed;

cusolverStatus_t CUSOLVERAPI cusolverMpCreate(cusolverMpHandle_t *handle,
                                            int device, cudaStream_t stream) {
  typedef cusolverStatus_t (*function_t)(cusolverMpHandle_t *, int, cudaStream_t);
  function_t next = (function_t)dlsym(RTLD_NEXT, "cusolverMpCreate");
  if (!next)
    abort();
  cusolverStatus_t status = next(handle, device, stream);
  if (status == CUSOLVER_STATUS_SUCCESS)
    ++created;
  return status;
}

cusolverStatus_t CUSOLVERAPI cusolverMpDestroy(cusolverMpHandle_t handle) {
  typedef cusolverStatus_t (*function_t)(cusolverMpHandle_t);
  function_t next = (function_t)dlsym(RTLD_NEXT, "cusolverMpDestroy");
  if (!next)
    abort();
  cusolverStatus_t status = next(handle);
  if (status == CUSOLVER_STATUS_SUCCESS)
    ++destroyed;
  return status;
}

__attribute__((destructor)) static void report(void) {
  const char *rank = getenv("OMPI_COMM_WORLD_RANK");
  if (!rank)
    abort();
  char path[128];
  snprintf(path, sizeof(path), "cusolver-context-counts-rank-%s.json", rank);
  FILE *f = fopen(path, "w");
  if (!f)
    abort();
  fprintf(f, "{\"rank\":%s,\"created\":%lu,\"destroyed\":%lu}\n", rank,
          created, destroyed);
  fclose(f);
}
