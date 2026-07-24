#!/bin/bash
set -euo pipefail

: "${CP2K_GPU_BINDING_DIR:?CP2K_GPU_BINDING_DIR is required}"

rank="${PMI_RANK:-${PMIX_RANK:-${OMPI_COMM_WORLD_RANK:-${SLURM_PROCID:-}}}}"
local_rank="${MPI_LOCALRANKID:-${PMIX_LOCAL_RANK:-${OMPI_COMM_WORLD_LOCAL_RANK:-${SLURM_LOCALID:-${rank}}}}}"
if [[ -z "${rank}" || -z "${local_rank}" ]]; then
  echo "ERROR: Could not determine the MPI rank for GPU binding" >&2
  exit 1
fi

device=$((local_rank % 2))
export CUDA_VISIBLE_DEVICES="${device}"
printf 'rank=%s local_rank=%s physical_device=%s host=%s\n' \
  "${rank}" "${local_rank}" "${device}" "$(hostname)" \
  > "${CP2K_GPU_BINDING_DIR}/rank-${rank}.log"

exec "$@"
