#!/bin/bash
set -euo pipefail

: "${GPU_SSH_HOST:?CP2K_GPU_SSH_HOST secret is required}"
: "${GITHUB_RUN_ID:?GITHUB_RUN_ID is required}"

readonly remote_root="${GPU_REMOTE_ROOT:-cp2k-cusolvermp-actions-${GITHUB_RUN_ID}}"
readonly assets_dir="${RUNNER_TEMP}/cp2k-cusolvermp-assets-${GITHUB_RUN_ID}"
readonly spack_version="1.2.1"
readonly spack_commit="131214174051056d434e43d34f7c7645a385a835"
readonly openblas_patch_sha="723ddc1553b6d27ff89d96985f7732695935c0d4d8df766987702689bdb750ac"
readonly ssh_control_path="${RUNNER_TEMP}/cp2k-ssh-${GITHUB_RUN_ID}"
readonly -a ssh_options=(
  -o ConnectTimeout=30
  -o ControlMaster=auto
  -o ControlPersist=60
  -o "ControlPath=${ssh_control_path}"
  -o ServerAliveInterval=30
  -o ServerAliveCountMax=3
)
remote_ready=0

retry() {
  local attempt
  for attempt in 1 2 3; do
    if "$@"; then
      return
    fi
    echo "Attempt ${attempt} failed: $*" >&2
    if ((attempt < 3)); then
      sleep 10
    fi
  done
  return 1
}

collect_results() {
  if ((remote_ready)); then
    for attempt in 1 2 3; do
      if scp "${ssh_options[@]}" -r \
        "${GPU_SSH_HOST}:${remote_root}/results/." artifacts/; then
        return
      fi
      echo "Result collection attempt ${attempt} failed" >&2
      sleep 5
    done
  fi
}
trap collect_results EXIT

mkdir -p "${assets_dir}" artifacts

curl --fail --location --retry 3 \
  "https://github.com/spack/spack/archive/refs/tags/v${spack_version}.tar.gz" \
  --output "${assets_dir}/spack.tar.gz"
curl --fail --location --retry 3 \
  "https://github.com/OpenMathLib/OpenBLAS/commit/88705a932831c0de1ed136b461c6c239802828b2.diff?full_index=1" \
  --output "${assets_dir}/${openblas_patch_sha}"
printf '%s  %s\n' "${openblas_patch_sha}" "${assets_dir}/${openblas_patch_sha}" |
  sha256sum --check

git clone --filter=blob:none --no-checkout \
  https://github.com/spack/spack-packages.git "${assets_dir}/spack-packages-git"
git -C "${assets_dir}/spack-packages-git" checkout "${spack_commit}"
git -C "${assets_dir}/spack-packages-git" archive --format=tar.gz \
  --output="${assets_dir}/spack-packages.tar.gz" "${spack_commit}"

tar --exclude=.git --exclude=artifacts --format=posix -czf \
  "${assets_dir}/cp2k-source.tar.gz" .

retry ssh "${ssh_options[@]}" -tt "${GPU_SSH_HOST}" \
  "mkdir -p \"\$HOME/${remote_root}/source\""
remote_ready=1
retry scp "${ssh_options[@]}" \
  "${assets_dir}/spack.tar.gz" "${assets_dir}/spack-packages.tar.gz" \
  "${assets_dir}/cp2k-source.tar.gz" "${assets_dir}/${openblas_patch_sha}" \
  "${GPU_SSH_HOST}:${remote_root}/"
retry ssh "${ssh_options[@]}" -tt "${GPU_SSH_HOST}" \
  "tar -xzf \"\$HOME/${remote_root}/cp2k-source.tar.gz\" -C \"\$HOME/${remote_root}/source\""

ssh "${ssh_options[@]}" -tt "${GPU_SSH_HOST}" \
  "srun --partition=16V100 --qos=flood-1o2gpu --nodes=1 --ntasks=1 \
   --gres=gpu:V100-SXM2:2 --time=02:30:00 --unbuffered \
   env WORK_ROOT=\"\$HOME/${remote_root}\" SOURCE_SHA=\"${GITHUB_SHA}\" \
   bash \"\$HOME/${remote_root}/source/.github/scripts/run_cusolvermp_no_ucc_v100.sh\""
