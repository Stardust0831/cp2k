#!/bin/bash
set -euo pipefail

: "${WORK_ROOT:?WORK_ROOT is required}"

readonly source_root="${WORK_ROOT}/source"
readonly build_root="${WORK_ROOT}/build-root"
readonly install_root="${WORK_ROOT}/install"
readonly results_root="${WORK_ROOT}/results"
readonly cuda_root="/opt/devtools/nvidia/cuda-12.4.1"
readonly nccl_root="/opt/devtools/nvidia/nccl_2.18.5_cuda12.4_sai_v2.18.5-1-sai.1"

mkdir -p "${build_root}/spack" "${results_root}"
if [[ ! -x "${build_root}/spack/spack/bin/spack" ]]; then
  tar -xzf "${WORK_ROOT}/spack.tar.gz" -C "${build_root}/spack"
  mv "${build_root}/spack/spack-1.2.1" "${build_root}/spack/spack"
fi
mkdir -p "${WORK_ROOT}/spack-packages"
tar -xzf "${WORK_ROOT}/spack-packages.tar.gz" -C "${WORK_ROOT}/spack-packages"

readonly openblas_patch_sha="723ddc1553b6d27ff89d96985f7732695935c0d4d8df766987702689bdb750ac"
readonly source_cache="${build_root}/spack/spack/var/spack/cache/_source-cache/archive"
mkdir -p "${source_cache}/${openblas_patch_sha:0:2}"
cp "${WORK_ROOT}/${openblas_patch_sha}" \
  "${source_cache}/${openblas_patch_sha:0:2}/${openblas_patch_sha}"

# The compute nodes cannot reach GitHub. Point the pinned builtin repository at
# the checkout staged by the Actions driver; package sources still use Spack's
# public source mirror and NVIDIA's redistribution server.
perl -0pi -e \
  's|  repos:\n    builtin:\n      commit: [^\n]+\n|  repos:\n    builtin: '"${WORK_ROOT}"'/spack-packages/repos/spack_repo/builtin\n|' \
  "${source_root}/tools/spack/cp2k_deps_p.yaml"

export BUILD_PATH="${build_root}"
export CP2K_ROOT="${source_root}"
export CUDA_HOME="${cuda_root}"
export INSTALL_PREFIX="${install_root}"
export PATH="${cuda_root}/bin:${PATH}"
export LD_LIBRARY_PATH="${nccl_root}/lib:${cuda_root}/targets/x86_64-linux/lib:${LD_LIBRARY_PATH:-}"

echo "=== Hardware ==="
hostname
nvidia-smi --query-gpu=index,name,compute_cap --format=csv,noheader
echo "=== CP2K revision ==="
printf 'CI revision: %s\n' "${SOURCE_SHA:-staged-source}"
printf 'PR source revision: %s\n' "${PR_SOURCE_SHA:-unknown}"
echo "=== Rewritten Spack repository configuration ==="
sed -n '/^  repos:/,/^  specs:/p' "${source_root}/tools/spack/cp2k_deps_p.yaml"

cd "${source_root}"
# Recreate the environment from this revision while retaining already installed
# packages in Spack's store from an earlier attempt.
# shellcheck source=/dev/null
source "${build_root}/spack/spack/share/spack/setup-env.sh"
spack env remove --yes-to-all cp2k_env 2>/dev/null || true
spack clean --failures
# Force make_cp2k.sh to regenerate the environment with the CUDA-compatible
# compiler selection while retaining Spack's installed package store.
completion_marker="${build_root}/spack/BUILD_DEPENDENCIES_COMPLETED"
if [[ -e "${completion_marker}" ]]; then
  unlink "${completion_marker}"
fi
./make_cp2k.sh --cp2k_version psmp --mpi_mode mpich --gpu_model V100 \
  --gcc_version 13 --rebuild_cp2k \
  --disable_feature all --enable_feature cusolver_mp --use_cache no \
  --install_path "${install_root}" --num_packages 2 -j 16 \
  2>&1 | tee "${results_root}/build.log"

# shellcheck source=/dev/null
source "${build_root}/spack/spack/share/spack/setup-env.sh"
spack -e cp2k_env find -dl | tee "${results_root}/spack-dag.log"
spack -e cp2k_env find --format '{name}' | sort -u |
  tee "${results_root}/spack-package-names.log"
grep -qx cusolvermp "${results_root}/spack-package-names.log"
grep -qx nccl "${results_root}/spack-package-names.log"
if grep -qx ucc "${results_root}/spack-package-names.log"; then
  echo "ERROR: UCC is present in the concretized Spack environment"
  exit 1
fi
if grep -qx libxstream "${results_root}/spack-package-names.log"; then
  echo "ERROR: libxstream is present in the CUDA-only Spack environment"
  exit 1
fi

# Activate the environment created by make_cp2k.sh and add CP2K's installed
# shared library to the runtime search path.
eval "$(spack env activate --sh cp2k_env)"
export LD_LIBRARY_PATH="${install_root}/lib:${LD_LIBRARY_PATH:-}"
readonly cp2k_bin="${install_root}/bin/cp2k.psmp"
"${cp2k_bin}" --version 2>&1 | tee "${results_root}/cp2k-version.log"
grep -q 'cusolvermp_nccl' "${results_root}/cp2k-version.log"

ldd "${cp2k_bin}" | tee "${results_root}/cp2k-ldd.log"
readonly cusolvermp_view="${build_root}/spack/spack/opt/spack/view"
cusolvermp_lib="$(find "${cusolvermp_view}" -name 'libcusolverMp.so.0' -print -quit)"
readonly cusolvermp_lib
: "${cusolvermp_lib:?libcusolverMp.so.0 was not found in the Spack view}"
ldd "${cusolvermp_lib}" | tee "${results_root}/cusolvermp-ldd.log"
if grep -Eqi 'libucc' "${results_root}/cp2k-ldd.log" \
  "${results_root}/cusolvermp-ldd.log"; then
  echo "ERROR: UCC is present in the runtime dependency graph"
  exit 1
fi
if grep -Eqi 'libxstream|libOpenCL' "${results_root}/cp2k-ldd.log"; then
  echo "ERROR: CUDA-only CP2K unexpectedly links against OpenCL/libxstream"
  exit 1
fi
if grep -Eqi 'libucs|libucp' "${results_root}/cusolvermp-ldd.log"; then
  echo "ERROR: cuSOLVERMp unexpectedly links against UCX"
  exit 1
fi

cd "${source_root}/tests/QS/regtest-cusolver"
export CUDA_VISIBLE_DEVICES=0,1
for input in Si8-generalized.inp Si8-generalized-complex.inp; do
  sed '/^&GLOBAL$/a\  &TIMINGS\n    THRESHOLD 0.0\n  &END TIMINGS' \
    "${input}" > "${results_root}/${input}"
done

export CP2K_GPU_BINDING_DIR="${results_root}/bindings-real"
mkdir -p "${CP2K_GPU_BINDING_DIR}"
mpiexec -n 2 "${source_root}/.github/scripts/cp2k_gpu_rank_wrapper.sh" \
  "${cp2k_bin}" -i "${results_root}/Si8-generalized.inp" \
  -o "${results_root}/Si8-generalized.out"
cat "${CP2K_GPU_BINDING_DIR}"/rank-*.log | sort |
  tee "${results_root}/gpu-bindings-real.log"

export CP2K_GPU_BINDING_DIR="${results_root}/bindings-complex"
mkdir -p "${CP2K_GPU_BINDING_DIR}"
mpiexec -n 2 "${source_root}/.github/scripts/cp2k_gpu_rank_wrapper.sh" \
  "${cp2k_bin}" -i "${results_root}/Si8-generalized-complex.inp" \
  -o "${results_root}/Si8-generalized-complex.out"
cat "${CP2K_GPU_BINDING_DIR}"/rank-*.log | sort |
  tee "${results_root}/gpu-bindings-complex.log"

for binding_log in "${results_root}"/gpu-bindings-*.log; do
  grep -q 'rank=0 .*physical_device=0 ' "${binding_log}"
  grep -q 'rank=1 .*physical_device=1 ' "${binding_log}"
done
grep -q 'cp_fm_general_cusolver' "${results_root}/Si8-generalized.out"
grep -q 'cp_cfm_general_cusolver' "${results_root}/Si8-generalized-complex.out"

grep -F 'ENERGY| Total FORCE_EVAL' "${results_root}/Si8-generalized.out" | tail -1 |
  tee "${results_root}/energies.log"
grep -F 'ENERGY| Total FORCE_EVAL' "${results_root}/Si8-generalized-complex.out" | tail -1 |
  tee -a "${results_root}/energies.log"

awk 'NR == 1 {d = $NF + 31.187602969867214; if (d < 0) d = -d; if (d > 1e-11) exit 1}
     NR == 2 {d = $NF + 31.461089087225638; if (d < 0) d = -d; if (d > 1e-10) exit 1}' \
  "${results_root}/energies.log"

echo "PASS: cuSOLVERMp 0.7.2 completed real and complex multi-rank tests without UCC."
