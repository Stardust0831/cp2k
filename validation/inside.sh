#!/usr/bin/env bash
set -eo pipefail
validation_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
dependency_prefix=$(python3 -c 'import json,sys; print(json.load(open(sys.argv[1]))["dependency_prefix"])' "$validation_root/request.json")
driver_ld=${LD_LIBRARY_PATH:-}
source "$dependency_prefix/share/sai/runtime-env.sh"
export CP2K_DATA_DIR="$validation_root/source/data"
export LD_LIBRARY_PATH="$validation_root/build/src:$LD_LIBRARY_PATH${driver_ld:+:$driver_ld}"
if [[ "$1" == build ]]; then
  shift
  export GIT_CEILING_DIRECTORIES="$validation_root"
  source /etc/profile.d/lmod.sh
  module purge
  module use /opt/modules/modulefiles/devtools
  module load cmake/3.31.6 openmpi/5.0.10-nvhpc26.3-gnu-cuda12-auto fftw/3.3.10 saiblas/2603-gnu-auto cuda/12.9.1 nvmplibs/26.7-tmp gcc/14.3.0 spglib/2.7.0
  export CC=/opt/devtools/gcc/14.3.0/bin/gcc CXX=/opt/devtools/gcc/14.3.0/bin/g++ FC=/opt/devtools/gcc/14.3.0/bin/gfortran
  export OMPI_CC="$CC" OMPI_CXX="$CXX" OMPI_FC="$FC"
  exec "$@"
fi
if [[ "$1" == inspect ]]; then
  shift
  exec "$@"
fi
if [[ "$1" == lifecycle ]]; then
  shift
  if [[ "${SAI_CUSOLVER_COUNT:-0}" == 1 ]]; then
    export LD_PRELOAD="$validation_root/bin/cusolver_count.so"
  fi
  exec "$validation_root/bin/lifecycle" "$@"
fi
[[ "$1" == run ]]
shift
exec "$validation_root/build/bin/cp2k.psmp" "$@"
