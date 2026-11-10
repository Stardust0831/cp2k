#!/bin/bash
# shellcheck disable=SC1091,SC2016
set -eo pipefail
export PATH=/usr/local/bin:/usr/bin:/bin SHELL=/bin/bash

task=$1
support=$2
source /etc/profile.d/lmod.sh
module purge
module use /opt/modules/modulefiles/devtools
module load apptainer/1.4.4 openmpi/5.0.10-nvhpc26.3-gnu-cuda12-auto
cd "$task"
scontrol show job -o "$SLURM_JOB_ID" > allocation.txt
nvidia-smi -L > devices.txt
image=$(python3 -c 'import json; print(json.load(open("request.json"))["environment"]["image"])')
prefix=$(python3 -c 'import json; print(json.load(open("request.json"))["environment"]["prefix"])')
build_tools=$(python3 -c 'import json; print(json.load(open("request.json"))["environment"]["build_tools"])')
auxiliary_cache=$(python3 -c 'import json; print(json.load(open("request.json"))["environment"]["auxiliary_cache"])')
mkdir -p build build-dependencies tmp apptainer-cache
export APPTAINER_TMPDIR="$task/tmp" APPTAINER_CACHEDIR="$task/apptainer-cache"
args=(apptainer exec --nv --cleanenv --no-home --no-mount 'bind-paths,home,cwd,tmp,hostfs')
for path in /usr /lib /lib64 /etc /opt/devtools /opt/modules; do
  args+=(--bind "$path:$path:ro")
done
args+=(--bind "$task:$task:rw" --bind "$task/source:$task/source:ro"
  --bind "$task/build-dependencies:$prefix/dependencies:rw"
  --bind "$build_tools:/sai-build-tools:ro" --bind "$auxiliary_cache:/sai-dependency-cache:ro"
  --bind "$support:$support:ro" --bind /var/tmp:/var/tmp:rw --pwd "$task"
  --env "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:?}"
  --env "SLURM_JOB_ID=$SLURM_JOB_ID" --env "TMPDIR=$task/tmp"
  --env SHELL=/bin/bash
  --env "GIT_CEILING_DIRECTORIES=$task" --env PYTHONDONTWRITEBYTECODE=1)
"${args[@]}" "$image" /bin/bash --noprofile --norc -c '
  set -eo pipefail
  prefix=$1; support=$2
  driver_ld=${LD_LIBRARY_PATH:-}
  source "$prefix/share/sai/runtime-env.sh"
  source /etc/profile.d/lmod.sh
  module purge
  module use /opt/modules/modulefiles/devtools
  module load cmake/3.31.6 openmpi/5.0.10-nvhpc26.3-gnu-cuda12-auto \
    fftw/3.3.10 saiblas/2603-gnu-auto cuda/12.9.1 nvmplibs/26.7-tmp gcc/14.3.0 spglib/2.7.0
  export LD_LIBRARY_PATH="$PWD/build/src:$LD_LIBRARY_PATH${driver_ld:+:$driver_ld}"
  export CP2K_DATA_DIR="$PWD/source/data" OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=1
  export OMPI_CC=/opt/devtools/gcc/14.3.0/bin/gcc OMPI_CXX=/opt/devtools/gcc/14.3.0/bin/g++
  export OMPI_FC=/opt/devtools/gcc/14.3.0/bin/gfortran
  exec /usr/bin/python3 "$support/build_and_test.py" "$prefix"
' cp2k-sai "$prefix" "$support"
