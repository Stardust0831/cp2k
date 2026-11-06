#!/usr/bin/env bash
#SBATCH --job-name=cp2k-cusolver-fixture-check
#SBATCH --partition=16V100
#SBATCH --qos=flood-gpu
#SBATCH --account=org-open-sai
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --gpus-per-node=4
#SBATCH --time=10
#SBATCH --export=NIL
#SBATCH --chdir=/org/org-open-sai/cp2k-ci/cp2k-cusolvermp-tester/state/development/cusolver-opt-20261106
#SBATCH --output=/org/org-open-sai/cp2k-ci/cp2k-cusolvermp-tester/state/development/cusolver-opt-20261106/diagnostics-full-spectrum/slurm-%j.log
set -eo pipefail
umask 077
ulimit -c 0
export USER=${SLURM_JOB_USER:?} LOGNAME=${SLURM_JOB_USER:?}
export PATH=/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin
export LD_LIBRARY_PATH="" LD_PRELOAD="" PYTHONDONTWRITEBYTECODE=1
unset APPTAINER_BIND APPTAINER_BINDPATH SINGULARITY_BIND SINGULARITY_BINDPATH
source /etc/profile.d/lmod.sh
module purge
module use /opt/modules/modulefiles/devtools /opt/modules/modulefiles/apps
module load apptainer/1.4.4 openmpi/5.0.10-nvhpc26.3-gnu-cuda12-auto
set -u
export SLURM_EXPORT_ENV=ALL
export OMPI_MCA_plm_slurm_args=--external-launcher PRTE_MCA_plm_slurm_args=--external-launcher
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 OMP_PLACES=cores OMP_PROC_BIND=close OMP_DYNAMIC=FALSE
export TMPDIR=/var/tmp
exec env HOME=/org/org-open-sai/cp2k-ci python3 diagnostics-full-spectrum/run.py
