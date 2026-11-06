"""Run the build and tests in the same isolated container and dependency mounts."""
import json
import os
from pathlib import Path
import sys

root = Path(__file__).resolve().parent
request = json.loads((root / "request.json").read_text())
args = ["apptainer", "exec", "--nv", "--cleanenv", "--no-home", "--no-mount",
        "bind-paths,home,cwd,tmp,hostfs"]
for path in ("/usr", "/lib", "/lib64", "/etc", "/opt/devtools", "/opt/modules"):
    args += ["--bind", path + ":" + path + ":ro"]
args += ["--bind", str(root) + ":" + str(root) + ":rw",
         "--bind", request["dependency_directory"] + ":" + request["dependency_prefix"] + ":ro",
         "--bind", "/var/tmp:/var/tmp:rw", "--pwd", str(Path.cwd())]
for key, value in os.environ.items():
    if key.startswith(("SLURM_", "OMPI_", "OPAL_", "PMIX_", "PMI_", "PRTE_", "UCX_", "NCCL_", "CUDA_", "FI_", "OMP_")) or key in ("TMPDIR", "OPENBLAS_NUM_THREADS", "SAI_CUSOLVER_COUNT"):
        args += ["--env", key + "=" + value]
args += [request["image"], "/usr/bin/bash", "--noprofile", "--norc",
         str(root / "inside.sh"), *sys.argv[1:]]
os.execvp(args[0], args)
