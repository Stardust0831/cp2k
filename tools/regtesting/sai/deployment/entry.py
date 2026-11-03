#!/usr/bin/env python3
"""The dedicated SSH key always enters this fixed SAI submission container."""
import os
from pathlib import Path
import re
import shlex
import sys

ROOT = Path("/org/org-open-sai/proj-sai-hpc-software/cp2k-cusolvermp-tester")
PROJECT = ROOT.parent
APPTAINER = "/opt/devtools/apptainer/1.4.4/bin/apptainer"


def parse_command(command):
    words = shlex.split(command)
    if words == ["latest"] or words == ["inspect"]:
        return words
    if len(words) == 2 and words[0] in ("status", "fetch"):
        if re.fullmatch(r"github-[0-9]{1,20}-[0-9]{1,6}", words[1]):
            return words
    if len(words) == 5 and words[0] == "start":
        if (re.fullmatch(r"github-[0-9]{1,20}-[0-9]{1,6}", words[1])
                and re.fullmatch(r"[0-9a-f]{64}", words[2])
                and re.fullmatch(r"[0-9a-f]{40}", words[3])
                and words[4] in ("once", "force")):
            return words
    raise ValueError("allowed operations: inspect, latest, start, status, fetch")


def main():
    words = parse_command(os.environ.get("SSH_ORIGINAL_COMMAND", ""))
    # No client command is interpreted by a host shell. No SSH_ORIGINAL_COMMAND
    # or inherited loader/Python/Apptainer settings enter the container.
    environment = dict(PATH="/usr/bin:/bin", HOME="/home/stardust", USER="stardust",
                       LOGNAME="stardust", LANG="C.UTF-8",
                       APPTAINER_TMPDIR=str(ROOT / "scratch"),
                       APPTAINER_CACHEDIR=str(ROOT / "cache"))
    args = [APPTAINER, "exec", "--cleanenv", "--containall", "--no-home",
            "--no-mount", "bind-paths,home,cwd,tmp,hostfs"]
    paths = [(str(ROOT / "code"), "ro"), (str(ROOT / "state"), "rw"),
             (str(PROJECT / "containers/software/cp2k/development"), "ro"),
             (str(PROJECT / "cache/repositories/cp2k"), "ro"),
             ("/etc/slurm", "ro"), ("/etc/hosts", "ro"), ("/run/munge", "ro"),
             ("/usr/lib/x86_64-linux-gnu/slurm", "ro")]
    for path, mode in paths:
        args += ["--bind", f"{path}:{path}:{mode}"]
    for name in ("sbatch", "squeue", "sacct", "scontrol", "scancel"):
        args += ["--bind", f"/usr/bin/{name}:/usr/bin/{name}:ro"]
    for name in ("libmunge.so.2", "libreadline.so.8", "libtinfo.so.6"):
        path = Path("/lib/x86_64-linux-gnu") / name
        args += ["--bind", f"{path.resolve()}:/usr/lib/x86_64-linux-gnu/{name}:ro"]
    args += ["--env", "SLURM_CONF=/etc/slurm/slurm.conf",
             "--env", "LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu/slurm",
             "--env", "PYTHONDONTWRITEBYTECODE=1",
             "--env", f"SAI_SOFTWARE_ROOT={ROOT / 'state'}",
             "--env", "PATH=/usr/bin:/bin", "--pwd", str(ROOT / "state"),
             str(ROOT / "ubuntu-submit.sif"), "/usr/bin/python3",
             str(ROOT / "code/server.py"), *words]
    os.execve(APPTAINER, args, environment)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError) as exc:
        print(f"SAI tester: {exc}", file=sys.stderr)
        raise SystemExit(2)
