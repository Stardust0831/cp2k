"""Pin one V100 and four physical cores per rank; retain actual mapping."""
import csv
import json
import os
from pathlib import Path
import socket
import subprocess
import sys

root = Path(__file__).resolve().parent
rank = int(os.environ["OMPI_COMM_WORLD_RANK"])
local_rank = int(os.environ["OMPI_COMM_WORLD_LOCAL_RANK"])
visible = os.environ["CUDA_VISIBLE_DEVICES"]
rows = list(csv.reader(subprocess.check_output([
    "nvidia-smi", "-i", visible, "--query-gpu=index,uuid,name,memory.total",
    "--format=csv,noheader,nounits"], text=True).splitlines()))
assert len(rows) == int(os.environ["OMPI_COMM_WORLD_LOCAL_SIZE"])
gpu = rows[local_rank][1].strip()
cores = set()
for cpu in os.sched_getaffinity(0):
    p = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
    cores.add(tuple(int((p / key).read_text()) for key in ("physical_package_id", "core_id")))
assert len(cores) == 4, cores
mapping = dict(rank=rank, local_rank=local_rank, host=socket.gethostname(), gpu=gpu,
               cores=sorted(cores), job=os.environ["SLURM_JOB_ID"], arguments=sys.argv[1:])
Path(f"rank-{rank}.json").write_text(json.dumps(mapping, indent=2) + "\n")
monitor = stream = None
if local_rank == 0:
    stream = Path("gpu-samples.csv").open("w")
    monitor = subprocess.Popen([
        "nvidia-smi", "-i", visible,
        "--query-gpu=timestamp,index,uuid,memory.total,memory.used,utilization.gpu,power.draw",
        "--format=csv", "-l", "5"], stdout=stream, stderr=subprocess.STDOUT)
os.environ["CUDA_VISIBLE_DEVICES"] = gpu
scratch = Path("/var/tmp") / f"cp2k-cusolver-opt-{os.getuid()}-{os.environ['SLURM_JOB_ID']}-{rank}"
scratch.mkdir(exist_ok=True)
os.environ.update(APPTAINER_TMPDIR=str(scratch), APPTAINER_CACHEDIR=str(scratch / "cache"))
try:
    code = subprocess.call([sys.executable, str(root / "container.py"),
                            os.environ.get("SAI_VALIDATION_MODE", "run"), *sys.argv[1:]])
finally:
    if monitor is not None:
        monitor.terminate()
        monitor.wait(timeout=10)
        stream.close()
raise SystemExit(code)
