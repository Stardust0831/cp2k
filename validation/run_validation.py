"""Build and compare independent patches in one four-V100 allocation."""
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

root = Path(__file__).resolve().parent
results = root / "results"
results.mkdir(exist_ok=True)
request = json.loads((root / "request.json").read_text())
source = root / "source"
container = ["python3", str(root / "container.py")]
mpi = ["mpirun", "-np", "4", "--map-by", "ppr:4:node:pe=4", "--bind-to", "core", "--report-bindings"]


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def cache_value(name):
    match = re.search(r"^" + re.escape(name) + r":[^=]+=(.*)$",
                      (root / "build/CMakeCache.txt").read_text(), re.M)
    assert match, name
    return match[1]


def run(name, command, *, cwd=root, timeout=1800, env=None):
    record = dict(command=[str(x) for x in command], cwd=str(cwd),
                  started_at=datetime.now(timezone.utc).isoformat())
    (results / "current-stage.json").write_text(json.dumps(dict(stage=name, **record), indent=2) + "\n")
    print("START", name, flush=True)
    start = time.monotonic()
    with (results / (name + ".log")).open("w") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT,
                                   cwd=cwd, timeout=timeout, env=env)
    record.update(exit_code=completed.returncode, seconds=time.monotonic() - start,
                  finished_at=datetime.now(timezone.utc).isoformat())
    (results / (name + ".json")).write_text(json.dumps(record, indent=2) + "\n")
    print("END", name, completed.returncode, round(record["seconds"], 2), flush=True)
    assert completed.returncode == 0, (name, completed.returncode)


assert digest(Path(request["image"])) == request["image_sha256"]
assert request["variants"][0]["name"] == "baseline" and not request["variants"][0].get("patch")
for filename, command in (
    ("allocation.txt", ["scontrol", "show", "job", "-o", os.environ["SLURM_JOB_ID"]]),
    ("lscpu.txt", ["lscpu"]),
    ("gpu-topology.txt", ["nvidia-smi", "topo", "-m"]),
):
    (results / filename).write_text(subprocess.check_output(command, text=True))
for name, expected in json.loads((root / "fixtures.json").read_text())["files"].items():
    assert digest(root / name) == expected, name

# Only the explicitly provided patch files can differ between builds.
originals = {}
for variant in request["variants"]:
    if variant.get("patch"):
        patch = root / variant["patch"]
        assert digest(patch) == variant["patch_sha256"]
        for name in re.findall(r"^\+\+\+ b/(.+)$", patch.read_text(), re.M):
            assert ".." not in Path(name).parts and (source / name).is_file(), name
            originals[name] = (source / name).read_bytes()
(results / "baseline-source-hashes.json").write_text(json.dumps(
    {name: hashlib.sha256(data).hexdigest() for name, data in originals.items()}, indent=2) + "\n")
(root / "bin").mkdir(exist_ok=True)
wrapper = root / "bin/cp2k.psmp"
wrapper.write_text('#!/bin/bash\nexec python3 "' + str(root / "rank.py") + '" "$@"\n')
wrapper.chmod(0o755)
run("configure", container + ["build", "python3", str(root / "configure.py")])

for variant in request["variants"]:
    name = variant["name"]
    for relative, data in originals.items():
        if (source / relative).read_bytes() != data:
            (source / relative).write_bytes(data)
    if variant.get("patch"):
        run(name + "-patch", ["patch", "-p1", "--batch", "-i", str(root / variant["patch"])], cwd=source)
    assert not (source / ".git").exists(), "Use the staged source archive, not a mutable checkout"
    (source / "REVISION").write_text(variant["source_sha"] + "\n")
    run(name + "-build", container + ["build", "cmake", "--build", str(root / "build"),
                                      "--target", "cp2k-bin", "--parallel", "16"])
    run(name + "-ldd", container + ["inspect", "ldd", str(root / "build/bin/cp2k.psmp")])
    ldd = (results / (name + "-ldd.log")).read_text()
    assert str(root / "build/src/libcp2k.so") in ldd and "not found" not in ldd
    run(name + "-version", container + ["inspect", str(root / "build/bin/cp2k.psmp"), "--version"])

    artifact = root / "artifacts" / name
    artifact.mkdir(parents=True, exist_ok=False)
    files = [root / "build/bin/cp2k.psmp", (root / "build/src/libcp2k.so").resolve(), root / "build/CMakeCache.txt"]
    hashes = {str(path.relative_to(root)): digest(path) for path in files}
    hashes.update({"source/" + relative: digest(source / relative) for relative in originals})
    (artifact / "hashes.json").write_text(json.dumps(hashes, indent=2) + "\n")
    for path in files:
        shutil.copy2(path, artifact / path.name)

    run(name + "-driver-build", container + ["build", "mpicc", "-O2", "-std=c11", "-Wall", "-Wextra", "-Werror",
                                             "-I" + str(source / "src/start"), "-L" + str(root / "build/src"),
                                             "-Wl,-rpath," + str(root / "build/src"), str(root / "lifecycle.c"),
                                             "-lcp2k", "-lm", "-o", str(root / "bin/lifecycle")])
    run(name + "-driver-ldd", container + ["inspect", "ldd", str(root / "bin/lifecycle")])
    ldd = (results / (name + "-driver-ldd.log")).read_text()
    assert str(root / "build/src/libcp2k.so") in ldd and "not found" not in ldd
    if request.get("count_context_calls", True):
        run(name + "-counter-build", container + ["build", "gcc", "-O2", "-shared", "-fPIC", "-Wall", "-Wextra", "-Werror",
            "-I" + str(Path(cache_value("CP2K_CUSOLVER_MP_ROOT")) / "include"),
            "-I" + str(Path(cache_value("CUDA_cudart_LIBRARY")).parents[1] / "include"),
            "-I" + str(Path(cache_value("CP2K_NCCL_ROOT")) / "include"),
            str(root / "cusolver_count.c"), "-ldl", "-o", str(root / "bin/cusolver_count.so")])
    lifecycle = results / name / "lifecycle"
    lifecycle.mkdir(parents=True)
    for path in (root / "inputs").iterdir():
        if path.is_file():
            shutil.copy2(path, lifecycle / path.name)
    run(name + "-lifecycle", mpi + ["python3", str(root / "rank.py")], cwd=lifecycle,
        timeout=300,
        env=dict(os.environ, SAI_VALIDATION_MODE="lifecycle",
                 SAI_CUSOLVER_COUNT="1" if request.get("count_context_calls", True) else "0"))

    suite = results / name / "regtest"
    suite.mkdir(parents=True, exist_ok=False)
    run(name + "-regtest", ["python3", str(source / "tests/do_regtest.py"), str(root / "bin"), "psmp",
                            "--mpiranks", "4", "--ompthreads", "4", "--maxtasks", "16", "--num_gpus", "4",
                            "--mpiexec", " ".join(mpi).replace("-np 4", "-np {N}"),
                            "--skip_unittests", "--testdirs", str(root / "TEST_DIRS"),
                            "--workbasedir", str(suite), "--cp2kdatadir", str(source / "data"), "--timeout", "180"])

    for family in ("ag-gamma", "diamond-gamma"):
        work = results / name / family
        shutil.copytree(root / "inputs" / family, work)
        run(name + "-" + family, mpi + ["python3", str(root / "rank.py"), "-i", "input.inp"], cwd=work)
        shutil.copy2(results / (name + "-" + family + ".log"), work / "cp2k.out")

run("audit", ["python3", str(root / "audit.py")])
(results / "completed.json").write_text(json.dumps(dict(
    variants=[v["name"] for v in request["variants"]], job=os.environ["SLURM_JOB_ID"],
    source_base=request["source_base"], finished_at=datetime.now(timezone.utc).isoformat()), indent=2) + "\n")
print("CP2K_CUSOLVER_OPT_VALIDATION_COMPLETE", flush=True)
