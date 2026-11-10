#!/usr/bin/env python3
"""Build the requested CP2K revision and run unchanged GPU regression inputs."""

import csv
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import time
import tomllib

SUITES = ("regtest-cusolver", "regtest-elpa", "regtest-elpa-2", "regtest-elpa-qr")


def run(command, log):
    print("Running:", shlex.join(map(str, command)), flush=True)
    start = time.monotonic()
    try:
        with Path(log).open("w") as output:
            subprocess.run(command, stdout=output, stderr=subprocess.STDOUT, check=True)
    except subprocess.CalledProcessError:
        print(
            "\n".join(Path(log).read_text(errors="replace").splitlines()[-100:]),
            flush=True,
        )
        raise
    print(f"Finished {log}: {time.monotonic() - start:.2f} s", flush=True)


def cache_options(path):
    result = {}
    for line in path.read_text().splitlines():
        if (
            not line
            or line.startswith(("#", "//"))
            or ":" not in line
            or "=" not in line
        ):
            continue
        key, value = line.split("=", 1)
        name, kind = key.split(":", 1)
        if kind != "INTERNAL":
            result[name] = (kind, value)
    return result


def restore_dependencies(prefix):
    # Runtime images omit headers/static archives. Reuse the matching accepted
    # build cache and its publisher's checksum/relocation implementation.
    sys.path.insert(0, "/sai-build-tools")
    from cp2k_auxiliary_cache import unpack

    state = json.loads((prefix / "share/sai/auxiliary-cache.json").read_text())
    evidence = prefix / "dependencies/.restore-evidence"
    evidence.mkdir()
    record = unpack(
        Path("/sai-dependency-cache") / state["key"],
        prefix,
        evidence,
        state["contract"],
    )
    if record["archive_sha256"] != state["cache_archive_sha256"]:
        raise RuntimeError("Build dependency cache does not match the image")
    receipt = {key: record[key] for key in ("key", "archive_sha256", "components")}
    Path("dependency-cache.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(
        "Restored accepted build dependencies:",
        ", ".join(record["components"]),
        flush=True,
    )


def configure(source, build, prefix):
    cache = cache_options(prefix / "share/sai/CMakeCache.txt")
    names = {
        "CMAKE_BUILD_TYPE",
        "CMAKE_C_COMPILER",
        "CMAKE_CXX_COMPILER",
        "CMAKE_Fortran_COMPILER",
        "CMAKE_C_FLAGS",
        "CMAKE_CXX_FLAGS",
        "CMAKE_Fortran_FLAGS",
        "CMAKE_CUDA_ARCHITECTURES",
        "CMAKE_PREFIX_PATH",
        "CMAKE_LIBRARY_PATH",
        "DBCSR_DIR",
        "Libint2_DIR",
        "libxsmm_DIR",
        "HDF5_ROOT",
        "HDF5_PREFER_PARALLEL",
        "HDF5_USE_STATIC_LIBRARIES",
        "ZLIB_ROOT",
        "BUILD_SHARED_LIBS",
        "CP2K_WITH_GPU",
        "CP2K_CUSOLVERMP_USE_NCCL",
        "CP2K_SCALAPACK_VENDOR",
        "CP2K_NVPL_SCALAPACK_MPI",
    }
    names.update(
        k for k in cache if k.startswith(("CP2K_USE_", "CP2K_ENABLE_", "CP2K_BLAS_"))
    )
    names.update(k for k in cache if k.startswith("CP2K_") and k.endswith("_ROOT"))
    names.update(
        k
        for k, (_, v) in cache.items()
        if k.endswith("_DIR") and v.startswith((str(prefix) + "/", "/opt/devtools/"))
    )
    args = ["cmake", "-S", str(source), "-B", str(build)]
    for name in sorted(names & cache.keys()):
        kind, value = cache[name]
        args.append(f"-D{name}:{'STRING' if kind == 'UNINITIALIZED' else kind}={value}")
    prerequisites = build / "prerequisites.cmake"
    prerequisites.write_text("find_package(ZLIB REQUIRED)\n")
    args += [
        f"-DCMAKE_PROJECT_cp2k_INCLUDE={prerequisites}",
        f"-DFYPP_EXECUTABLE={source / 'tools/build_utils/fypp'}",
    ]
    receipt = json.loads((prefix / "share/sai/spack-dependencies.json").read_text())
    if not receipt["dependencies_installed"]:
        raise RuntimeError("Incomplete dependency snapshot")
    for name in ("libint", "libxsmm", "hdf5", "plumed"):
        os.environ["CP2K_" + name.upper() + "_ROOT"] = receipt["roots"][name]
    run(args, "configure.log")
    current = cache_options(build / "CMakeCache.txt")
    enabled = lambda values: {
        k
        for k, (_, v) in values.items()
        if k.startswith(("CP2K_USE_", "CP2K_ENABLE_")) and v in ("ON", "CUDA")
    }
    missing = enabled(cache) - enabled(current)
    if missing:
        raise RuntimeError("Lost dependency features: " + ", ".join(sorted(missing)))


def rank(binary, arguments):
    local_rank = int(os.environ.get("OMPI_COMM_WORLD_LOCAL_RANK", "0"))
    devices = json.loads(os.environ["SAI_GPU_DEVICES"])
    os.environ["CUDA_VISIBLE_DEVICES"] = devices[local_rank]
    if "OMPI_COMM_WORLD_RANK" in os.environ:
        rank_id = int(os.environ["OMPI_COMM_WORLD_RANK"])
        name = Path(arguments[-1]).name
        record = dict(
            rank=rank_id,
            size=int(os.environ["OMPI_COMM_WORLD_SIZE"]),
            gpu=devices[local_rank],
            cpus=sorted(os.sched_getaffinity(0)),
        )
        Path(f"rank-{name}-{rank_id}.json").write_text(json.dumps(record) + "\n")
    os.execv(binary, [binary, *arguments])


def flags(text):
    return set(re.search(r"^\s*cp2kflags:\s*(.*)$", text, re.M)[1].split())


def main(prefix):
    task = Path.cwd()
    source, build = task / "source", task / "build"
    request = json.loads((task / "request.json").read_text())
    prefix = Path(prefix)
    restore_dependencies(prefix)
    run([prefix / "bin/cp2k.psmp", "--version"], "dependency-version.txt")
    configure(source, build, prefix)
    run(
        ["cmake", "--build", str(build), "--target", "cp2k-bin", "--parallel", "16"],
        "build.log",
    )
    binary = build / "bin/cp2k.psmp"
    run(["ldd", binary], "ldd.txt")
    linkage = (task / "ldd.txt").read_text()
    if "not found" in linkage or f"=> {build}/src/libcp2k" not in linkage:
        raise RuntimeError(
            "The test must load the newly built libcp2k and all dependencies"
        )
    run([binary, "--version"], "version.txt")
    version = (task / "version.txt").read_text()
    if request["revision"] not in version:
        raise RuntimeError("Built CP2K revision does not match the requested SHA")
    required = flags((task / "dependency-version.txt").read_text())
    required |= {"cusolvermp", "cusolvermp_nccl", "elpa", "dbcsr_acc", "offload_cuda"}
    if required - flags(version):
        raise RuntimeError(
            "Built CP2K lost enabled features: " + str(required - flags(version))
        )
    devices = [
        row[0].strip()
        for row in csv.reader(
            subprocess.check_output(
                [
                    "nvidia-smi",
                    "-i",
                    os.environ["CUDA_VISIBLE_DEVICES"],
                    "--query-gpu=uuid",
                    "--format=csv,noheader",
                ],
                text=True,
            ).splitlines()
        )
    ]
    if len(devices) != 4 or len(set(devices)) != 4:
        raise RuntimeError("Expected four distinct GPUs")
    os.environ["SAI_GPU_DEVICES"] = json.dumps(devices)
    testbin = build / "testbin"
    testbin.mkdir()
    wrapper = testbin / "cp2k.psmp"
    wrapper.write_text(
        "#!/bin/sh\nexec "
        + shlex.join(
            [sys.executable, str(Path(__file__).resolve()), "--rank", str(binary)]
        )
        + ' "$@"\n'
    )
    wrapper.chmod(0o755)
    args = [
        sys.executable,
        str(source / "tests/do_regtest.py"),
        str(testbin),
        "psmp",
        "--mpiranks",
        "4",
        "--ompthreads",
        "4",
        "--maxtasks",
        "16",
        "--skip_unittests",
        "--timeout",
        "300",
        "--workbasedir",
        str(task / "regtesting"),
        "--mpiexec",
        "mpirun -n {N} --map-by ppr:4:node:pe=4 --bind-to core",
    ]
    for suite in SUITES:
        args += ["--restrictdir", "QS/" + suite]
    run(args, "regtest.log")
    report = (task / "regtest.log").read_text()
    print(report, end="", flush=True)
    if re.findall(r"^Status: (\S+)$", report, re.M) != ["OK"]:
        raise RuntimeError("Regression harness did not report OK")
    workdirs = list((task / "regtesting").glob("TEST-*"))
    if len(workdirs) != 1:
        raise RuntimeError("Expected one regression run")
    timers, cases = set(), 0
    for suite in SUITES:
        definitions = tomllib.loads(
            (source / "tests/QS" / suite / "TEST_FILES.toml").read_text()
        )
        for name in definitions:
            directory = workdirs[0] / "QS" / suite
            output = (directory / (name + ".out")).read_text()
            if "PROGRAM ENDED AT" not in output:
                raise RuntimeError("Incomplete case: " + name)
            timers.update(
                re.findall(
                    r"^\s*(cp_(?:c?fm_general|fm_diag)_cusolver|cp_fm_diag_elpa)\s+[1-9][0-9]*\s",
                    output,
                    re.M,
                )
            )
            assigned = []
            for rank_id in range(4):
                record = json.loads(
                    (directory / f"rank-{name}-{rank_id}.json").read_text()
                )
                if record["rank"] != rank_id or record["size"] != 4:
                    raise RuntimeError("Incorrect MPI mapping")
                assigned.append(record["gpu"])
            if set(assigned) != set(devices):
                raise RuntimeError("Incorrect GPU mapping")
            cases += 1
    expected = {"cp_fm_general_cusolver", "cp_cfm_general_cusolver", "cp_fm_diag_elpa"}
    if expected - timers:
        raise RuntimeError("Missing solver execution: " + str(expected - timers))
    summary = (
        f"{cases} cases passed on 4 V100; real/complex cuSOLVERMp and ELPA verified"
    )
    result = dict(
        status="OK",
        revision=request["revision"],
        summary=summary,
        cases=cases,
        timers=sorted(timers),
        flags=sorted(flags(version)),
        devices=devices,
    )
    (task / "result.json").write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    if sys.argv[1] == "--rank":
        rank(sys.argv[2], sys.argv[3:])
    else:
        main(sys.argv[1])
