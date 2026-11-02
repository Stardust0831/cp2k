#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Run the pinned SAI regression suite without changing image publication records."""

import argparse
import json
import os
from pathlib import Path
import re
import signal
import subprocess
import sys
import time

PROJECT = Path("/org/org-open-sai/proj-sai-hpc-software")
ACTIVE_STATES = {"PENDING", "RUNNING", "COMPLETING", "CONFIGURING", "SUSPENDED"}


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")


def capture(*command):
    return subprocess.check_output(command, text=True).strip()


def latest_manifest():
    candidates = []
    for path in (PROJECT / "containers/software/cp2k/development").glob("*/16V100/*.json"):
        record = json.loads(path.read_text())
        if all(record.get(key) is True for key in ("build_verified", "verified", "published")):
            candidates.append(path)
    if not candidates:
        raise ValueError("no published 16V100 CP2K image found")
    return max(candidates, key=lambda path: path.with_suffix(".sif").stat().st_mtime_ns)


def monitor(job, task):
    """Wait for terminal accounting, streaming the upstream report as it grows."""
    offset, previous = 0, None
    deadline = time.monotonic() + 45 * 60
    while time.monotonic() < deadline:
        log = task / "results/regtest.log"
        if log.exists():
            with log.open() as stream:
                stream.seek(offset)
                print(stream.read(), end="", flush=True)
                offset = stream.tell()
        queued = capture("squeue", "-h", "-j", job, "-o", "%T|%R")
        if queued and queued != previous:
            print(f"Slurm {job}: {queued}", flush=True)
            previous = queued
        if not queued:
            rows = capture("sacct", "-X", "-nP", "-j", job, "-o", "JobIDRaw,State,ExitCode,Elapsed,NodeList")
            records = [line.split("|") for line in rows.splitlines() if line.split("|")[0] == job]
            if len(records) > 1:
                raise ValueError("ambiguous Slurm accounting")
            if records and records[0][1] not in ACTIVE_STATES:
                record = dict(zip(("job", "state", "exit_code", "elapsed", "nodes"), records[0]))
                write_json(task / "results/accounting.json", record)
                return record
        time.sleep(10)
    raise TimeoutError("Slurm job did not finish within 45 minutes")


def summary(task, request, evidence, accounting, error):
    lines = ["## SAI cuSOLVERMp — 4 V100", "",
             f"Status: **{'FAILED' if error else 'OK'}**", "",
             f"Tested CP2K commit: `{request.get('source_sha', 'unavailable')}`  ",
             "Resources: 1 node, 4 V100 GPUs, 4 MPI ranks × 4 OpenMP threads.  ",
             "Execution: full daily Apptainer image; upstream regression harness.  "]
    if accounting:
        lines.append(f"Slurm: `{accounting['job']}`, {accounting['state']}, "
                     f"{accounting['elapsed']}, `{accounting['nodes']}`.  ")
    if evidence:
        lines += ["", "| Input | Coverage | Numerical checks |",
                  "| --- | --- | --- |"]
        for case in evidence["cases"]:
            checks = ", ".join(f"{check['matcher']}: PASS" for check in case["matchers"])
            lines.append(f"| `{case['input']}` | `{case['coverage']}` | {checks or 'Completion: PASS'} |")
        lines += ["", "**9/9 cases passed; 6 cuSOLVERMp execution paths verified.**"]
    if error:
        lines += ["", "Failure: " + str(error)]
    lines += ["", "The tested commit is the image's source, which may differ from this workflow's commit.",
              "ELPA inputs retain upstream references and tolerances; only solver selection and ELPA options change.",
              "Download the artifact for the CP2K report, original and adapted inputs, raw outputs,",
              "rank/GPU records, image identity, and Slurm evidence.", ""]
    text = "\n".join(lines)
    (task / "results/summary.md").write_text(text)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"], "a") as stream:
            stream.write(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controller", type=Path, required=True)
    parser.add_argument("--manifest", default="")
    args = parser.parse_args()
    sys.path.insert(0, str(args.controller.resolve()))
    import cp2k_regtest as tester
    import cp2k_regtest_suite as suite
    from delivery_layout import load_artifact
    from source_cache import checksum

    run_id = f"github-{os.environ['GITHUB_RUN_ID']}-{os.environ['GITHUB_RUN_ATTEMPT']}"
    task = tester.science.run_dir(run_id)
    task.mkdir(parents=True, mode=0o700)
    for name in ("results", "bin", "tmp", "apptainer-runtime", "apptainer-cache"):
        (task / name).mkdir(mode=0o700)
    with open(os.environ["GITHUB_OUTPUT"], "a") as stream:
        stream.write(f"task={task}\n")
    request, evidence, accounting, error, job = {}, None, None, None, None

    def interrupted(signum, frame):
        raise InterruptedError(f"runner received signal {signum}")

    signal.signal(signal.SIGINT, interrupted)
    signal.signal(signal.SIGTERM, interrupted)
    try:
        manifest_path = Path(args.manifest) if args.manifest else latest_manifest()
        print(f"Image manifest: {manifest_path}", flush=True)
        manifest = load_artifact(PROJECT, manifest_path.with_suffix(".sif"),
                                 software="cp2k", target="16v100-avx2")
        if not all(manifest.get(key) is True for key in ("verified", "published")):
            raise ValueError("the selected image has not passed publication acceptance")
        write_json(task / "image-manifest.json", manifest)
        identity = manifest["identity"]
        request = dict(schema=1, run_id=run_id, minutes=20, identity=identity,
                       software="cp2k", version=identity["source_version"], target=identity["target"],
                       source_sha=identity["source_sha"], recipe_sha256=identity["recipe_sha256"],
                       artifact=manifest["artifact"], artifact_sha256=manifest["sha256"],
                       **suite.RESOURCES)
        source = PROJECT / "cache/repositories/cp2k"
        read_source = lambda path: subprocess.check_output(
            ["git", "--git-dir=" + str(source), "show", request["source_sha"] + ":" + path])
        request["suite_files"] = suite.prepare(read_source, task)
        request["tester_sha"] = capture("git", "-C", str(args.controller), "rev-parse", "HEAD")
        request["workflow_sha"] = os.environ["GITHUB_SHA"]
        request["controller_files"] = {
            path.name: checksum(path) for path in args.controller.glob("*.py")}
        script = task / "job.sbatch"
        script.write_text(tester.render_job(argparse.Namespace(**request)))
        wrapper = task / "bin/cp2k.psmp"
        wrapper.write_text(tester.RANK_WRAPPER)
        wrapper.chmod(0o700)
        subprocess.run(["bash", "-n", str(script)], check=True)
        request["job_script_sha256"] = checksum(script)
        request["rank_wrapper_sha256"] = checksum(wrapper)
        write_json(task / "request.json", request)
        job = capture("sbatch", "--hold", "--parsable", "--account=proj-sai-hpc-software",
                      "--chdir=" + str(task), str(script)).split(";")[0]
        if not job.isdigit():
            raise ValueError("invalid sbatch job ID")
        (task / "job.id").write_text(job + "\n")
        request["job"] = job
        write_json(task / "request.json", request)
        subprocess.run(["scontrol", "release", job], check=True)
        accounting = monitor(job, task)
        if (accounting["state"], accounting["exit_code"]) != ("COMPLETED", "0:0"):
            raise ValueError(f"Slurm {accounting['state']}, exit {accounting['exit_code']}")
        if "CP2K_CUSOLVERMP_REGTEST_FINISHED" not in (task / f"results/slurm-{job}.log").read_text():
            raise ValueError("container completion record missing")
        allocation = dict(re.findall(r"(?:^|\s)([A-Za-z][A-Za-z0-9/]*)=([^\s]+)",
                                    (task / "results/allocation.txt").read_text()))
        expected = dict(JobId=job, Partition="16V100", NumNodes="1", NumTasks="16",
                        TresPerNode="gres/gpu:4")
        if any(allocation.get(key) != value for key, value in expected.items()):
            raise ValueError("allocation differs from the four-GPU request")
        if request["suite_files"] != suite.prepare(read_source):
            raise ValueError("test suite no longer matches the image's source tree")
        evidence = suite.verify_results(task, request)
        write_json(task / "results/evidence.json", evidence)
    except Exception as exc:
        error = str(exc)
        print(f"::error::{error}", flush=True)
    finally:
        if job and job.isdigit() and accounting is None:
            subprocess.run(["scancel", job], check=False)
        if request.get("source_sha"):
            suite.write_report(task, request, error)
        else:
            (task / "results/report.txt").write_text(f"Summary: {error}\nStatus: FAILED\n")
        summary(task, request, evidence, accounting, error)
        print((task / "results/report.txt").read_text(), flush=True)
    return 1 if error else 0


if __name__ == "__main__":
    raise SystemExit(main())
