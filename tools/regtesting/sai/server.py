#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Trusted SAI operations: fixed suite, resources, paths, and job template."""
import argparse
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import tarfile

ROOT = Path("/org/org-open-sai/proj-sai-hpc-software/cp2k-cusolvermp-tester")
PROJECT = ROOT.parent
STATE = ROOT / "state"
CONTROL = ROOT / "code/controller"
ACTIVE = {"PENDING", "RUNNING", "COMPLETING", "CONFIGURING", "SUSPENDED", "REQUEUED"}
sys.path.insert(0, str(CONTROL))


def write(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def capture(*args):
    return subprocess.check_output(args, text=True, stderr=subprocess.PIPE).strip()


def records():
    result = []
    for path in (PROJECT / "containers/software/cp2k/development").glob("*/16V100/*.json"):
        item = json.loads(path.read_text())
        if all(item.get(k) is True for k in ("build_verified", "verified", "published")):
            result.append((path.with_suffix(".sif").stat().st_mtime_ns, path, item))
    if not result:
        raise ValueError("no accepted development/16V100 image")
    return sorted(result, key=lambda item: (item[0], str(item[1])))


def task_dir(run_id):
    if not re.fullmatch(r"github-[0-9]{1,20}-[0-9]{1,6}", run_id):
        raise ValueError("invalid run ID")
    task = STATE / "runtime-tests" / run_id
    if task.resolve() != task or task.is_symlink():
        raise ValueError("redirected task directory")
    return task


def latest():
    _, path, item = records()[-1]
    return dict(source_sha=item["identity"]["source_sha"], image_sha256=item["sha256"],
                manifest=str(path), server_version=json.loads((ROOT / "code/version.json").read_text()))


def source_reader(sha):
    return lambda path: subprocess.check_output([
        "git", "-c", "safe.directory=" + str(PROJECT / "cache/repositories/cp2k"),
        "--git-dir=" + str(PROJECT / "cache/repositories/cp2k"), "show", sha + ":" + path])


def describe(task):
    request = json.loads((task / "request.json").read_text())
    record = dict(run_id=task.name, source_sha=request["source_sha"],
                  image_sha256=request["artifact_sha256"])
    if (task / "results/final.json").exists():
        return {**record, **json.loads((task / "results/final.json").read_text())}
    if not (task / "job.id").exists():
        if not (task / "submission.started").exists():
            return finalize(task, request, initial_error="Setup interrupted before Slurm submission")
        # sbatch is held until its ID is durable. Recover a lost SSH response
        # by the fixed job name and comment; never guess or submit another job.
        queued = capture("squeue", "-h", "-u", "stardust", "-o", "%i|%j|%k")
        matches = [row.split("|")[0] for row in queued.splitlines()
                   if row.split("|")[1:] == ["cp2k-cusolvermp-" + task.name, task.name]]
        if len(matches) != 1 or not matches[0].isdigit():
            raise ValueError("Unresolved Slurm submission for " + task.name + "; inspect before retrying")
        (task / "job.id").write_text(matches[0] + "\n")
    job = (task / "job.id").read_text().strip()
    if not job.isdigit():
        raise ValueError("invalid recorded Slurm job")
    if request.get("job") != job:
        request["job"] = job
        write(task / "request.json", request)
    queued = capture("squeue", "-h", "-j", job, "-o", "%T|%R")
    if queued:
        if queued == "PENDING|JobHeldUser" and not (task / "job.released").exists():
            release(task, job)
            return dict(record, job=job, state="RELEASED", terminal=False)
        return dict(record, job=job, state=queued.split("|", 1)[0], detail=queued, terminal=False)
    rows = capture("sacct", "-X", "-nP", "-j", job, "-o",
                   "JobIDRaw,State,ExitCode,Elapsed,NodeList,User,AllocTRES")
    rows = [line.split("|") for line in rows.splitlines() if line.split("|")[0] == job]
    if not rows or rows[0][1] in ACTIVE:
        return dict(record, job=job, state="ACCOUNTING_PENDING", terminal=False)
    if len(rows) != 1:
        raise ValueError("ambiguous Slurm accounting")
    accounting = dict(zip(("job", "state", "exit_code", "elapsed", "nodes", "user", "allocation"), rows[0]))
    write(task / "results/accounting.json", accounting)
    return finalize(task, request, accounting)


def release(task, job):
    subprocess.run(["scontrol", "release", job], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    (task / "job.released").write_text(job + "\n")


def finalize(task, request, accounting=None, initial_error=None):
    import cp2k_regtest_suite as suite
    error, evidence = initial_error, None
    if error is None:
        try:
            if (accounting["state"], accounting["exit_code"], accounting["user"]) != ("COMPLETED", "0:0", "stardust"):
                raise ValueError(f"Slurm {accounting['state']}, exit {accounting['exit_code']}")
            if "CP2K_CUSOLVERMP_REGTEST_FINISHED" not in (task / f"results/slurm-{accounting['job']}.log").read_text():
                raise ValueError("container completion marker missing")
            allocation = dict(re.findall(r"(?:^|\s)([A-Za-z][A-Za-z0-9/]*)=([^\s]+)",
                                         (task / "results/allocation.txt").read_text()))
            expected = dict(JobId=accounting["job"], Partition="16V100", NumNodes="1",
                            NumTasks="16", TresPerNode="gres/gpu:4")
            if any(allocation.get(k) != v for k, v in expected.items()):
                raise ValueError("allocation differs from the four-GPU request")
            for name, key in (("job.sbatch", "job_script_sha256"), ("bin/cp2k.psmp", "rank_wrapper_sha256")):
                if hashlib.sha256((task / name).read_bytes()).hexdigest() != request[key]:
                    raise ValueError("submitted job or rank wrapper changed")
            if request["suite_files"] != suite.prepare(source_reader(request["source_sha"])):
                raise ValueError("test fixtures differ from the image source")
            for name, digest in request["controller_files"].items():
                if hashlib.sha256((CONTROL / name).read_bytes()).hexdigest() != digest:
                    raise ValueError("trusted controller changed")
            evidence = suite.verify_results(task, request)
            write(task / "results/evidence.json", evidence)
        except Exception as exc:
            error = str(exc)
    suite.write_report(task, request, error)
    now = datetime.now(timezone.utc).isoformat()
    with (task / "results/report.txt").open("a") as stream:
        stream.write(f"ReportTime: {now}\nActionURL: {request['action_url']}\n")
    final = dict(run_id=task.name, source_sha=request["source_sha"],
                 image_sha256=request["artifact_sha256"], terminal=True,
                 status="FAILED" if error else "OK", error=error, completed_at=now,
                 action_url=request["action_url"], accounting=accounting,
                 job=(accounting or {}).get("job"), server_version=request["server_version"])
    lines = ["## SAI cuSOLVERMp — 4 V100", "", f"Status: **{final['status']}**", "",
             f"CP2K: `{request['source_sha']}`", "",
             "Full daily SIF; 4 MPI ranks × 4 OpenMP threads; no explicit memory request.", ""]
    if accounting:
        lines += [f"Slurm `{accounting['job']}`: {accounting['state']}, "
                  f"{accounting['elapsed']}, `{accounting['nodes']}`.", ""]
    if evidence:
        lines += [f"{len(evidence['cases'])}/9 cases passed; {evidence['numerical_checks']} numerical checks; "
                  "6 cuSOLVERMp execution paths verified.", "",
                  "| Input | Backend coverage |", "| --- | --- |"]
        lines += [f"| `{case['input']}` | `{case['coverage']}` |" for case in evidence["cases"]]
    if error:
        lines += ["Failure: " + str(error)]
    lines += ["", f"[Action]({request['action_url']})", "",
              "Raw artifacts include original/adapted inputs, outputs, backend timers, GPU/rank mapping, "
              "image and source identity, and final Slurm accounting.", ""]
    (task / "results/summary.md").write_text("\n".join(lines))
    write(task / "results/final.json", final)
    (STATE / "completed").mkdir(exist_ok=True)
    write(STATE / "completed" / (request["artifact_sha256"] + ".json"), final)
    return final


def start(run_id, image_sha, workflow_sha, mode):
    import cp2k_regtest as tester
    import cp2k_regtest_suite as suite
    from delivery_layout import load_artifact
    from source_cache import checksum
    task = task_dir(run_id)
    if task.exists():
        return describe(task)
    # An interrupted SSH/Action connection never causes a second allocation.
    for other in sorted((STATE / "runtime-tests").glob("github-*")):
        if not (other / "results/final.json").exists():
            previous = describe(other)
            if not previous["terminal"]:
                return dict(previous, resumed=True)
    version = json.loads((ROOT / "code/version.json").read_text())
    previous_path = STATE / "completed" / (image_sha + ".json")
    if mode == "once" and previous_path.exists():
        previous = json.loads(previous_path.read_text())
        if previous["server_version"] == version:
            return dict(previous, skipped=True)
    choices = [item for _, _, item in records() if item["sha256"] == image_sha]
    if not choices:
        raise ValueError("requested image is not an accepted 16V100 daily image")
    manifest = load_artifact(PROJECT, choices[-1]["artifact"], software="cp2k", target="16v100-avx2")
    identity = manifest["identity"]
    task.mkdir(mode=0o700)
    for name in ("results", "bin", "tmp", "apptainer-runtime", "apptainer-cache"):
        (task / name).mkdir(mode=0o700)
    request = dict(schema=1, run_id=run_id, minutes=20, identity=identity, software="cp2k",
                   version=identity["source_version"], target=identity["target"],
                   source_sha=identity["source_sha"], recipe_sha256=identity["recipe_sha256"],
                   artifact=manifest["artifact"], artifact_sha256=manifest["sha256"],
                   workflow_sha=workflow_sha, server_version=version,
                   action_url="https://github.com/Stardust0831/cp2k/actions/runs/" + run_id.split("-")[1],
                   **suite.RESOURCES)
    write(task / "request.json", request)
    write(task / "image-manifest.json", manifest)
    try:
        request["suite_files"] = suite.prepare(source_reader(request["source_sha"]), task)
        request["controller_files"] = {p.name: checksum(p) for p in CONTROL.glob("*.py")}
        script = task / "job.sbatch"
        script.write_text(tester.render_job(argparse.Namespace(**request)))
        wrapper = task / "bin/cp2k.psmp"
        wrapper.write_text(tester.RANK_WRAPPER)
        wrapper.chmod(0o700)
        subprocess.run(["bash", "-n", str(script)], check=True)
        request["job_script_sha256"] = checksum(script)
        request["rank_wrapper_sha256"] = checksum(wrapper)
        write(task / "request.json", request)
        (task / "submission.started").write_text(datetime.now(timezone.utc).isoformat() + "\n")
        try:
            job = capture("sbatch", "--hold", "--parsable", "--account=proj-sai-hpc-software",
                          "--comment=" + run_id, "--chdir=" + str(task), str(script)).split(";")[0]
        except subprocess.CalledProcessError as exc:
            (task / "results/submission-error.txt").write_text(exc.stderr or str(exc))
            (task / "submission.started").unlink()
            raise
        if not job.isdigit():
            raise ValueError("invalid Slurm response; inspect held allocation")
        (task / "job.id").write_text(job + "\n")
        request["job"] = job
        write(task / "request.json", request)
        release(task, job)
    except Exception as exc:
        if (task / "submission.started").exists():
            raise  # Preserve the allocation and resume by its recorded ID.
        return finalize(task, request, initial_error=str(exc))
    return dict(run_id=run_id, job=job, source_sha=request["source_sha"],
                image_sha256=image_sha, state="SUBMITTED", terminal=False)


def fetch(run_id):
    task = task_dir(run_id)
    if not (task / "results/final.json").is_file():
        raise ValueError("run is not finalized")
    archive = STATE / "archives" / (run_id + ".tar.gz")
    archive.parent.mkdir(exist_ok=True)
    if not archive.exists():
        temporary = archive.with_suffix(".tmp")
        with tarfile.open(temporary, "w:gz") as handle:
            for path in sorted(task.rglob("*")):
                relative = path.relative_to(task)
                if relative.parts[0] in ("tmp", "apptainer-runtime", "apptainer-cache"):
                    continue
                if path.is_file() and not path.is_symlink() and path.resolve() == path:
                    handle.add(path, arcname=str(relative), recursive=False)
        temporary.replace(archive)
    with archive.open("rb") as source:
        import shutil
        shutil.copyfileobj(source, sys.stdout.buffer)


def main(argv):
    if argv == ["inspect"]:
        return dict(uid=os.getuid(), container=os.environ.get("APPTAINER_CONTAINER"),
                    home_ssh_visible=Path("/home/stardust/.ssh").exists(),
                    code_writable=os.access(ROOT / "code", os.W_OK),
                    slurm=capture("sbatch", "--version"), state_root=str(STATE),
                    server_version=json.loads((ROOT / "code/version.json").read_text()))
    with (STATE / "gateway.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if argv == ["latest"]:
            return latest()
        if len(argv) == 5 and argv[0] == "start":
            return start(*argv[1:])
        if len(argv) == 2 and argv[0] == "status":
            return describe(task_dir(argv[1]))
        if len(argv) == 2 and argv[0] == "fetch":
            fetch(argv[1])
            return None
    raise ValueError("invalid operation")


if __name__ == "__main__":
    try:
        result = main(sys.argv[1:])
        if result is not None:
            print(json.dumps(result), flush=True)
    except Exception as exc:
        print(f"SAI tester: {exc}", file=sys.stderr, flush=True)
        raise SystemExit(1)
