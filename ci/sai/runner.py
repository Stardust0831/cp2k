#!/usr/bin/env python3
"""Bridge cp2k-ci's SSH target protocol to the SAI Slurm service account."""

import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import shutil
import signal
import subprocess
import sys
import tarfile
import tempfile
import time

HERE = Path(__file__).resolve().parent
REPOSITORY = "https://github.com/cp2k/cp2k.git"


def check_request(branch, revision):
    if not re.fullmatch(r"master|pull/[1-9][0-9]*/merge", branch):
        raise ValueError("Expected master or pull/NUMBER/merge")
    if not re.fullmatch(r"[0-9a-f]{40}", revision):
        raise ValueError("Expected a full commit SHA")


def request_arguments(arguments, original_command):
    receiving = arguments[:1] == ["--receive"]
    arguments = arguments[1:] if receiving else arguments
    if not arguments and original_command:
        command = shlex.split(original_command)
        expected = (
            [str(HERE / "runner.py"), "--receive"]
            if receiving
            else ["ci/sai/runner.py"]
        )
        if command[: len(expected)] != expected:
            raise ValueError("Only the configured CI command is allowed")
        arguments = command[len(expected) :]
    if len(arguments) != 2:
        raise ValueError("Usage: runner.py [--receive] BRANCH COMMIT_SHA")
    check_request(*arguments)
    return receiving, arguments


def capture(*args):
    return subprocess.check_output(args, text=True).strip()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n")


def final_status(text):
    statuses = re.findall(r"^Status: (\S+)$", text, re.M)
    return statuses[-1] if statuses else "UNKNOWN"


def relay(settings, branch, revision):
    print(f"CommitSHA: {revision}", flush=True)
    root = Path(settings["work_root"]).resolve()
    root.mkdir(parents=True, exist_ok=True)
    cache = root / "source.git"
    with (root / "source.lock").open("w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if not cache.exists():
            subprocess.run(["git", "init", "--quiet", "--bare", str(cache)], check=True)
        exists = (
            subprocess.run(
                ["git", "-C", str(cache), "cat-file", "-e", revision + "^{commit}"],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            ).returncode
            == 0
        )
        if not exists:
            subprocess.run(
                ["git", "-C", str(cache), "fetch", "--depth=1", REPOSITORY, revision],
                check=True,
                timeout=300,
            )
        if (
            capture("git", "-C", str(cache), "rev-parse", revision + "^{commit}")
            != revision
        ):
            raise ValueError("Fetched revision does not match the request")
    # Fetch by immutable SHA: the PR merge ref may advance while a job is queued.
    command = shlex.join([settings["remote_command"], "--receive", branch, revision])
    log = root / (revision + "-" + str(time.time_ns()) + ".log")
    with subprocess.Popen(
        [
            "ssh",
            "-o",
            "BatchMode=yes",
            "-o",
            "ConnectTimeout=30",
            "-o",
            "ServerAliveInterval=30",
            "-o",
            "ServerAliveCountMax=3",
            settings["ssh_host"],
            command,
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
    ) as remote:
        first = remote.stdout.readline()
        print(first, end="", flush=True)
        if first == "Source: upload\n":
            with tempfile.TemporaryFile(dir=root) as archive:
                archive_source(cache, revision, archive)
                print(f"Source archive: {archive.tell()} bytes", flush=True)
                archive.seek(0)
                shutil.copyfileobj(archive, remote.stdin.buffer)
        remote.stdin.close()
        with log.open("w") as output:
            output.write(first)
            for line in remote.stdout:
                output.write(line)
                output.flush()
                print(line, end="", flush=True)
        code = remote.wait()
    if code or final_status(log.read_text()) != "OK":
        raise RuntimeError(f"SAI run failed; relay log: {log}")


def archive_source(repository, revision, output):
    # The configured site cache can be shared with the daily build account.
    subprocess.run(
        [
            "git",
            "-c",
            f"safe.directory={repository}",
            "-C",
            str(repository),
            "archive",
            "--format=tar.gz",
            revision,
            "CMakeLists.txt",
            "CMakePresets.json",
            "LICENSE",
            "cmake",
            "src",
            "data",
            "tests",
            "tools/build_utils",
        ],
        stdout=output,
        check=True,
    )


def extract_source(stream, destination, revision):
    with tarfile.open(fileobj=stream, mode="r|*") as archive:
        if archive.pax_headers.get("comment") != revision:
            raise ValueError("Git archive does not identify the requested commit")
        archive.extractall(destination, filter="data")
    if not (destination / "CMakeLists.txt").is_file():
        raise ValueError("Incomplete CP2K source archive")
    (destination / "REVISION").write_text(revision + "\n")


def select_image(catalog):
    candidates = []
    for path in Path(catalog).glob("*/16V100/*.json"):
        item = json.loads(path.read_text())
        if all(
            item.get(k) is True for k in ("build_verified", "verified", "published")
        ):
            candidates.append((path.stat().st_mtime_ns, path, item))
    if not candidates:
        raise ValueError("No accepted SAI dependency image")
    _, path, item = max(candidates, key=lambda x: x[0])
    image = Path(item["artifact"])
    with image.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    if digest != item["sha256"]:
        raise ValueError("Dependency image checksum mismatch")
    return dict(
        image=str(image),
        prefix=item["identity"]["install_prefix"],
        sha256=digest,
        manifest=str(path),
        build_tools=item["controller"],
    )


def interrupted(signum, frame):
    raise InterruptedError(f"Interrupted by signal {signum}")


def receive(settings, branch, revision):
    root = Path(settings["work_root"]).resolve()
    root.mkdir(parents=True, exist_ok=True)
    task = Path(tempfile.mkdtemp(prefix=revision[:12] + "-", dir=root))
    source = task / "source"
    source.mkdir()
    repository = settings.get("source_repository")
    cached = (
        repository
        and subprocess.run(
            [
                "git",
                "-c",
                f"safe.directory={repository}",
                "-C",
                repository,
                "cat-file",
                "-e",
                revision + "^{commit}",
            ],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode
        == 0
    )
    print("Source: cache" if cached else "Source: upload", flush=True)
    if cached:
        with tempfile.TemporaryFile(dir=task) as archive:
            archive_source(repository, revision, archive)
            archive.seek(0)
            extract_source(archive, source, revision)
    else:
        extract_source(sys.stdin.buffer, source, revision)
    request = dict(
        branch=branch,
        revision=revision,
        environment=select_image(settings["image_catalog"]),
    )
    request["environment"]["auxiliary_cache"] = settings["auxiliary_cache"]
    write_json(task / "request.json", request)
    print(f"SAIRun: {task.name}\nRawFiles: {task}", flush=True)
    for signum in (signal.SIGHUP, signal.SIGINT, signal.SIGTERM):
        signal.signal(signum, interrupted)
    job, finished = None, False
    try:
        job = capture(
            "sbatch",
            "--hold",
            "--parsable",
            "--account=org-open-sai",
            "--partition=16V100",
            "--qos=rush-gpu",
            "--nodes=1",
            "--ntasks=16",
            "--gpus-per-node=4",
            "--time=00:30:00",
            "--export=NIL",
            "--job-name=cp2k-sai-ci",
            "--chdir=" + str(task),
            "--output=" + str(task / "slurm.log"),
            str(HERE / "job.sh"),
            str(task),
            str(HERE),
        ).split(";")[0]
        if not job.isdigit():
            job = None
            raise ValueError("Invalid sbatch response")
        (task / "job.id").write_text(job + "\n")
        subprocess.run(["scontrol", "release", job], check=True)
        print(f"SlurmJob: {job}", flush=True)
        deadline = time.monotonic() + 7200
        offset = 0
        while True:
            queued = capture("squeue", "-h", "-j", job, "-o", "%T")
            if (task / "slurm.log").exists():
                with (task / "slurm.log").open() as stream:
                    stream.seek(offset)
                    for line in stream:
                        # Only the aggregate result below is a dashboard status.
                        if line.startswith(("Status:", "Summary:")):
                            line = "Regtest| " + line
                        print(line, end="", flush=True)
                    offset = stream.tell()
            if not queued:
                break
            if time.monotonic() >= deadline:
                raise TimeoutError("SAI queue/run deadline exceeded")
            time.sleep(10)
        accounting = capture(
            "sacct",
            "-X",
            "-nP",
            "-j",
            job,
            "-o",
            "JobIDRaw,State,ExitCode,ElapsedRaw,AllocTRES,NodeList",
        )
        (task / "accounting.txt").write_text(accounting + "\n")
        rows = [
            line.split("|")
            for line in accounting.splitlines()
            if line.split("|")[0] == job
        ]
        if len(rows) != 1 or rows[0][1:3] != ["COMPLETED", "0:0"]:
            raise RuntimeError("Slurm did not complete successfully: " + accounting)
        finished = True
        result = json.loads((task / "result.json").read_text())
        if result["status"] != "OK" or result["revision"] != revision:
            raise RuntimeError("Build/test validation failed")
        print(f"Summary: {result['summary']}\nStatus: OK", flush=True)
    finally:
        if job and not finished:
            subprocess.run(["scancel", job], check=False)


def main():
    try:
        receiving, args = request_arguments(
            sys.argv[1:], os.environ.get("SSH_ORIGINAL_COMMAND", "")
        )
        branch, revision = args
        settings = json.loads((HERE / "config.json").read_text())
        (receive if receiving else relay)(settings, branch, revision)
    except Exception as exc:
        print(f"Summary: SAI CI failed: {exc}\nStatus: FAILED", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
