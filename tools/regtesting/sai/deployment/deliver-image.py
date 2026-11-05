#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Run on the daily producer: deliver accepted images into the tester's home."""
import argparse
from datetime import datetime, timezone
import errno
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys


def write(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def deliver(source, destination):
    from delivery_layout import load_artifact
    from source_cache import checksum
    candidates = []
    for path in (source / "containers/software/cp2k/development").glob("*/16V100/*.json"):
        record = json.loads(path.read_text())
        if not path.is_symlink() and all(record.get(k) is True for k in ("build_verified", "verified", "published")):
            candidates.append((path.with_suffix(".sif").stat().st_mtime_ns, str(path), path))
    if not candidates:
        raise ValueError("no accepted development/16V100 image")
    path = max(candidates)[2]
    original = path.read_bytes()
    record = json.loads(original)
    target = destination / path.relative_to(source)
    producer = dict(manifest=str(path), sha256=hashlib.sha256(original).hexdigest(), record=record)
    if target.exists():
        previous = json.loads(target.read_text())
        if previous.get("producer") != producer:
            raise ValueError("accepted producer manifest changed: " + str(path))
        return dict(manifest=str(target), source_sha=record["identity"]["source_sha"], image_sha256=record["sha256"])

    # The sidecar becomes visible only after both the image and source are ready.
    image = path.with_suffix(".sif")
    if load_artifact(source, image, software="cp2k", target="16v100-avx2") != record:
        raise ValueError("producer manifest changed during delivery")
    source_repo = source / "cache/repositories/cp2k"
    repository = destination / "cache/repositories/cp2k"
    if not repository.exists():
        repository.parent.mkdir(parents=True, exist_ok=True)
        # Fetch only accepted revisions. New files inherit the delivery ACL;
        # no producer object permissions or external alternates are retained.
        subprocess.run(["git", "init", "--bare", "--shared=group", str(repository)], check=True,
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    revision = record["identity"]["source_sha"]
    subprocess.run(["git", "--git-dir=" + str(repository), "fetch", "--depth=1", "--no-tags",
                    "--no-auto-maintenance", "--update-shallow", str(source_repo),
                    revision + ":refs/heads/accepted/" + revision], check=True)
    subprocess.run(["git", "--git-dir=" + str(repository), "cat-file", "-e",
                    revision + ":tests/do_regtest.py"], check=True)
    if (repository / "objects/info/alternates").exists():
        raise ValueError("delivered source must not depend on external object storage")

    target.parent.mkdir(parents=True, exist_ok=True)
    delivered_image = target.with_suffix(".sif")
    if not delivered_image.exists():
        temporary = delivered_image.with_suffix(".partial")
        temporary.unlink(missing_ok=True)
        try:
            os.link(image, temporary)
        except OSError as exc:
            if exc.errno != errno.EXDEV:
                raise
            shutil.copyfile(image, temporary)
            temporary.chmod(0o444)
        if checksum(temporary) != record["sha256"]:
            raise ValueError("delivered image checksum differs")
        temporary.replace(delivered_image)
    elif checksum(delivered_image) != record["sha256"]:
        raise ValueError("existing delivered image checksum differs")
    if path.read_bytes() != original:
        raise ValueError("producer manifest changed during delivery")
    write(target, dict(record, artifact=str(delivered_image), producer=producer))
    return dict(manifest=str(target), source_sha=revision, image_sha256=record["sha256"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("source", "destination", "controller"):
        parser.add_argument("--" + name, type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.controller))
    with (args.destination / "sync.lock").open("a+") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            result = dict(status="OK", **deliver(args.source, args.destination))
        except Exception as exc:
            result = dict(status="FAILED", error=str(exc))
        result["checked_at"] = datetime.now(timezone.utc).isoformat()
        write(args.destination / "sync.json", result)
        print(json.dumps(result), flush=True)
        return 0 if result["status"] == "OK" else 1


if __name__ == "__main__":
    raise SystemExit(main())
