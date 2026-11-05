#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-2.0-or-later
"""Drive the restricted SAI endpoint from the VPS GitHub runner."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile
import time


def remote(*args, binary=False):
    command = ["ssh", "-F", os.environ["SAI_SSH_CONFIG"], "sai-cp2k-tester", " ".join(args)]
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=180)
    if result.returncode:
        raise RuntimeError(result.stderr.decode(errors="replace").strip())
    return result.stdout if binary else json.loads(result.stdout)


def outputs(**items):
    if os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a") as stream:
            for key, value in items.items():
                stream.write(f"{key}={value}\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    destination = Path("artifacts")
    destination.mkdir(exist_ok=True)
    run_id = "github-" + os.environ["GITHUB_RUN_ID"] + "-" + os.environ["GITHUB_RUN_ATTEMPT"]
    state = None
    latest = None
    try:
        inspect = remote("inspect")
        if inspect["uid"] != 1478400824 or not inspect["container"] or inspect["home_ssh_visible"] or inspect["code_writable"]:
            raise ValueError("SAI submission container isolation check failed")
        (destination / "submission-container.json").write_text(json.dumps(inspect, indent=2) + "\n")
        latest = remote("latest")
        print(f"Latest accepted CP2K image: {latest['source_sha']} ({latest['image_sha256']})", flush=True)
        state = remote("start", run_id, latest["image_sha256"], os.environ["GITHUB_SHA"],
                       "force" if args.force else "once")
        if state.get("skipped"):
            print(f"Image already tested by {state['run_id']}: {state['status']}; no GPU job submitted.")
            outputs(publish="false", skipped="true")
            if os.environ.get("GITHUB_STEP_SUMMARY"):
                Path(os.environ["GITHUB_STEP_SUMMARY"]).write_text(
                    f"No new image. Previous test: [{state['status']}]({state['action_url']}).\n")
            return 0
        run_id = state["run_id"]
        deadline = time.monotonic() + 45 * 60
        previous = None
        while not state.get("terminal"):
            if time.monotonic() > deadline:
                raise TimeoutError(f"Slurm run {run_id} is still active; a later Action will resume it")
            text = f"Slurm {state.get('job')}: {state.get('state')} {state.get('detail', '')}"
            if text != previous:
                print(text, flush=True)
                previous = text
            time.sleep(15)
            # A lost observation does not resubmit/cancel the live allocation.
            for attempt in range(3):
                try:
                    state = remote("status", run_id)
                    break
                except (RuntimeError, subprocess.TimeoutExpired):
                    if attempt == 2:
                        raise
                    time.sleep(10)
        raw = destination / "raw-files.tar.gz"
        raw.write_bytes(remote("fetch", run_id, binary=True))
        with tarfile.open(raw) as archive:
            for source, target in (("results/report.txt", "report.txt"),
                                   ("results/summary.md", "summary.md"),
                                   ("results/final.json", "latest.json")):
                member = archive.getmember(source)
                if not member.isfile() or member.size > 3_000_000:
                    raise ValueError("invalid report artifact")
                (destination / target).write_bytes(archive.extractfile(member).read())
        final = json.loads((destination / "latest.json").read_text())
        report = (destination / "report.txt").read_text()
        if final != state or re.findall(r"^Status: (\S+)$", report, re.M)[-1:] != [state["status"]]:
            raise ValueError("downloaded report differs from final server state")
        print(report, flush=True)
        outputs(publish="true", skipped="false")
        if os.environ.get("GITHUB_STEP_SUMMARY"):
            Path(os.environ["GITHUB_STEP_SUMMARY"]).write_text((destination / "summary.md").read_text())
        return 0 if state["status"] == "OK" else 1
    except Exception as exc:
        print(f"::error::{exc}", flush=True)
        # Infrastructure failures replace stale success, without claiming that
        # an untested image passed or that a still-running allocation failed.
        source = (state or latest or {}).get("source_sha")
        report = f"CommitSHA: {source}\n" if source else ""
        report += f"Summary: Tester infrastructure: {' '.join(str(exc).split())}\nStatus: UNKNOWN\n"
        (destination / "report.txt").write_text(report)
        final = dict(status="UNKNOWN", source_sha=source, run_id=run_id,
                     completed_at=datetime.now(timezone.utc).isoformat(), error=str(exc),
                     action_url=f"{os.environ['GITHUB_SERVER_URL']}/{os.environ['GITHUB_REPOSITORY']}"
                                f"/actions/runs/{os.environ['GITHUB_RUN_ID']}")
        (destination / "latest.json").write_text(json.dumps(final, indent=2) + "\n")
        (destination / "summary.md").write_text("SAI tester infrastructure error: " + str(exc) + "\n")
        outputs(publish="true", skipped="false")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
