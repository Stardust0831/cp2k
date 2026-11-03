#!/usr/bin/env python3
"""Extend the verified Ubuntu Base with the host's matching Ubuntu Python/git."""
import hashlib
import os
from pathlib import Path
import re
import shutil
import subprocess
import tarfile

ROOT = Path("/org/org-open-sai/proj-sai-hpc-software/cp2k-cusolvermp-tester")
BASE_SHA = "e77b6f10c2590cef872b33ee9f635a0e3fd1f57fb074c0e52b5c7f56147a0c86"


def main():
    archive = ROOT / "ubuntu-base-24.04.5-base-amd64.tar.gz"
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == BASE_SHA
    filesystem = ROOT / "rootfs"
    filesystem.mkdir(mode=0o700)
    with tarfile.open(archive) as handle:
        # Official base tarball is hash-pinned and contains normal Unix symlinks.
        handle.extractall(filesystem, filter="tar")
    for name in ("passwd", "group"):
        with (filesystem / "etc" / name).open("a") as stream:
            for account in ("stardust", "slurm", "munge"):
                stream.write(subprocess.check_output(["getent", name, account], text=True))
    shutil.copytree("/usr/lib/python3.12", filesystem / "usr/lib/python3.12",
                    ignore=shutil.ignore_patterns("__pycache__"))
    binaries = [Path("/usr/bin/python3.12"), Path("/usr/bin/git")]
    libraries = list(Path("/usr/lib/python3.12").rglob("*.so"))
    for source in binaries:
        shutil.copy2(source, filesystem / str(source).lstrip("/"))
    (filesystem / "usr/bin/python3").symlink_to("python3.12")
    for binary in binaries + libraries:
        output = subprocess.run(["ldd", str(binary)], capture_output=True, text=True, check=True).stdout
        for path in re.findall(r"(?:=>\s+|^\s*)(/[^\s]+)", output, re.M):
            source = Path(path)
            target = filesystem / str(source).lstrip("/")
            if not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source.resolve(), target)
    # Empty the copied writable caches; this image is executed read-only.
    environment = {key: value for key, value in os.environ.items()
                   if not key.startswith(("APPTAINER_", "SINGULARITY_"))}
    environment.update(APPTAINER_TMPDIR=str(ROOT / "scratch"),
                       APPTAINER_CACHEDIR=str(ROOT / "cache"), TMPDIR=str(ROOT / "scratch"))
    subprocess.run(["/opt/devtools/apptainer/1.4.4/bin/apptainer", "build",
                    "--mksquashfs-args", "-processors 2 -mem 256M",
                    str(ROOT / "ubuntu-submit.sif"), str(filesystem)], check=True, env=environment)
    digest = hashlib.sha256((ROOT / "ubuntu-submit.sif").read_bytes()).hexdigest()
    (ROOT / "ubuntu-submit.sif.sha256").write_text(digest + "  ubuntu-submit.sif\n")
    print(digest)


if __name__ == "__main__":
    main()
