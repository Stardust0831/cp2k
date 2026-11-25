#!/usr/bin/env python3
"""Cooperatively reserve physical cores for a local Open MPI / OpenMP run."""

import argparse
import ctypes
import fcntl
import json
import os
from pathlib import Path
import re
import select
import shlex
import shutil
import signal
import stat
import subprocess
import sys
import time


def cpu_list(text):
    cpus = set()
    for part in text.strip().split(","):
        bounds = [int(i) for i in part.split("-")]
        if len(bounds) == 1:
            bounds *= 2
        if len(bounds) != 2 or not 0 <= bounds[0] <= bounds[1]:
            raise ValueError(f"Invalid CPU list: {text!r}")
        cpus.update(range(bounds[0], bounds[1] + 1))
    return cpus


def topology(root=Path("/sys/devices/system/cpu")):
    """Use complete sibling groups, including siblings outside our affinity."""
    online = cpu_list((root / "online").read_text())
    groups = {}
    for cpu in sorted(online):
        siblings = tuple(
            sorted(
                cpu_list((root / f"cpu{cpu}/topology/thread_siblings_list").read_text())
            )
        )
        if cpu not in siblings:
            raise ValueError(f"CPU {cpu} is missing from its sysfs sibling group")
        groups[cpu] = siblings
    for cpu, siblings in groups.items():
        if any(groups[s] != siblings for s in set(siblings) & online):
            raise ValueError(f"Inconsistent sysfs sibling groups at CPU {cpu}")
    return online, sorted(set(groups.values()))


class Pool:
    """Persistent lock inodes; closing descriptors releases a reservation.

    EX locks reserve calculation cores; SH locks keep system cores unused by
    calculations. The short global mutex makes a whole request atomic. Never
    unlink lock files: another process may still have the old inode open.
    """

    def __init__(self, directory):
        directory = Path(directory)
        try:
            directory.mkdir(mode=0o700)
        except FileExistsError:
            pass
        self.directory = os.open(
            directory, os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW
        )
        self.fds = []
        info = os.fstat(self.directory)
        if info.st_uid != os.getuid() or stat.S_IMODE(info.st_mode) != 0o700:
            self.close()
            raise ValueError(
                f"Pool must be owned by this user with mode 0700: {directory}"
            )

    def open_lock(self, name):
        fd = os.open(
            name,
            os.O_CREAT | os.O_RDWR | os.O_NOFOLLOW | os.O_CLOEXEC,
            0o600,
            dir_fd=self.directory,
        )
        info = os.fstat(fd)
        if (
            not stat.S_ISREG(info.st_mode)
            or info.st_uid != os.getuid()
            or info.st_nlink != 1
            or stat.S_IMODE(info.st_mode) != 0o600
        ):
            os.close(fd)
            raise ValueError(f"Unsafe pool lock: {name}")
        return fd

    def close(self):
        for fd in self.fds:
            os.close(fd)
        self.fds.clear()
        if self.directory is not None:
            os.close(self.directory)
            self.directory = None

    def reserve(self, groups, allowed, count, spare):
        if spare >= len(groups):
            raise ValueError(
                f"Cannot leave {spare} of {len(groups)} physical cores unused"
            )
        mutex = self.open_lock("allocation.lock")
        free = []
        reserved = groups[:spare]
        try:
            fcntl.flock(mutex, fcntl.LOCK_EX)
            for group in reserved + [g for g in groups[spare:] if set(g) & allowed]:
                fd = self.open_lock("core-" + "-".join(map(str, group)) + ".lock")
                system = group in reserved
                try:
                    fcntl.flock(
                        fd, (fcntl.LOCK_SH if system else fcntl.LOCK_EX) | fcntl.LOCK_NB
                    )
                except BlockingIOError:
                    os.close(fd)
                    if system:
                        raise ValueError(
                            f"System reserve {list(group)} is in use; retry after that run ends"
                        )
                    continue
                self.fds.append(fd)
                if not system:
                    free.append(group)
            if len(free) < count:
                raise ValueError(
                    f"Need {count} physical cores, only {len(free)} available "
                    f"after {spare} system cores and other participants' reservations"
                )
            # Keep the first count calculation locks and all system locks.
            for fd in self.fds[spare + count :]:
                os.close(fd)
            del self.fds[spare + count :]
            return reserved, free, [min(set(g) & allowed) for g in free[:count]]
        except BaseException:
            for fd in self.fds:
                os.close(fd)
            self.fds.clear()
            raise
        finally:
            os.close(mutex)


def rank_exec(layout, command):
    """Set OS CPU numbers before loading the application's OpenMP runtime."""
    rank = int(os.environ["OMPI_COMM_WORLD_RANK"])
    if (
        int(os.environ["OMPI_COMM_WORLD_SIZE"]) != len(layout)
        or int(os.environ["OMPI_COMM_WORLD_LOCAL_SIZE"]) != len(layout)
        or int(os.environ["OMPI_COMM_WORLD_LOCAL_RANK"]) != rank
    ):
        raise ValueError("Expected one local Open MPI application")
    cpus = layout[rank]
    if not set(cpus) <= os.sched_getaffinity(0):
        raise ValueError(f"Rank {rank}: selected CPUs are outside inherited affinity")
    os.sched_setaffinity(0, cpus)
    env = os.environ.copy()
    for key in ("GOMP_CPU_AFFINITY", "KMP_AFFINITY", "KMP_HW_SUBSET"):
        env.pop(key, None)
    env.update(
        OMP_NUM_THREADS=str(len(cpus)),
        OMP_THREAD_LIMIT=str(len(cpus)),
        OMP_MAX_ACTIVE_LEVELS="1",
        OMP_DYNAMIC="FALSE",
        OMP_PLACES=",".join("{%d}" % cpu for cpu in cpus),
        OMP_PROC_BIND="close",
        OPENBLAS_NUM_THREADS="1",
        MKL_NUM_THREADS="1",
        MKL_DYNAMIC="FALSE",
        BLIS_NUM_THREADS="1",
    )
    print(
        f"cp2k-local: rank {rank} CPUs {cpus}, OMP_PLACES={env['OMP_PLACES']}",
        file=sys.stderr,
        flush=True,
    )
    os.execvpe(command[0], command, env)


def signal_children(signum):
    """Orphaned descendants become direct children of this subreaper."""
    children = Path(f"/proc/self/task/{os.getpid()}/children").read_text()
    for pid in map(int, children.split()):
        try:
            fd = os.pidfd_open(pid)
            try:
                # Check parentage *after* opening the pidfd. A recycled PID must
                # never cause us to signal another user's or unrelated process.
                fields = Path(f"/proc/{pid}/stat").read_text().rsplit(")", 1)[1].split()
                if int(fields[1]) == os.getpid():
                    signal.pidfd_send_signal(fd, signum)
            finally:
                os.close(fd)
        except (ProcessLookupError, FileNotFoundError):
            pass


def supervise(command, cpus, parent_pipe, cancelled):
    # Adopt orphaned MPI children so reservations outlive *all* descendants.
    libc = ctypes.CDLL(None, use_errno=True)
    if libc.prctl(36, 1, 0, 0, 0) != 0:  # PR_SET_CHILD_SUBREAPER
        raise OSError(ctypes.get_errno(), "Cannot become a child subreaper")
    proc = None
    result = None
    stopping = None
    try:
        proc = subprocess.Popen(
            command,
            start_new_session=True,
            preexec_fn=lambda: os.sched_setaffinity(0, cpus),
        )
        while True:
            try:
                while True:
                    pid, status = os.waitpid(-1, os.WNOHANG)
                    if pid == 0:
                        break
                    if pid == proc.pid:
                        proc.returncode = os.waitstatus_to_exitcode(status)
                        result = (
                            proc.returncode
                            if proc.returncode >= 0
                            else 128 - proc.returncode
                        )
            except ChildProcessError:
                return 128 + cancelled[0] if cancelled[0] else result
            if not cancelled[0] and select.select([parent_pipe], [], [], 0)[0]:
                cancelled[0] = signal.SIGTERM  # launcher died, including SIGKILL
            if cancelled[0] or result is not None:
                if stopping is None:
                    stopping = time.monotonic()
                signal_children(
                    signal.SIGKILL
                    if time.monotonic() - stopping > 2
                    else signal.SIGTERM
                )
            time.sleep(0.05)
    finally:
        # Exceptions must not release CPU locks while MPI children still run.
        while True:
            signal_children(signal.SIGKILL)
            try:
                pid, status = os.waitpid(-1, os.WNOHANG)
            except ChildProcessError:
                break
            if proc is not None and pid == proc.pid:
                proc.returncode = os.waitstatus_to_exitcode(status)
            if pid == 0:
                time.sleep(0.05)


def run_guarded(command, cpus, pool):
    """An independent supervisor keeps locks even if the launcher is killed."""
    signals = {signal.SIGINT, signal.SIGTERM, signal.SIGHUP}
    oldmask = signal.pthread_sigmask(signal.SIG_BLOCK, signals)
    reader, writer = os.pipe2(os.O_CLOEXEC)
    sys.stdout.flush()
    sys.stderr.flush()
    pid = os.fork()
    if pid == 0:
        os.close(writer)
        os.setsid()
        cancelled = [0]

        def cancel(signum, frame):
            cancelled[0] = signum

        for signum in signals:
            signal.signal(signum, cancel)
        signal.pthread_sigmask(signal.SIG_SETMASK, oldmask)
        code = 2
        try:
            code = supervise(command, cpus, reader, cancelled)
        except BaseException as exc:
            print(f"cp2k-local: supervisor: {exc}", file=sys.stderr, flush=True)
        finally:
            pool.close()
            os.close(reader)
            os._exit(code if code is not None else 2)

    os.close(reader)
    pool.close()  # close, never LOCK_UN: the supervisor shares these descriptions

    def forward(signum, frame):
        try:
            os.kill(pid, signum)
        except ProcessLookupError:
            pass

    previous = {s: signal.signal(s, forward) for s in signals}
    signal.pthread_sigmask(signal.SIG_SETMASK, oldmask)
    try:
        _, status = os.waitpid(pid, 0)
        code = os.waitstatus_to_exitcode(status)
        return code if code >= 0 else 128 - code
    finally:
        os.close(writer)
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-n", "--ranks", type=int, required=True)
    parser.add_argument("-t", "--threads", type=int, required=True)
    parser.add_argument(
        "--reserve-cores",
        type=int,
        default=1,
        help="whole system cores left unused by participating calculations (default: 1)",
    )
    parser.add_argument("--mpirun", default="mpirun", help="Open MPI 4/5 executable")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print a plan and release it without launching",
    )
    parser.add_argument(
        "command",
        nargs=argparse.REMAINDER,
        help="-- cp2k.psmp -i input.inp -o output.out",
    )
    args = parser.parse_args(argv)
    if args.command[:1] == ["--"]:
        args.command.pop(0)
    if args.ranks < 1 or args.threads < 1 or args.reserve_cores < 0 or not args.command:
        parser.error(
            "positive ranks/threads, nonnegative reserve-cores and a command are required"
        )
    if (
        sys.platform != "linux"
        or not hasattr(os, "pidfd_open")
        or not hasattr(signal, "pidfd_send_signal")
    ):
        parser.error(
            "Linux >= 5.3 and Python >= 3.9 built with pidfd support are required"
        )
    if os.getuid() == 0 or os.getuid() != os.geteuid():
        parser.error("run as an ordinary, non-setuid user")
    if any(
        k in os.environ
        for k in ("SLURM_JOB_ID", "PBS_JOBID", "LSB_JOBID", "OMPI_COMM_WORLD_RANK")
    ):
        parser.error(
            "this launcher is for a local workstation, outside a scheduler or MPI job"
        )
    mpirun = shutil.which(args.mpirun)
    if not mpirun or not shutil.which(args.command[0]):
        parser.error("mpirun and the application must be executable")
    version = subprocess.run(
        [mpirun, "--version"], text=True, capture_output=True, timeout=10, check=True
    ).stdout
    if not re.search(r"\(Open MPI\) [45]\.", version):
        parser.error("only Open MPI 4 and 5 are supported")
    online, groups = topology()
    allowed = online & os.sched_getaffinity(0)
    pool = Pool(Path(f"/tmp/cp2k-local-{os.getuid()}"))
    try:
        reserved, free, selected = pool.reserve(
            groups, allowed, args.ranks * args.threads, args.reserve_cores
        )
        layout = [
            selected[i : i + args.threads]
            for i in range(0, len(selected), args.threads)
        ]
        command = [
            mpirun,
            "--host",
            f"localhost:{args.ranks}",
            "-n",
            str(args.ranks),
            "--map-by",
            "slot",
            "--bind-to",
            "none",
            sys.executable,
            str(Path(__file__).resolve()),
            "--_rank",
            json.dumps(layout),
            *args.command,
        ]
        plan = dict(
            online=sorted(online),
            allowed=sorted(allowed),
            physical_cores=groups,
            system_cores=reserved,
            available_cores=free,
            available_cpus=sorted(allowed & {c for g in free for c in g}),
            selected_cpus=selected,
            rank_cpus=layout,
            dry_run=args.dry_run,
        )
        print("cp2k-local: plan " + json.dumps(plan), file=sys.stderr, flush=True)
        print("cp2k-local: command " + shlex.join(command), file=sys.stderr, flush=True)
        if args.dry_run:
            return 0
        return run_guarded(command, selected, pool)
    finally:
        pool.close()


if __name__ == "__main__":
    try:
        if sys.argv[1:2] == ["--_rank"]:
            rank_exec(json.loads(sys.argv[2]), sys.argv[3:])
        else:
            sys.exit(main())
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as exc:
        print(f"cp2k-local: {exc}", file=sys.stderr)
        sys.exit(2)
