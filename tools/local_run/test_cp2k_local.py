#!/usr/bin/env python3
"""Small Linux-only regression tests; no MPI installation is needed."""

import multiprocessing
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from cp2k_local import Pool, cpu_list, topology


def compete(directory, allowed, start, results, release):
    pool = Pool(directory)
    try:
        start.wait()
        try:
            reserved, free, chosen = pool.reserve(
                [(0, 4), (1, 5), (2, 6)], allowed, 2, 1
            )
        except ValueError:
            results.put(None)
        else:
            results.put(chosen)
            release.wait(10)
    finally:
        pool.close()


class LocalRunTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.directory = Path(self.temp.name) / "pool"
        self.addCleanup(self.temp.cleanup)

    def pool(self):
        pool = Pool(self.directory)
        self.addCleanup(pool.close)
        return pool

    def test_sysfs_with_offline_and_nonadjacent_siblings(self):
        root = Path(self.temp.name)
        (root / "online").write_text("0-2,4-5\n")
        for cpu, siblings in [
            (0, "0,4"),
            (4, "0,4"),
            (1, "1,5"),
            (5, "1,5"),
            (2, "2,6"),
        ]:
            p = root / f"cpu{cpu}/topology"
            p.mkdir(parents=True)
            (p / "thread_siblings_list").write_text(siblings)
        online, groups = topology(root)
        self.assertEqual(online, {0, 1, 2, 4, 5})
        self.assertEqual(groups, [(0, 4), (1, 5), (2, 6)])
        _, _, selected = self.pool().reserve(groups, online & {4, 5, 6}, 1, 1)
        self.assertEqual(selected, [5])
        (root / "cpu4/topology/thread_siblings_list").write_text("4")
        with self.assertRaisesRegex(ValueError, "Inconsistent"):
            topology(root)

    def test_atomic_competing_sibling_views(self):
        ctx = multiprocessing.get_context("spawn")
        start, release, results = ctx.Event(), ctx.Event(), ctx.Queue()
        workers = [
            ctx.Process(
                target=compete, args=(self.directory, allowed, start, results, release)
            )
            for allowed in ({0, 1, 2}, {4, 5, 6})
        ]
        try:
            for worker in workers:
                worker.start()
            start.set()
            choices = [results.get(timeout=10) for _ in workers]
            self.assertEqual(sum(c is not None for c in choices), 1)
        finally:
            release.set()
            for worker in workers:
                worker.join(10)
                if worker.is_alive():
                    worker.kill()
                    worker.join()
        self.assertEqual([p.exitcode for p in workers], [0, 0])

    def test_shared_system_reserve_and_no_sibling_overlap(self):
        groups = [(i, i + 8) for i in range(8)]
        a = self.pool()
        system, _, selected_a = a.reserve(groups, set(range(8)), 3, 2)
        _, _, selected_b = self.pool().reserve(groups, set(range(8, 16)), 3, 0)
        self.assertEqual(system, groups[:2])
        self.assertEqual(selected_a, [2, 3, 4])
        self.assertEqual(selected_b, [13, 14, 15])
        with self.assertRaisesRegex(ValueError, "only 0 available"):
            self.pool().reserve(groups, set(range(16)), 1, 2)
        a.close()
        _, _, selected_c = self.pool().reserve(groups, set(range(16)), 3, 2)
        self.assertEqual(selected_c, [2, 3, 4])

    def test_larger_reserve_cannot_displace_live_calculation(self):
        groups = [(i,) for i in range(4)]
        self.pool().reserve(groups, set(range(4)), 1, 0)
        with self.assertRaisesRegex(ValueError, "System reserve"):
            self.pool().reserve(groups, set(range(4)), 1, 1)

    def test_failed_request_releases_partial_allocation(self):
        groups = [(i,) for i in range(4)]
        with self.assertRaisesRegex(ValueError, "Need 4 physical cores"):
            self.pool().reserve(groups, set(range(4)), 4, 1)
        self.assertEqual(self.pool().reserve(groups, set(range(4)), 3, 1)[2], [1, 2, 3])

    def test_private_directory_and_symlink_checks(self):
        self.directory.mkdir(mode=0o755)
        self.directory.chmod(0o755)  # Independent of the invoking shell's umask.
        with self.assertRaisesRegex(ValueError, "0700"):
            Pool(self.directory)
        self.directory.chmod(0o700)
        (self.directory / "allocation.lock").symlink_to(Path(self.temp.name) / "victim")
        with self.assertRaises(OSError):
            self.pool().reserve([(0,)], {0}, 1, 0)
        self.assertFalse((Path(self.temp.name) / "victim").exists())

    def test_owner_check(self):
        # No chown or second OS user is needed to exercise the UID rejection.
        from unittest.mock import patch

        self.directory.mkdir(mode=0o700)
        with patch("cp2k_local.os.getuid", return_value=os.getuid() + 1):
            with self.assertRaisesRegex(ValueError, "owned by this user"):
                Pool(self.directory)

    def test_cpu_list(self):
        self.assertEqual(cpu_list("0-2,7,9-10\n"), {0, 1, 2, 7, 9, 10})
        with self.assertRaises(ValueError):
            cpu_list("4-2")

    def test_orphan_in_separate_session_is_reaped_before_release(self):
        marker = Path(self.temp.name) / "orphan.pid"
        child = (
            "import os, signal, time; from pathlib import Path; "
            "os.setsid(); signal.signal(signal.SIGTERM, signal.SIG_IGN); "
            f"Path({str(marker)!r}).write_text(str(os.getpid())); time.sleep(60)"
        )
        application = (
            "import subprocess, sys, time; from pathlib import Path\n"
            f"subprocess.Popen([sys.executable, '-c', {child!r}])\n"
            f"while not Path({str(marker)!r}).exists(): time.sleep(.01)\n"
            "sys.exit(17)\n"
        )
        supervisor = (
            "import os, sys; from cp2k_local import Pool, run_guarded\n"
            f"pool = Pool({str(self.directory)!r})\n"
            "pool.reserve([(0,)], {0}, 1, 0)\n"
            f"sys.exit(run_guarded([sys.executable, '-c', {application!r}], "
            "list(os.sched_getaffinity(0)), pool))\n"
        )
        completed = subprocess.run(
            [sys.executable, "-c", supervisor],
            cwd=Path(__file__).resolve().parent,
            timeout=15,
        )
        self.assertEqual(completed.returncode, 17)
        self.assertFalse(Path(f"/proc/{marker.read_text()}").exists())
        self.assertEqual(self.pool().reserve([(0,)], {0}, 1, 0)[2], [0])


if __name__ == "__main__":
    unittest.main()
