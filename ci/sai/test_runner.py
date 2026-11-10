#!/usr/bin/env python3
"""Check the SSH request, source archive, and report boundaries without Slurm."""

import io
import json
from pathlib import Path
import subprocess
import tarfile
import tempfile
import unittest

import runner


class RunnerTests(unittest.TestCase):
    revision = "a" * 40

    def test_master_and_pull_merge(self):
        for branch in ("master", "pull/6456/merge"):
            runner.check_request(branch, self.revision)

    def test_reject_other_refs_and_shell_text(self):
        for branch in ("main", "pull/2/head", "master; true", "../master", "-x"):
            with self.subTest(branch=branch), self.assertRaises(ValueError):
                runner.check_request(branch, self.revision)

    def test_require_full_immutable_sha(self):
        for revision in ("HEAD", "abc1234", "a" * 39, "a" * 41, "a" * 39 + ";"):
            with self.subTest(revision=revision), self.assertRaises(ValueError):
                runner.check_request("master", revision)

    def test_forced_ssh_command(self):
        receiving, args = runner.request_arguments(
            ["--receive"], f"{runner.HERE}/runner.py --receive master {self.revision}"
        )
        self.assertTrue(receiving)
        self.assertEqual(args, ["master", self.revision])
        with self.assertRaises(ValueError):
            runner.request_arguments(["--receive"], "sh -c whoami")
        with self.assertRaises(ValueError):
            runner.request_arguments(
                [], f"ci/sai/runner.py master {self.revision}; true"
            )

    def archive(self, revision, member="CMakeLists.txt", link=None):
        stream = io.BytesIO()
        with tarfile.open(
            fileobj=stream,
            mode="w:gz",
            format=tarfile.PAX_FORMAT,
            pax_headers={"comment": revision},
        ) as archive:
            entry = tarfile.TarInfo(member)
            if link:
                entry.type, entry.linkname = tarfile.SYMTYPE, link
            archive.addfile(entry)
        stream.seek(0)
        return stream

    def test_git_archive_revision_is_preserved(self):
        with tempfile.TemporaryDirectory() as tmp:
            runner.extract_source(self.archive(self.revision), Path(tmp), self.revision)
            self.assertEqual((Path(tmp) / "REVISION").read_text(), self.revision + "\n")

    def test_cached_archive_ignores_working_tree_changes(self):
        with tempfile.TemporaryDirectory() as tmp:
            repository, source = Path(tmp) / "repo", Path(tmp) / "source"
            repository.mkdir()
            git = ["git", "-C", str(repository)]
            subprocess.run([*git, "init", "--quiet"], check=True)
            for name in (
                "CMakeLists.txt",
                "CMakePresets.json",
                "LICENSE",
                "cmake/test",
                "src/test",
                "data/test",
                "tests/test",
                "tools/build_utils/test",
            ):
                path = repository / name
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("committed\n")
            subprocess.run([*git, "add", "."], check=True)
            subprocess.run(
                [
                    *git,
                    "-c",
                    "user.name=CI test",
                    "-c",
                    "user.email=ci@example.invalid",
                    "commit",
                    "--quiet",
                    "-m",
                    "test",
                ],
                check=True,
            )
            revision = subprocess.check_output(
                [*git, "rev-parse", "HEAD"], text=True
            ).strip()
            (repository / "CMakeLists.txt").write_text("uncommitted\n")
            with tempfile.TemporaryFile() as archive:
                runner.archive_source(repository, revision, archive)
                archive.seek(0)
                runner.extract_source(archive, source, revision)
            self.assertEqual((source / "CMakeLists.txt").read_text(), "committed\n")

    def test_archive_of_wrong_commit_rejected(self):
        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(ValueError):
            runner.extract_source(self.archive("b" * 40), Path(tmp), self.revision)

    def test_archive_cannot_escape_run_directory(self):
        for member, link in (("../escaped", None), ("escape", "/etc")):
            with self.subTest(member=member), tempfile.TemporaryDirectory() as tmp:
                with self.assertRaises(tarfile.FilterError):
                    runner.extract_source(
                        self.archive(self.revision, member, link),
                        Path(tmp),
                        self.revision,
                    )

    def test_incomplete_source_is_not_a_build(self):
        with tempfile.TemporaryDirectory() as tmp, self.assertRaises(ValueError):
            runner.extract_source(
                self.archive(self.revision, "README.md"), Path(tmp), self.revision
            )

    def test_missing_or_later_failed_report_cannot_pass(self):
        self.assertEqual(runner.final_status("Summary: compiled\n"), "UNKNOWN")
        self.assertEqual(runner.final_status("Status: OK\nStatus: FAILED\n"), "FAILED")
        self.assertEqual(runner.final_status("Status: OK\n"), "OK")

    def test_unaccepted_image_is_not_selected(self):
        with tempfile.TemporaryDirectory() as tmp:
            image_dir = Path(tmp) / "candidate/16V100"
            image_dir.mkdir(parents=True)
            (image_dir / "image.json").write_text(json.dumps(dict(verified=True)))
            with self.assertRaises(ValueError):
                runner.select_image(tmp)


if __name__ == "__main__":
    unittest.main()
