"""Exercise the SSH boundary, interrupted submissions, and failure reporting."""
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

BASE = Path(__file__).resolve().parents[1]


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, BASE / path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


entry = module("entry", "deployment/entry.py")
server = module("server", "server.py")
relay = module("relay", "relay.py")


class GatewayTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name).resolve()
        self.task = self.root / "state/runtime-tests/github-123-1"
        (self.task / "results").mkdir(parents=True)
        (self.root / "code").mkdir()
        (self.root / "code/version.json").write_text('{"version": 1}')
        self.request = dict(source_sha="a" * 40, artifact_sha256="b" * 64,
                            recipe_sha256="c" * 64, action_url="https://example.test/run/123",
                            server_version={"version": 1})
        server.write(self.task / "request.json", self.request)
        self.enterContext(patch.object(server, "ROOT", self.root))
        self.enterContext(patch.object(server, "STATE", self.root / "state"))

    def test_fixed_protocol(self):
        for command in ("inspect", "latest", "status github-123-1", "fetch github-123-1",
                        "start github-123-1 " + "b" * 64 + " " + "a" * 40 + " once"):
            with self.subTest(command=command):
                self.assertEqual(entry.parse_command(command), command.split())
        for command in ("", "bash", "inspect; id", "latest && id", "sftp", "scp -t /tmp/x",
                        "fetch ../../.ssh/id_ed25519", "status github-123-1 extra",
                        "start github-123-1 " + "b" * 64 + " " + "a" * 40 + " arbitrary"):
            with self.subTest(command=command), self.assertRaises(ValueError):
                entry.parse_command(command)

    def test_lost_submit_response_recovers_only_matching_held_job(self):
        (self.task / "submission.started").touch()
        replies = ["987|other-job|github-123-1\n456|cp2k-cusolvermp-github-123-1|github-123-1",
                   "PENDING|JobHeldUser"]
        with patch.object(server, "capture", side_effect=replies), patch.object(server.subprocess, "run") as run:
            result = server.describe(self.task)
        self.assertEqual(result["job"], "456")
        self.assertEqual(result["state"], "RELEASED")
        self.assertEqual((self.task / "job.id").read_text(), "456\n")
        self.assertEqual(run.call_args.args[0], ["scontrol", "release", "456"])

    def test_ambiguous_submission_does_not_resubmit(self):
        (self.task / "submission.started").touch()
        with patch.object(server, "capture", return_value=""), self.assertRaisesRegex(ValueError, "Unresolved"):
            server.describe(self.task)
        self.assertFalse((self.task / "job.id").exists())

    def test_later_manual_hold_is_not_released(self):
        (self.task / "job.id").write_text("456\n")
        (self.task / "job.released").write_text("456\n")
        with patch.object(server, "capture", return_value="PENDING|JobHeldUser"), patch.object(server, "release") as release:
            result = server.describe(self.task)
        self.assertFalse(result["terminal"])
        release.assert_not_called()

    def test_slurm_failure_cannot_publish_success(self):
        def report(task, request, error):
            (task / "results/report.txt").write_text("Status: " + ("FAILED" if error else "OK") + "\n")
        suite = types.SimpleNamespace(write_report=report)
        with patch.dict(sys.modules, cp2k_regtest_suite=suite):
            result = server.finalize(self.task, self.request, dict(job="456", state="FAILED",
                exit_code="1:0", user="stardust", elapsed="00:01:00", nodes="node"))
        self.assertEqual(result["status"], "FAILED")
        self.assertEqual(json.loads((self.task / "results/final.json").read_text()), result)
        self.assertEqual(json.loads((self.root / "state/completed" / ("b" * 64 + ".json")).read_text()), result)

    def test_same_image_skips_without_allocating(self):
        previous = dict(self.request, status="OK", terminal=True, run_id=self.task.name)
        server.write(self.task / "results/final.json", previous)
        (self.root / "state/completed").mkdir()
        server.write(self.root / "state/completed" / ("b" * 64 + ".json"), previous)
        modules = dict(cp2k_regtest=types.SimpleNamespace(), cp2k_regtest_suite=types.SimpleNamespace(),
                       delivery_layout=types.SimpleNamespace(load_artifact=None),
                       source_cache=types.SimpleNamespace(checksum=None))
        with patch.dict(sys.modules, modules), patch.object(server, "capture") as capture:
            result = server.start("github-124-1", "b" * 64, "a" * 40, "once")
        self.assertTrue(result["skipped"])
        capture.assert_not_called()

    def test_relay_replaces_stale_success_on_transport_error(self):
        previous_cwd = Path.cwd()
        self.addCleanup(os.chdir, previous_cwd)
        os.chdir(self.root)
        environment = dict(GITHUB_RUN_ID="124", GITHUB_RUN_ATTEMPT="1", GITHUB_SHA="a" * 40,
                           GITHUB_SERVER_URL="https://github.com", GITHUB_REPOSITORY="Stardust0831/cp2k")
        inspect = dict(uid=1478400690, container="ubuntu-submit.sif", home_ssh_visible=False, code_writable=False)
        latest = dict(source_sha="d" * 40, image_sha256="b" * 64)
        with patch.dict(os.environ, environment, clear=True), patch.object(sys, "argv", ["relay.py"]), \
                patch.object(relay, "remote", side_effect=[inspect, latest, RuntimeError("SSH unavailable")]), \
                patch("builtins.print"):
            self.assertEqual(relay.main(), 1)
        final = json.loads(Path("artifacts/latest.json").read_text())
        self.assertEqual(final["status"], "UNKNOWN")
        self.assertEqual(final["source_sha"], latest["source_sha"])
        self.assertIn("Status: UNKNOWN", Path("artifacts/report.txt").read_text())


if __name__ == "__main__":
    unittest.main()
