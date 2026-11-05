"""Check that delivery survives loss of access to the producer directory."""
from datetime import datetime, timedelta, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

BASE = Path(__file__).resolve().parents[1]


def module(name, relative):
    spec = importlib.util.spec_from_file_location(name, BASE / relative)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


delivery = module("delivery", "deployment/deliver-image.py")
server = module("delivery_server", "server.py")


class DeliveryTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        self.source = self.root / "producer"
        self.destination = self.root / "service"
        self.destination.mkdir()
        repository = self.source / "cache/repositories/cp2k"
        repository.mkdir(parents=True)
        self.git("init", "-q", str(repository))
        (repository / "tests").mkdir()
        (repository / "tests/do_regtest.py").write_text("# pinned harness\n")
        self.git("-C", str(repository), "add", ".")
        self.git("-C", str(repository), "-c", "user.name=Test", "-c", "user.email=test@example.org",
                 "commit", "-qm", "fixture")
        self.revision = self.git("-C", str(repository), "rev-parse", "HEAD").strip()
        self.image = self.source / "containers/software/cp2k/development/build/16V100/run.sif"
        self.image.parent.mkdir(parents=True)
        self.image.write_bytes(b"immutable image")
        self.record = dict(identity=dict(source_sha=self.revision), artifact=str(self.image),
                           sha256=self.checksum(self.image), build_verified=True, verified=True, published=True)
        self.image.with_suffix(".json").write_text(json.dumps(self.record))
        modules = dict(delivery_layout=types.SimpleNamespace(load_artifact=lambda *a, **k: self.record),
                       source_cache=types.SimpleNamespace(checksum=self.checksum))
        self.enterContext(patch.dict(sys.modules, modules))

    @staticmethod
    def git(*args):
        return subprocess.check_output(["git", *args], text=True, stderr=subprocess.PIPE)

    @staticmethod
    def checksum(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()

    def test_delivered_source_and_image_need_no_producer_access(self):
        result = delivery.deliver(self.source, self.destination)
        sidecar = Path(result["manifest"])
        record = json.loads(sidecar.read_text())
        self.assertEqual(record["producer"]["record"], self.record)
        self.assertEqual(record["producer"]["sha256"], self.checksum(self.image.with_suffix(".json")))
        self.assertEqual(delivery.deliver(self.source, self.destination), result)
        self.source.rename(self.root / "unavailable")
        self.assertEqual(self.checksum(record["artifact"]), self.record["sha256"])
        repository = self.destination / "cache/repositories/cp2k"
        self.assertFalse((repository / "objects/info/alternates").exists())
        self.assertEqual(self.git("--git-dir=" + str(repository), "show",
                                  self.revision + ":tests/do_regtest.py"), "# pinned harness\n")

    def test_missing_source_does_not_publish_sidecar(self):
        self.record["identity"]["source_sha"] = "0" * 40
        self.image.with_suffix(".json").write_text(json.dumps(self.record))
        with self.assertRaises(subprocess.CalledProcessError):
            delivery.deliver(self.source, self.destination)
        self.assertEqual(list(self.destination.glob("containers/**/*.json")), [])

    def test_failed_or_stale_delivery_cannot_reuse_old_success(self):
        for status, age in (("FAILED", 0), ("OK", 3 * 3600)):
            checked = datetime.now(timezone.utc) - timedelta(seconds=age)
            delivery.write(self.destination / "sync.json", dict(status=status, checked_at=checked.isoformat()))
            with patch.object(server, "PROJECT", self.destination), patch.object(server, "records") as records:
                with self.assertRaisesRegex(ValueError, "delivery failed or is stale"):
                    server.latest()
                records.assert_not_called()


if __name__ == "__main__":
    unittest.main()
