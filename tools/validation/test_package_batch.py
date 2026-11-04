"""Small non-R archive and validation-policy regression tests."""
import io
import json
from pathlib import Path
import tarfile
import tempfile
import unittest

from package_batch import archive_tests, sha
from verify_package_batch import prior_test_log


class ArchiveTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="msr-archive-unit-")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.package = self.root / "candidate"
        (self.package / "R").mkdir(parents=True)
        (self.package / "R/code.R").write_bytes(b"identity(1)\n")
        self.sources = {"R/code.R": sha(self.package / "R/code.R")}
        self.members = {"multiScaleR/R/code.R": b"identity(1)\n",
            "multiScaleR/tests/testthat.R": b"stop('unit fixture never executes R')\n",
            "multiScaleR/build/vignette.rds": b"index fixture",
            "multiScaleR/inst/doc/feature-compatibility.html": b"html fixture"}

    def archive(self, members=None, link=None):
        archive = self.root / "package.tar.gz"
        with tarfile.open(archive, "w:gz") as stream:
            for name, data in (members or self.members).items():
                item = tarfile.TarInfo(name); item.size = len(data)
                stream.addfile(item, io.BytesIO(data))
            if link:
                item = tarfile.TarInfo(link); item.type = tarfile.SYMTYPE; item.linkname = "outside"
                stream.addfile(item)
        return archive

    def test_exact_source_and_tests_only(self):
        result = archive_tests(self.archive(), self.package, self.sources, self.root / "tests")
        self.assertTrue(result["source_identity_passed"])
        self.assertEqual(list((self.root / "tests").iterdir()), [self.root / "tests/testthat.R"])

    def test_traversal_and_links_are_rejected_before_extraction(self):
        for name in ("multiScaleR/tests/../../escape", "other/tests/test.R", "multiScaleR/tests/C:/x", "multiScaleR/tests/a\\b"):
            with self.subTest(name=name):
                with self.assertRaises(ValueError):
                    archive_tests(self.archive({name: b"bad"}), self.package, {}, self.root / "tests")
                self.assertFalse((self.root / "tests").exists())
        with self.assertRaises(ValueError):
            archive_tests(self.archive(link="multiScaleR/tests/link"), self.package, self.sources, self.root / "tests")

    def test_source_drift_or_missing_vignette_blocks_extraction(self):
        for change in ("source", "vignette"):
            members = dict(self.members)
            if change == "source": members["multiScaleR/R/code.R"] = b"identity(2)\n"
            else: del members["multiScaleR/build/vignette.rds"]
            with self.assertRaises(ValueError):
                archive_tests(self.archive(members), self.package, self.sources, self.root / "tests")
            self.assertFalse((self.root / "tests").exists())

    def test_full_profile_preserves_all_tests_and_gates(self):
        profile = json.loads(Path(__file__).with_name("profiles.json").read_text())
        self.assertEqual(profile["full"]["exact_archive_installed_tests"], "all")
        self.assertTrue(profile["full"]["normal_vignette_build"])
        self.assertEqual(profile["full"]["scientific_markers_required"], 14)
        self.assertEqual(profile["full"]["maximum_owned_socket_workers"], 2)

    def test_parent_log_supports_named_and_reusable_batches(self):
        for metadata, phase in (({"test_phase": "tests04"}, "tests04"),
                                ({"tests_complete": True}, "installed01")):
            with self.subTest(metadata=metadata):
                (self.root / "complete01.json").write_text(json.dumps(metadata))
                log = self.root / (phase + ".log")
                log.write_text("retained test ledger", encoding="utf-8-sig")
                self.assertEqual(prior_test_log(self.root), "retained test ledger")


if __name__ == "__main__":
    unittest.main()
