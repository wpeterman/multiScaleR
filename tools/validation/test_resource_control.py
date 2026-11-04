"""Exercise short/normal/unknown-inventory controller lifecycles without R."""
from pathlib import Path
import subprocess
import tempfile
import unittest


class ResourceLifecycleTests(unittest.TestCase):
    def test_recording_and_existing_unknown_inventory_guard(self):
        code = Path(__file__).resolve().parent
        for case in ("completed", "running", "unknown"):
            with self.subTest(case=case), tempfile.TemporaryDirectory(prefix="msr-resource-unit-") as output:
                result = subprocess.run(["pwsh", "-NoProfile", "-File", str(code / "test_resource_control.ps1"),
                    "-Engine", str(code / "resource_control.ps1"), "-Output", output, "-Case", case],
                    capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
