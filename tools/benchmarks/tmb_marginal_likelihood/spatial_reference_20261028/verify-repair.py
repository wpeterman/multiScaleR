"""Audit both the original frozen study and the pre-sampling repair."""
import json
from pathlib import Path
import subprocess
import sys
import hashlib

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def verify():
    baseline = read(HERE / "REPAIR_BASELINE.json")
    original = subprocess.run([sys.executable, str(HERE / "verify-study.py")], capture_output=True, text=True)
    changed = [name for name, expected in baseline["producer_sha256"].items() if digest(HERE / name) != expected]
    failure_changed = [name for name, expected in baseline["retained_failure_sha256"].items()
                       if digest(HERE / name) != expected]
    result = {"original_identity_gate_passed": original.returncode == 0,
              "repair_producer_changed": changed, "retained_failure_changed": failure_changed,
              "dry_repair_gate_passed": baseline["dry_repair_gate_passed"]}
    result["passed"] = original.returncode == 0 and not changed and not failure_changed and baseline["dry_repair_gate_passed"]
    return result


if __name__ == "__main__":
    try:
        output = verify()
        print(json.dumps(output, indent=2))
        sys.exit(0 if output["passed"] else 1)
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        print(json.dumps({"passed": False, "error": str(error)}))
        sys.exit(1)
