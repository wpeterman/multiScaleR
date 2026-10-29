"""Freeze/verify a supervisor continuation without rerunning completed work."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def original_gate():
    return subprocess.run([sys.executable, str(HERE / "verify-repair.py")], capture_output=True, text=True).returncode == 0


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--freeze", action="store_true")
args = parser.parse_args()
target = HERE / "CONTINUATION_BASELINE.json"
if args.freeze:
    assert not target.exists() and original_gate()
    completed = read(HERE / "null-a-v2.json")
    assert completed["completed"] and completed["sampling_completed"] and not completed["reference_qualified"]
    for name in ("null-b", "smooth-a", "smooth-b"):
        assert not list(HERE.glob(name + "-v2*")), "Pending phase already has outputs: " + name
    assert read(HERE / "unique-checkpoints-test01.json")["passed"]
    files = ["run-pending-v3.ps1", "CONTINUATION_PROTOCOL.md", "continuation.py",
        "null-a-v2.json", "null-a-v2-exit.json", "checkpoint-v2-running.json",
        "panel-v2-started.json", "REPAIR_BASELINE.json", "unique-checkpoints-test01.json"]
    record = {"sha256": {name: digest(HERE / name) for name in files},
              "remaining": ["null-b", "smooth-a", "smooth-b"],
              "prior_wall_seconds": read(HERE / "null-a-exit.json")["elapsed_seconds"] +
                                    read(HERE / "null-a-v2-exit.json")["elapsed_seconds"],
              "completed_first_fit_will_not_replay": True}
    assert all(record["sha256"].values())
    temporary = target.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2), encoding="utf-8")
    temporary.replace(target)
record = read(target)
changed = [name for name, expected in record["sha256"].items() if digest(HERE / name) != expected]
output = {"passed": not changed and original_gate(), "continuation_or_completed_evidence_changed": changed,
          "remaining_declared_slots": record["remaining"]}
print(json.dumps(output, indent=2))
sys.exit(0 if output["passed"] else 1)
