"""Retain the failed pre-sampling execution, then freeze its operational repair."""
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path, algorithm="sha256"):
    return hashlib.new(algorithm, path.read_bytes()).hexdigest()


target = HERE / "REPAIR_BASELINE.json"
assert not target.exists(), "Retain existing repair baseline"
subprocess.run([sys.executable, str(HERE / "verify-study.py")], check=True)
dry = read(HERE / "repair-unit01.json")
assert dry["passed"] and dry["all_four_numeric_starts_identical"] and dry["tmbstan_initial_mapping_valid"]
assert digest(HERE / "run-fit-v2.R", "md5") == dry["repaired_driver_md5"]
assert read(HERE / "atomic-writer-validation02.json")["replacement_passed"]
failed = read(HERE / "null-a.json")
assert not failed["completed"] and not failed.get("sampling_completed", False)
assert failed["source_unchanged"] and failed["producer_unchanged"]
producers = ["OPS_REPAIR_PROTOCOL.md", "run-fit-v2.R", "resource-controller-v2.ps1",
             "run-panel-v2.ps1", "probe-repair.R", "verify-repair.py", "freeze-repair.py"]
retained = ["BASELINE.json", "null-a.json", "null-a-started.json", "null-a-exit.json",
            "null-a-process.json", "null-a.log", "null-a.err", "null-a-resources.jsonl",
            "panel-started.json", "checkpoint-running.json", "checkpoint-running.json.tmp",
            "unit01.json", "unit01.log", "unit02.json", "unit02.log",
            "atomic-writer-probe01-failure.json", "atomic-writer-validation.json",
            "atomic-writer-validation.json.tmp", "pre-repair-probe/resource-controller-v2.ps1",
            "repair-unit01.json", "atomic-writer-validation02.json"]
record = {"frozen_at_UTC": datetime.datetime.now(datetime.timezone.utc).isoformat(),
          "dry_repair_gate_passed": True, "prior_sampling_calls_completed": 0,
          "producer_sha256": {name: digest(HERE / name) for name in producers},
          "retained_failure_sha256": {name: digest(HERE / name) for name in retained},
          "prior_failed_wall_seconds": read(HERE / "null-a-exit.json")["elapsed_seconds"],
          "planned_scientific_slots": ["null-a", "null-b", "smooth-a", "smooth-b"],
          "scope": "Pre-sampling operational repair; target, priors, numeric starts, seeds, controls and statistical gates unchanged."}
temporary = target.with_suffix(".json.tmp")
temporary.write_text(json.dumps(record, indent=2), encoding="utf-8")
temporary.replace(target)
print(json.dumps({"frozen": True, "repair_producers": len(producers), "retained_failure_files": len(retained)}))
