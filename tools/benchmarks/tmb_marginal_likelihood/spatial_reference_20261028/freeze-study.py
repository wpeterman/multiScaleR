"""Freeze identities before admitting the four predeclared fits."""
import datetime
import hashlib
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path, algorithm="sha256"):
    return hashlib.new(algorithm, path.read_bytes()).hexdigest()


target = HERE / "BASELINE.json"
assert not target.exists(), "Retain existing frozen study"
preceding = HERE.parent / "functional_precision_20261028"
subprocess.run(["C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe",
                str(preceding / "verify-candidate.py")], check=True)
unit = read(HERE / "unit04.json")
assert unit["passed"] and unit["failures"] == unit["errors"] == unit["warnings"] == unit["skips"] == 0
for name, expected in unit["source"].items():
    assert digest(ROOT / name, "md5") == expected, name
plan = read(HERE / "RUNLIST.json")
assert digest(HERE / "PROTOCOL.md") == plan["protocol_sha256"].lower()
qualified = read(preceding / "check01.json")
snapshot = read(preceding / "snapshot01-input.json")
installed = Path(qualified["workspace"]) / "multiScaleR.Rcheck/multiScaleR"
producers = ["adapter.R", "run-fit.R", "PROTOCOL.md", "RUNLIST.json", "test-adapter.R",
             "validate-adapter.R", "resource-controller.ps1", "run-panel.ps1", "verify-study.py", "freeze-study.py"]
record = {"frozen_at_UTC": datetime.datetime.now(datetime.timezone.utc).isoformat(),
          "head": snapshot["head"], "source_sha256": snapshot["source_sha256"],
          "producer_sha256": {name: digest(HERE / name) for name in producers},
          "installed_package": str(installed), "archive_SHA256": qualified["archive_SHA256"],
          "installed_sha256": {p.relative_to(installed).as_posix(): digest(p) for p in sorted(installed.rglob("*"))
                               if p.is_file() and "doc" not in p.relative_to(installed).parts},
          "adapter_assertions": unit["assertions"], "planned_runs": [row["id"] for row in plan["runs"]]}
temporary = target.with_suffix(".json.tmp")
temporary.write_text(json.dumps(record, indent=2), encoding="utf-8")
temporary.replace(target)
print(json.dumps({"frozen": True, "source_files": len(record["source_sha256"]),
                  "installed_files": len(record["installed_sha256"]), "adapter_assertions": unit["assertions"]}))
