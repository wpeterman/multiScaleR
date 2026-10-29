"""Freeze the completed four-slot study; never rerun a numerical producer."""
import datetime
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


assert not (HERE / "final-source-evidence.json").exists(), "Retain closed study"
subprocess.run([sys.executable, str(HERE / "continuation.py")], check=True)
analysis = read(HERE / "paired-analysis.json")
assert analysis["completed_runs"] == 4 and analysis["pending_runs"] == analysis["failed_runs"] == 0
assert analysis["identity_gate_passed"]
unit = read(HERE / "analysis-unit02.json")
assert unit["passed"] and digest(HERE / "analyze.py").upper() == unit["analysis_sha256"]
assert digest(HERE / "test-analysis.py").upper() == unit["test_sha256"]
rows = []
for name in ("null-a", "null-a-v2", "null-b-v2", "smooth-a-v2", "smooth-b-v2"):
    rows.extend(json.loads(line) for line in (HERE / (name + "-resources.jsonl")).read_text(encoding="utf-8-sig").splitlines())
assert rows and all(row["unknown"] == row["foreign"] == 0 and row["actual_interpreters"] <= 4
                    and row["RAM_bytes"] > 12 * 2**30 and row["disk_bytes"] > 5 * 2**30 for row in rows)
frozen = read(HERE / "BASELINE.json")
record = {"closed_at_UTC": datetime.datetime.now(datetime.timezone.utc).isoformat(),
    "source_sha256": frozen["source_sha256"], "installed_sha256": frozen["installed_sha256"],
    "archive_SHA256": frozen["archive_SHA256"], "all_four_completed": True,
    "original_reference_flags_retained": True, "owned_R_absence_confirmed_before_closeout": True,
    "resource_samples": len(rows), "maximum_actual_interpreters": max(r["actual_interpreters"] for r in rows),
    "minimum_RAM_GiB": min(r["RAM_bytes"] for r in rows) / 2**30,
    "minimum_disk_GiB": min(r["disk_bytes"] for r in rows) / 2**30,
    "wall_seconds": read(HERE / "panel-v3-complete.json")["cumulative_seconds"],
    "result_sha256": {name: digest(HERE / name) for name in
        ["null-a.json", "null-a-v2.json", "null-b-v2.json", "smooth-a-v2.json", "smooth-b-v2.json", "paired-analysis.json"]},
    "scope": "Frozen original study only; later source/report correction requires separate provenance and never overwrites these outcomes."}
temporary = HERE / "final-source-evidence.json.tmp"
temporary.write_text(json.dumps(record, indent=2), encoding="utf-8")
temporary.replace(HERE / "final-source-evidence.json")
manifest = [{"file": path.relative_to(HERE).as_posix(), "bytes": path.stat().st_size, "sha256": digest(path)}
    for path in sorted(HERE.rglob("*")) if path.is_file() and "__pycache__" not in path.parts
    and path.name != "ARTIFACT_MANIFEST.json"]
temporary = HERE / "ARTIFACT_MANIFEST.json.tmp"
temporary.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
temporary.replace(HERE / "ARTIFACT_MANIFEST.json")
print(json.dumps({"closed": True, "artifacts": len(manifest),
    "artifact_MiB": sum(row["bytes"] for row in manifest) / 2**20,
    "resource_samples": len(rows), "minimum_RAM_GiB": record["minimum_RAM_GiB"],
    "minimum_disk_GiB": record["minimum_disk_GiB"], "wall_seconds": record["wall_seconds"]}))
