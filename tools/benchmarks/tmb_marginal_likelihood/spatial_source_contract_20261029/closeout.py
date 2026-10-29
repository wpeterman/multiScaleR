"""Freeze the spatial source preflight once, including its original failure."""
import datetime
import difflib
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
spec = importlib.util.spec_from_file_location("candidate_audit", HERE / "verify-candidate.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
result = audit.verify()
assert result["passed"], result
assert not (HERE / "final-source-evidence.json").exists()
snapshot = audit.read(HERE / "snapshot01-input.json")
patch = []
for name in snapshot["changed_from_baseline"] + snapshot["new_source"]:
    predecessor = HERE / "predecessor" / name
    if name in snapshot["changed_from_baseline"]:
        assert predecessor.is_file(), name
    before = predecessor.read_text(encoding="utf-8-sig").splitlines(keepends=True) if predecessor.is_file() else []
    after = (ROOT / name).read_text(encoding="utf-8-sig").splitlines(keepends=True)
    patch.extend(difflib.unified_diff(before, after, fromfile="predecessor/" + name, tofile=name))
(HERE / "review.patch").write_text("".join(patch), encoding="utf-8")
unit = audit.read(HERE / "unit02-test-spatial-source-report.json")
assert unit["assertions"] == 64 and all(unit[k] == 0 for k in ("failures", "errors", "warnings", "skips"))
compatibility = audit.read(HERE / "compatibility01.json")
assert compatibility["passed"] and compatibility["source_stable"]
assert sum(row["assertions"] for row in compatibility["results"]) == 263
assert audit.read(HERE / "documentation01.json")["namespace_unchanged"]
admission = audit.read(HERE / "closeout-admission.json")
assert admission["R_process_count"] == 0 and admission["index_empty"]
assert admission["branch"] == "codex/inference-dev" and admission["head"] == snapshot["head"]
age = datetime.datetime.now(datetime.timezone.utc) - datetime.datetime.fromisoformat(admission["time"])
assert 0 <= age.total_seconds() < 300
resources = [json.loads(line) for path in HERE.glob("*-resources.jsonl")
             for line in path.read_text(encoding="utf-8-sig").splitlines()]
assert resources and all(row["unknown"] == 0 and row["foreign"] == 0 and row["actual_interpreters"] <= 4
    and row["RAM_bytes"] > 12 * 2**30 and row["disk_bytes"] > 5 * 2**30 for row in resources)
record = {"closed_at_UTC": datetime.datetime.now(datetime.timezone.utc).isoformat(),
          "audit": result, "source_sha256": snapshot["source_sha256"],
          "owned_R_absent_confirmed_before_closeout": True,
          "scope": "Frozen 317-file source candidate; later source edits require separate validation.",
          "focused_assertions": unit["assertions"], "compatibility_assertions": 263,
          "fit_eligibility": False, "historical_raster_lineage": "unverified_from_compact_inputs",
          "historical_distance_unit_agreement": "unverified_from_compact_inputs",
          "actual_valid_cell_and_kernel_coverage": "not_evaluated",
          "new_sampler_fit_or_calibration_panel": False,
          "existing_fitter_prior_gate_public_interface_and_prediction_sources_unchanged": True,
          "all_phase_resource_samples": len(resources),
          "all_phase_maximum_actual_interpreters": max(row["actual_interpreters"] for row in resources),
          "all_phase_minimum_RAM_GiB": min(row["RAM_bytes"] for row in resources) / 2**30,
          "all_phase_minimum_disk_GiB": min(row["disk_bytes"] for row in resources) / 2**30}
target = HERE / "final-source-evidence.json"
temporary = target.with_suffix(".json.tmp")
temporary.write_text(json.dumps(record, indent=2), encoding="utf-8")
temporary.replace(target)
manifest = [{"file": path.relative_to(HERE).as_posix(), "bytes": path.stat().st_size,
             "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for path in sorted(HERE.rglob("*")) if path.is_file()
            and "__pycache__" not in path.parts and not path.name.endswith(".tmp")
            and path.name != "ARTIFACT_MANIFEST.json"]
(HERE / "ARTIFACT_MANIFEST.json.tmp").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
(HERE / "ARTIFACT_MANIFEST.json.tmp").replace(HERE / "ARTIFACT_MANIFEST.json")
print(json.dumps({"passed": True, "artifacts": len(manifest),
                  "artifact_MiB": sum(row["bytes"] for row in manifest) / 2**20, "audit": result}))
