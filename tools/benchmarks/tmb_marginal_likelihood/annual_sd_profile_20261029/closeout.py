"""Freeze the annual-SD primitive after its fresh complete candidate gate."""
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
trace = audit.read(HERE / "unit01-trace.json")
assert not trace["global_fit_eligible"] and not trace["uncertainty"]["available"]
assert all(p["available"] and not p["global_fit_eligible"] and len(p["attempts"]) == 4 for p in trace["points"])
assert all(a["eligible_nuisance"] for p in trace["points"] for a in p["attempts"])
assert all(p["limits"]["completed"] and not p["limits"]["witnesses"] for p in trace["points"])
assert all(len(limit["attempts"]) == 3 for p in trace["points"] for limit in p["limits"]["references"].values())
assert all(a["sd_coordinate"] == "structurally_omitted_kernel_limit"
    for p in trace["points"] for limit in p["limits"]["references"].values() for a in limit["attempts"])
unit = audit.read(HERE / "unit01-test-factorized-scale-sd-profile.json")
assert unit["assertions"] == 235 and all(unit[k] == 0 for k in ("failures", "errors", "warnings", "skips"))
parity = audit.read(HERE / "parity01.json")
assert parity["passed"] and parity["RNG_unchanged"] and len(parity["fixed_vector_cases"]) == 6
assert all(parity["original_formals_unchanged"].values())
compatibility = audit.read(HERE / "compatibility01.json")
assert compatibility["passed"] and compatibility["source_stable"]
assert sum(row["assertions"] for row in compatibility["results"]) == 483
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
          "scope": "Frozen 322-file private annual log-scale SD trace candidate; any future source edits require separate validation.",
          "owned_R_absent_confirmed_before_closeout": True, "focused_assertions": 235,
          "compatibility_assertions": 483, "exact_response_fixed_vector_parity_cases": 6,
          "response_helper_formals_preserved": True, "global_fit_eligible": False, "interval_available": False,
          "all_declared_fixture_points_available": True, "all_twelve_nuisance_starts_eligible": True,
          "all_eighteen_limit_attempts_retained": True,
          "annual_SD_structurally_irrelevant_in_exact_kernel_limits": True,
          "public_fitter_prior_diagnostic_prediction_sources_unchanged": True,
          "new_sampler_or_calibration_panel": False,
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
