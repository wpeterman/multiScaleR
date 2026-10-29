"""Freeze guide validation with explicit unchanged-code coverage provenance."""
import datetime
import difflib
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
spec = importlib.util.spec_from_file_location("guide_audit", HERE / "verify-candidate.py")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)
result = audit.verify()
assert result["passed"], result
assert not (HERE / "final-source-evidence.json").exists()
admission = audit.read(HERE / "closeout-admission.json")
assert admission["R_process_count"] == 0 and admission["index_empty"]
age = datetime.datetime.now(datetime.timezone.utc) - datetime.datetime.fromisoformat(admission["time"])
assert 0 <= age.total_seconds() < 300
snapshot = audit.read(HERE / "snapshot01-input.json")
assert admission["branch"] == "codex/inference-dev" and admission["head"] == snapshot["head"]
patch = []
for name in snapshot["changed_from_baseline"] + snapshot["new_source"]:
    predecessor = HERE / "predecessor" / name
    if name in snapshot["changed_from_baseline"]:
        assert predecessor.is_file(), name
    before = predecessor.read_text(encoding="utf-8-sig").splitlines(keepends=True) if predecessor.is_file() else []
    after = (ROOT / name).read_text(encoding="utf-8-sig").splitlines(keepends=True)
    patch.extend(difflib.unified_diff(before, after, fromfile="predecessor/" + name, tofile=name))
(HERE / "review.patch").write_text("".join(patch), encoding="utf-8")
record = {"closed_at_UTC": datetime.datetime.now(datetime.timezone.utc).isoformat(),
          "audit": result, "source_sha256": snapshot["source_sha256"],
          "scope": "Frozen 318-file documentation candidate; new-archive examples/vignettes checked, complete numerical tests retained only for byte-identical code/test/data sources.",
          "owned_R_absent_confirmed_before_closeout": True,
          "new_numerical_fit_chain_dataset_or_predictive_replicate": False,
          "installed_tests_rerun": False, "visual_inspection_claimed": False,
          "original_render_harness_failure_retained": True}
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
