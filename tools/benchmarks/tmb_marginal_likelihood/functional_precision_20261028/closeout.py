"""Write one compact, final, immutable source/evidence checkpoint."""
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
    before = predecessor.read_text(encoding="utf-8-sig").splitlines(keepends=True) if predecessor.is_file() else []
    after = (ROOT / name).read_text(encoding="utf-8-sig").splitlines(keepends=True)
    patch.extend(difflib.unified_diff(before, after, fromfile="predecessor/" + name, tofile=name))
(HERE / "review.patch").write_text("".join(patch), encoding="utf-8")
record = {"closed_at_UTC": datetime.datetime.now(datetime.timezone.utc).isoformat(),
          "audit": result, "source_sha256": snapshot["source_sha256"],
          "owned_R_absent_confirmed_before_closeout": True,
          "scope": "Frozen 311-file source only; later deliberate edits require new validation."}
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
