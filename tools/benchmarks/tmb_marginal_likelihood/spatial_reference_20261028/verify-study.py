"""Read-only identities for the declared spatial reference study."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def verify():
    frozen = read(HERE / "BASELINE.json")
    previous_gate = subprocess.run(["C:/Users/peterman.73/AppData/Local/Programs/Python/Python312/python.exe",
        str(HERE.parent / "functional_precision_20261028/verify-candidate.py")],
        capture_output=True, text=True)
    changed = [name for name, expected in frozen["source_sha256"].items() if digest(ROOT / name) != expected]
    producer_changed = [name for name, expected in frozen["producer_sha256"].items() if digest(HERE / name) != expected]
    installed = Path(frozen["installed_package"])
    runtime_changed = [name for name, expected in frozen["installed_sha256"].items()
                       if digest(installed / name) != expected]
    preceding = HERE.parent / "functional_precision_20261028"
    artifacts = read(preceding / "ARTIFACT_MANIFEST.json")
    previous_changed = [row["file"] for row in artifacts if digest(preceding / row["file"]) != row["sha256"]]
    branch = subprocess.check_output(["git", "-C", str(ROOT), "branch", "--show-current"], text=True).strip()
    head = subprocess.check_output(["git", "-C", str(ROOT), "rev-parse", "HEAD"], text=True).strip()
    index = subprocess.check_output(["git", "-C", str(ROOT), "diff", "--cached", "--name-only"], text=True).strip()
    output = {"source_changed": changed, "producer_changed": producer_changed,
              "installed_runtime_changed": runtime_changed, "preceding_evidence_changed": previous_changed,
              "branch": branch, "head": head, "index_empty": not index,
              "preceding_gate_and_additional_source_inventory_passed": previous_gate.returncode == 0}
    output["passed"] = not any((changed, producer_changed, runtime_changed, previous_changed, index)) and \
        branch == "codex/inference-dev" and head == frozen["head"] and previous_gate.returncode == 0
    return output


if __name__ == "__main__":
    try:
        result = verify()
        print(json.dumps(result, indent=2))
        sys.exit(0 if result["passed"] else 1)
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        print(json.dumps({"passed": False, "error": str(error)}))
        sys.exit(1)
