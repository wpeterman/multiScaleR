"""Read-only audit of the active precision batch; never reruns a producer."""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

from snapshot import HERE, ROOT, digest, git, inventory


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def verify():
    frozen = read(HERE / "snapshot01-input.json")
    source = frozen["source_sha256"]
    changed = sorted(name for name, expected in source.items() if digest(ROOT / name) != expected)
    additional = sorted(set(inventory()) - set(source))
    previous = HERE.parent / "prediction_contracts_20261028"
    closed = read(previous / "ARTIFACT_MANIFEST.json")
    evidence_changed = [row["file"] for row in closed if digest(previous / row["file"]) != row["sha256"]]
    unified = read(HERE.parent / "unified_development_20261028/checkpoint.json")
    historic_changed = [name for name, expected in unified["pinned_sha256"].items()
                       if name.startswith("tools/") and digest(ROOT / name) != expected]
    output = {"source_stable": not changed and not additional, "source_files": len(source),
              "changed": changed, "additional": additional, "closed_contract_evidence_changed": evidence_changed,
              "historic_evidence_changed": historic_changed, "branch": git("branch", "--show-current"),
              "head": git("rev-parse", "HEAD"), "index_empty": not git("diff", "--cached", "--name-only")}
    if not (HERE / "check01.json").is_file():
        output.update(passed=False, state="running_or_attention_required")
        return output
    completed = read(HERE / "check01.json")
    built = read(HERE / "build01.json")
    phases = [read(HERE / (name + "-exit.json")) for name in ("build01", "check01", "installed01")]
    rows = [json.loads(line) for name in ("build01", "check01", "installed01")
            for line in (HERE / (name + "-resources.jsonl")).read_text(encoding="utf-8-sig").splitlines()]
    resource_passed = bool(rows) and all(p["operational_passed"] and p["exit_code"] == 0 for p in phases) and all(
        r["unknown"] == 0 and r["foreign"] == 0 and r["actual_interpreters"] <= 4
        and r["RAM_bytes"] > 12 * 2**30 and r["disk_bytes"] > 5 * 2**30 for r in rows)
    tests = (HERE / "installed01.log").read_text(encoding="utf-8-sig")
    counts = re.findall(r"\[ FAIL (\d+) \| WARN (\d+) \| SKIP (\d+) \| PASS (\d+) \]", tests)
    if not counts:
        raise ValueError("Complete installed-test summary missing")
    failed, warnings, skips, passed = map(int, counts[-1])
    check_log = (HERE / "R-CMD-check.log").read_text(encoding="utf-8-sig")
    test_dir = Path(completed["installed_test_directory"])
    test_manifest = read(test_dir / "archive-test-manifest.json")
    tests_unchanged = all(digest(test_dir / row["file"]) == row["sha256"] for row in test_manifest)
    output.update(state="completed", check_passed=completed["passed"],
        check_status_OK=bool(re.search(r"(?m)^Status: OK\s*$", check_log)),
        installed_tests={"passed": passed, "failures": failed, "warnings": warnings, "skips": skips},
        NOT_CRAN=completed["NOT_CRAN"], archive_tests_unchanged=tests_unchanged,
        resource_passed=resource_passed, resource_samples=len(rows),
        maximum_actual_interpreters=max(r["actual_interpreters"] for r in rows),
        minimum_RAM_GiB=min(r["RAM_bytes"] for r in rows) / 2**30,
        minimum_disk_GiB=min(r["disk_bytes"] for r in rows) / 2**30,
        tests_run_separately=completed["tests_run_separately"],
        archive_identity=digest(Path(built["archive"])).upper() == built["SHA256"] == completed["archive_SHA256"])
    output["passed"] = all((output["source_stable"], output["branch"] == "codex/inference-dev",
        output["head"] == frozen["head"], output["index_empty"], not evidence_changed, not historic_changed,
        output["check_passed"], output["check_status_OK"], resource_passed, failed == 0,
        output["archive_identity"], tests_unchanged, completed["NOT_CRAN"]))
    return output


if __name__ == "__main__":
    try:
        result = verify()
        print(json.dumps(result, indent=2))
        sys.exit(0 if result["passed"] else 1)
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        print(json.dumps({"passed": False, "error": str(error)}))
        sys.exit(1)
