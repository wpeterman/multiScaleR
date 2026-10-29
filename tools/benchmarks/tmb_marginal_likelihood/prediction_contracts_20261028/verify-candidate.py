"""Read-only source/check/resource audit for this exact prediction-contract batch."""
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).rstrip("\r\n")


def digest(path, algorithm="sha256"):
    return hashlib.new(algorithm, path.read_bytes()).hexdigest() if path.is_file() else None


def verify():
    with (HERE / "check01/checked-source-md5.csv").open(encoding="utf-8-sig", newline="") as stream:
        checked = {row["file"].replace("\\", "/"): row["md5"] for row in csv.DictReader(stream)}
    changed = sorted(name for name, expected in checked.items()
                     if digest(ROOT / name, "md5") != expected)
    inventory = {name for name in git("ls-files").splitlines()
                 if not name.startswith("tools/") and (ROOT / name).is_file()}
    for folder in ("R", "man", "tests"):
        inventory.update(path.relative_to(ROOT).as_posix() for path in (ROOT / folder).rglob("*")
                         if path.is_file() and path.suffix.lower() not in (".o", ".dll", ".so", ".a"))
    inventory.update(name for name in git("ls-files", "--others", "--exclude-standard",
                                          "vignettes", "inst", "src").splitlines()
                     if (ROOT / name).is_file() and Path(name).suffix.lower() not in (".o", ".dll", ".so", ".a"))
    additional = sorted(inventory - set(checked))
    result = {"source_stable": not changed and not additional,
              "checked_files": len(checked), "changed": changed, "additional": additional,
              "branch": git("branch", "--show-current"), "head": git("rev-parse", "HEAD"),
              "index_empty": not git("diff", "--cached", "--name-only")}
    closed = ROOT / "tools/benchmarks/tmb_marginal_likelihood/unified_development_20261028/checkpoint.json"
    previous = json.loads(closed.read_text(encoding="utf-8"))
    evidence = {name: expected for name, expected in previous["pinned_sha256"].items()
                if name.startswith("tools/")}
    result["prior_evidence_changed"] = sorted(name for name, expected in evidence.items()
                                            if digest(ROOT / name) != expected)
    if not (HERE / "check04.json").is_file():
        result.update(state="running_or_attention_required", passed=False)
        return result
    completed = json.loads((HERE / "check04.json").read_text(encoding="utf-8-sig"))
    supplied = json.loads((HERE / "check04-input.json").read_text(encoding="utf-8-sig"))
    phases = [json.loads((HERE / (name + "-exit.json")).read_text(encoding="utf-8-sig"))
              for name in ("check04", "installed04", "supplement05")]
    rows = [json.loads(line) for name in ("check04", "installed04", "supplement05")
            for line in (HERE / (name + "-resources.jsonl")).read_text(encoding="utf-8-sig").splitlines()]
    logs = (HERE / "check04/R-CMD-check.log").read_text(encoding="utf-8-sig")
    tests = (HERE / "check04/installed-tests.log").read_text(encoding="utf-8-sig")
    counts = re.findall(r"\[ FAIL (\d+) \| WARN (\d+) \| SKIP (\d+) \| PASS (\d+) \]", tests)
    if not counts:
        raise ValueError("Installed test summary is missing")
    failed, warnings, skips, passed = map(int, counts[-1])
    supplement = json.loads((HERE / "supplement05.json").read_text(encoding="utf-8-sig"))
    cran_skips = re.findall(r"On CRAN \((\d+)\)", tests)
    recovered = bool(cran_skips) and int(cran_skips[-1]) == 7 and supplement["passed"] and supplement["archive_tests_unchanged"]
    additional_passes = supplement["additional_unique_assertions"] if recovered else 0
    resource_passed = bool(rows) and all(phase["exit_code"] == 0 and phase["operational_passed"]
        for phase in phases) and all(row["unknown"] == 0 and row["foreign"] == 0
        and row["actual_interpreters"] <= 4 and row["RAM_bytes"] > 12 * 2**30
        and row["disk_bytes"] > 5 * 2**30 for row in rows)
    archive_hash = digest(Path(supplied["archive"]))
    result.update(state="completed", check_passed=completed["passed"],
        check_status_OK=bool(re.search(r"(?m)^Status: OK\s*$", logs)),
        installed_tests={"passed": passed + additional_passes, "failures": failed + supplement["failures"],
          "warnings": warnings, "skips": skips - (7 if recovered else 0)},
        base_CRAN_mode_assertions=passed, supplementary_unique_assertions=additional_passes,
        supplemental_total_assertions=supplement["supplemental_assertions"], CRAN_skip_coverage_recovered=recovered,
        resource_passed=resource_passed, resource_samples=len(rows),
        maximum_actual_interpreters=max(row["actual_interpreters"] for row in rows),
        minimum_RAM_GiB=min(row["RAM_bytes"] for row in rows) / 2**30,
        minimum_disk_GiB=min(row["disk_bytes"] for row in rows) / 2**30,
        tests_run_separately=completed["tests_run_separately"],
        archive_identity=bool(archive_hash and archive_hash.upper() == supplied["SHA256"] == completed["archive_SHA256"]
            and digest(Path(supplied["check_root"]) / "multiScaleR_0.7.4.tar.gz").upper() == supplied["SHA256"]))
    result["passed"] = all((result["source_stable"], result["branch"] == "codex/inference-dev",
        result["head"] == previous["head"], result["index_empty"], not result["prior_evidence_changed"],
        result["check_passed"], result["check_status_OK"], resource_passed, failed == 0,
        result["archive_identity"], recovered))
    return result


if __name__ == "__main__":
    try:
        output = verify()
        print(json.dumps(output, indent=2))
        sys.exit(0 if output["passed"] else 1)
    except (OSError, ValueError, KeyError, subprocess.SubprocessError) as error:
        print(json.dumps({"passed": False, "error": str(error)}))
        sys.exit(1)
