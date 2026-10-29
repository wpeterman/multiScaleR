"""Read-only identity assertions for the unified development checkpoint.

This verifies retained coverage; it never fits a model, refreshes a baseline,
changes the index, or writes into a closed study. Intentional development drift
requires a new, separately validated checkpoint rather than automatic acceptance.
"""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
COMPILED = {".o", ".dll", ".so", ".a", ".pyc"}


def git(*arguments):
    return subprocess.check_output(
        ["git", "-C", str(ROOT), *arguments], text=True, encoding="utf-8"
    ).rstrip("\r\n")


def digest(path, algorithm="sha256"):
    if not path.is_file():
        return None
    value = hashlib.new(algorithm)
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def differences(expected, actual):
    """Keep missing, changed and extra identities distinct for useful failures."""
    return {
        "missing": sorted(set(expected) - set(actual)),
        "changed": sorted(name for name in expected.keys() & actual.keys()
                          if expected[name] != actual[name]),
        "additional": sorted(set(actual) - set(expected)),
    }


def package_inventory():
    # Mirror the checked candidate's inventory. Build-excluded tools are outside
    # runtime coverage; pinned evidence/handoff files are audited separately.
    names = set(git("ls-files").splitlines())
    names = {name for name in names
             if not name.startswith("tools/") and (ROOT / name).is_file()}
    for folder in ("R", "man", "tests"):
        names.update(path.relative_to(ROOT).as_posix()
                     for path in (ROOT / folder).rglob("*")
                     if path.is_file() and path.suffix.lower() not in COMPILED)
    names.update(name for name in git(
        "ls-files", "--others", "--exclude-standard", "vignettes", "inst", "src"
    ).splitlines() if (ROOT / name).is_file()
                 and Path(name).suffix.lower() not in COMPILED)
    return names


def self_test():
    cases = [
        ({"a": "1"}, {"a": "1"}, {"missing": [], "changed": [], "additional": []}),
        ({"a": "1"}, {}, {"missing": ["a"], "changed": [], "additional": []}),
        ({"a": "1"}, {"a": "2"}, {"missing": [], "changed": ["a"], "additional": []}),
        ({"a": "1"}, {"a": "1", "b": "2"},
         {"missing": [], "changed": [], "additional": ["b"]}),
    ]
    for expected, actual, answer in cases:
        if differences(expected, actual) != answer:
            raise AssertionError("Identity-drift detection failed")
    return len(cases)


def verify():
    checkpoint = json.loads((HERE / "checkpoint.json").read_text(encoding="utf-8"))
    assertions = {}
    assertions["branch"] = git("branch", "--show-current") == checkpoint["branch"]
    assertions["head"] = git("rev-parse", "HEAD") == checkpoint["head"]
    assertions["index_empty"] = not git("diff", "--cached", "--name-only")
    version = next(line.split(":", 1)[1].strip()
                   for line in (ROOT / "DESCRIPTION").read_text().splitlines()
                   if line.startswith("Version:"))
    assertions["version"] = version == checkpoint["version"]
    refs = dict(line.split(" ", 1) for line in git(
        "for-each-ref", "--format=%(refname:short) %(objectname)", "refs/heads"
    ).splitlines())
    assertions["old_branch_refs_preserved"] = all(
        refs.get(name) == revision
        for name, revision in checkpoint["preserved_refs"].items())

    manifest = ROOT / checkpoint["checked_manifest"]
    with manifest.open(encoding="utf-8-sig", newline="") as stream:
        checked = {row["file"].replace("\\", "/"): row["md5"]
                   for row in csv.DictReader(stream)}
    current = {name: digest(ROOT / name, "md5") for name in package_inventory()}
    coverage = differences(checked, current)
    assertions["checked_source_count"] = len(checked) == checkpoint["checked_source_count"]
    assertions["checked_source_identity"] = not any(coverage.values())

    pinned = checkpoint["pinned_sha256"]
    pin_actual = {name: digest(ROOT / name) for name in pinned if (ROOT / name).is_file()}
    pin_changes = differences(pinned, pin_actual)
    assertions["pinned_source_and_evidence"] = not any(pin_changes.values())

    summary = (ROOT / checkpoint["check_summary"]).read_text(encoding="utf-8-sig")
    tests = (ROOT / checkpoint["installed_tests"]).read_text(encoding="utf-8-sig")
    assertions["retained_check_status"] = all(
        phrase in summary for phrase in ("0 errors", "0 warnings", "0 notes"))
    assertions["retained_test_counts"] = (
        "[ FAIL 0 | WARN 12 | SKIP 10 | PASS 4763 ]" in tests)
    assertions["drift_detection_self_tests"] = self_test() == 4

    # Existing dirty changes are preserved. New handoff paths are the only
    # intended addition. Untracked research directory contents are not inferred
    # to be safe to delete just because they are outside runtime coverage.
    before = set(checkpoint["baseline_git_status"])
    now = set(git("status", "--porcelain=v1", "--untracked-files=normal").splitlines())
    status_removed = sorted(before - now)
    status_added = sorted(now - before)
    allowed = "?? " + HERE.relative_to(ROOT).as_posix() + "/"
    assertions["prior_working_status_preserved"] = not status_removed and all(
        row == allowed for row in status_added)

    return {
        "passed": all(assertions.values()),
        "assertions": assertions,
        "checked_files": len(checked),
        "pinned_files": len(pinned),
        "source_drift": coverage,
        "pin_drift": pin_changes,
        "working_status_removed": status_removed,
        "working_status_added": status_added,
        "scope": "Read-only checkpoint identity; no numerical studies rerun",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--self-test", action="store_true",
                        help="Test missing/changed/additional detection without any file writes")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps({"passed": True, "drift_detection_cases": self_test()}))
        return 0
    result = verify()
    print(json.dumps(result, indent=2))
    return 0 if result["passed"] else 1


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, StopIteration, subprocess.SubprocessError) as error:
        print(json.dumps({"passed": False, "error": str(error)}))
        sys.exit(1)
