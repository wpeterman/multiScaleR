"""Independently accept frozen source/archive/test identities and reviewed ledger changes."""
import argparse
from collections import Counter
import datetime as dt
import json
from pathlib import Path
import re
import subprocess

from make_plan import sha
from package_batch import read, write
from phase_runner import verify_inputs


def skipped(text):
    """Compare skip identities and reasons, excluding reporter decorations."""
    if "══ Skipped" not in text:
        return ""
    value = text.split("══ Skipped", 1)[1].splitlines()[1:]
    records = []
    for line in value:
        if line.startswith("══") or "[ FAIL " in line:
            break
        if line.strip():
            records.append(line.strip())
    return "\n".join(records)


def prior_test_log(parent):
    """Read named historical phases or the reusable batch's installed suite."""
    parent = Path(parent)
    previous = read(parent / "complete01.json")
    phase = previous.get("test_phase", "installed01")
    return (parent / (phase + ".log")).read_text(encoding="utf-8-sig")


def run(path):
    config = read(path)
    repo, output = Path(config["repo"]), Path(config["output"])
    parent = Path(config["parent_output"])
    baseline, complete, archive, normal = [read(output / name) for name in
        ("baseline01.json", "complete01.json", "archive01.json", "normal-check-complete01.json")]
    assert complete["passed"] and complete["tests_complete"] and complete["check_accepted"]
    assert normal["normal_check_completed"] and normal["normal_vignette_enabled"]
    assert subprocess.check_output(["git", "branch", "--show-current"], cwd=repo).decode().strip() == "inference-dev"
    assert subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo).decode().strip() == baseline["head"]
    assert not subprocess.check_output(["git", "diff", "--cached", "--name-only"], cwd=repo).strip()
    assert all(sha(repo / name) == value for name, value in baseline["root_source_sha256"].items())
    candidate = Path(baseline["candidate_package"])
    assert all(sha(candidate / name) == value for name, value in baseline["source_sha256"].items())
    assert archive["source_identity_passed"] and archive["normal_vignette_index_present"]
    assert sha(archive["archive"]) == archive["archive_sha256"] == normal["archive_sha256"] == complete["archive_sha256"]
    for phase in ("build01", "check01", "installed01"):
        plan = read(output / (phase + "-plan01.json"))
        verify_inputs(plan)
        assert read(output / (phase + "-verified.json"))["inputs_verified"]
        assert read(output / (phase + "-owned-seal01.json"))["verified_phase_owned_processes_remaining"] == 0
        resources = [json.loads(line) for line in (output / (phase + "-resources.jsonl")).read_text(encoding="utf-8-sig").splitlines()]
        assert resources
        for item in resources:
            assert item["unknown"] == 0 and item["owned_interpreters"] <= 3 and item["threads"] == 1
            assert item["disk_bytes"] > 5 * 2**30 and item["RAM_bytes"] > 8 * 2**30
            assert item["owned_priority"] == "Normal"
            assert item["actual_interpreters"] == item["owned_interpreters"] + item["foreign"]
    assert all(sha(Path(normal["installed_library"]) / name) == value for name, value in normal["installed_code_sha256"].items())
    assert all(sha(Path(normal["exact_archive_tests"]) / name) == value for name, value in normal["test_sha256"].items())
    log = (output / "installed01.log").read_text(encoding="utf-8-sig")
    values = re.findall(r"\[ FAIL (\d+) \| WARN (\d+) \| SKIP (\d+) \| PASS (\d+) \]", log)[-1]
    counts = dict(zip(("failed", "warnings", "skips", "passed"), map(int, values)))
    assert counts["failed"] == counts["warnings"] == 0
    previous_log = prior_test_log(parent)
    assert skipped(log) == skipped(previous_log)
    prior = read(config["parent_ledger"])
    current = read(output / "test-counts01.json")
    assert set(current) == set(prior) | set(config["new_test_files"])
    for name, before in prior.items():
        expected = dict(before)
        for key, delta in config["reviewed_ledger_deltas"].get(name, {}).items():
            expected[key] += delta
        assert current[name] == expected, (name, current[name], expected)
    focused = read(config["focused_result"])
    assert focused["counts"]["failed"] == focused["counts"]["errors"] == focused["counts"]["warnings"] == 0
    for name in config["new_test_files"]:
        assert current[name] == focused["by_file"][name]
    total = Counter()
    for value in current.values():
        total.update(value)
    for name in ("passed", "failed", "warnings", "skips"):
        assert total[name] == counts[name]
    assert total["errors"] == 0
    for marker in config["scientific_markers"]:
        assert log.count(marker) == 1, marker
    assert len(config["scientific_markers"]) == 14
    assert sum(log.count(value) for value in ("UNMARKED_NATIVE_PIPE_GUARD retained", "UNMARKED_NATIVE_PIPE_ACCEPTED native adapter checked")) == 1
    result = dict(accepted=True, verified_at=dt.datetime.now().astimezone().isoformat(),
        verifier_sha256=sha(Path(__file__)), configuration_sha256=sha(path),
        branch="inference-dev", head=baseline["head"], check_status=complete["check_status"],
        normal_vignette_enabled=True, full_validation=counts, archive_sha256=archive["archive_sha256"],
        source_files=len(baseline["source_sha256"]), source_archive_installed_tests_dependency_verified=True,
        all14scientific_markers_preserved=True, skip_identities_preserved=True,
        reviewed_ledger_deltas=config["reviewed_ledger_deltas"], verified_owned_R_remaining=0,
        private_admission_guards_preserved=True, timing_comparison_eligible=False,
        as_cran_incoming_status="unresolved_external_metadata_timeouts_from_full67")
    write(config["result"], result)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    run(parser.parse_args().input)
