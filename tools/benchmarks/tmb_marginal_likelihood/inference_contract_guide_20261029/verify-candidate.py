"""Audit documentation separately from retained identical-code test coverage."""
import hashlib
import html
import json
from pathlib import Path
import re
import sys
import tarfile

from snapshot import HERE, ROOT, digest, git, inventory


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def protected(name):
    return name in ("DESCRIPTION", "NAMESPACE") or name.startswith((
        "R/", "src/", "tests/", "data/", "man/")) or (name.startswith("inst/") and not name.startswith("inst/doc/"))


def archived_protected(path):
    result = {}
    with tarfile.open(path, "r:gz") as archive:
        for member in archive.getmembers():
            name = member.name.removeprefix("multiScaleR/")
            if not member.isfile() or not protected(name):
                continue
            value = archive.extractfile(member).read()
            if name == "DESCRIPTION":
                # R CMD build adds a creation timestamp; source DESCRIPTION is
                # checked byte-for-byte separately. Normalize only that field.
                value = b"\n".join(line for line in value.splitlines() if not line.startswith(b"Packaged:"))
            result[name] = hashlib.sha256(value).hexdigest()
    return result


def verify():
    snapshot = read(HERE / "snapshot01-input.json")
    baseline = read(HERE / "baseline.json")
    source = snapshot["source_sha256"]
    changed = [name for name, expected in source.items() if digest(ROOT / name) != expected]
    additional = sorted(set(inventory()) - set(source))
    protected_changed = [name for name, expected in baseline["source_sha256"].items()
                         if protected(name) and digest(ROOT / name) != expected]
    evidence_changed = []
    for folder, expected_manifest in baseline["closed_manifest_sha256"].items():
        previous = HERE.parent / folder
        if digest(previous / "ARTIFACT_MANIFEST.json") != expected_manifest:
            evidence_changed.append(folder + "/ARTIFACT_MANIFEST.json")
        for row in read(previous / "ARTIFACT_MANIFEST.json"):
            if digest(previous / row["file"]) != row["sha256"]:
                evidence_changed.append(folder + "/" + row["file"])
    unified = read(HERE.parent / "unified_development_20261028/checkpoint.json")
    historic_changed = [name for name, expected in unified["pinned_sha256"].items()
                       if name.startswith("tools/") and digest(ROOT / name) != expected]
    chunks_unchanged = {}
    for name in ("vignettes/quickstart.Rmd", "vignettes/bayesian-scales.Rmd"):
        before = (HERE / "predecessor" / name).read_text(encoding="utf-8-sig")
        after = (ROOT / name).read_text(encoding="utf-8-sig")
        pattern = r"(?ms)^```\{r[^\n]*\}\n.*?^```[ \t]*$"
        old_chunks, new_chunks = re.findall(pattern, before), re.findall(pattern, after)
        chunks_unchanged[name] = bool(old_chunks) and old_chunks == new_chunks
    render = read(HERE / "render02.json")
    rendered = Path(render["rendered"])
    text = " ".join(html.unescape(re.sub(r"<[^>]+>", " ", rendered.read_text(encoding="utf-8"))).split())
    content_valid = all(phrase in text for phrase in (
        "Inference Targets and Prediction Uncertainty", "not_computed", "expected_response",
        "Keep group policies and uncertainty explicit", "Check units, support and numerical precision"))
    output = {"source_stable": not changed and not additional, "source_files": len(source),
              "changed": changed, "additional": additional, "protected_source_changed": protected_changed,
              "closed_evidence_changed": evidence_changed, "historic_evidence_changed": historic_changed,
              "existing_vignette_code_unchanged": chunks_unchanged,
              "render_passed": render["passed"] and render["source_stable"],
              "rendered_content_valid": content_valid, "visual_inspection_claimed": False,
              "branch": git("branch", "--show-current"), "head": git("rev-parse", "HEAD"),
              "index_empty": not git("diff", "--cached", "--name-only")}
    if not (HERE / "check01.json").exists():
        output.update(passed=False, state="running_or_attention_required")
        return output
    complete = read(HERE / "check01.json")
    built = read(HERE / "build01.json")
    previous = HERE.parent / "spatial_source_contract_20261029"
    previous_build = read(previous / "build01.json")
    previous_complete = read(previous / "final-source-evidence.json")
    old_archive = Path(previous_build["archive"])
    new_archive = Path(built["archive"])
    archive_protected_same = archived_protected(old_archive) == archived_protected(new_archive)
    phases = [read(HERE / (name + "-exit.json")) for name in ("build01", "check01", "render02")]
    rows = [json.loads(line) for name in ("build01", "check01", "render01", "render02")
            for line in (HERE / (name + "-resources.jsonl")).read_text(encoding="utf-8-sig").splitlines()]
    resources_passed = bool(rows) and all(p["operational_passed"] and p["exit_code"] == 0 for p in phases) and all(
        r["unknown"] == 0 and r["foreign"] == 0 and r["actual_interpreters"] <= 4 and
        r["RAM_bytes"] > 12 * 2**30 and r["disk_bytes"] > 5 * 2**30 for r in rows)
    check_log = (HERE / "R-CMD-check.log").read_text(encoding="utf-8-sig")
    output.update(state="completed", check_status_OK=bool(re.search(r"(?m)^Status: OK\s*$", check_log)),
        check_passed=complete["passed"], resource_passed=resources_passed, resource_samples=len(rows),
        maximum_actual_interpreters=max(r["actual_interpreters"] for r in rows),
        minimum_RAM_GiB=min(r["RAM_bytes"] for r in rows) / 2**30,
        minimum_disk_GiB=min(r["disk_bytes"] for r in rows) / 2**30,
        new_archive_SHA256=digest(new_archive).upper(),
        archive_identity=digest(new_archive).upper() == built["SHA256"] == complete["archive_SHA256"],
        previous_archive_identity=digest(old_archive).upper() == baseline["previous_archive_SHA256"] == previous_build["SHA256"],
        archived_code_test_data_help_and_namespace_unchanged=archive_protected_same,
        retained_installed_test_evidence=previous_complete["audit"]["installed_tests"],
        retained_installed_tests_NOT_CRAN=previous_complete["audit"]["NOT_CRAN"],
        installed_tests_rerun=complete["installed_tests_rerun"],
        coverage_scope="New archive examples/vignettes checked; prior complete suite covers byte-identical executable/test/data sources, not a new-archive test rerun.")
    output["passed"] = all((output["source_stable"], not protected_changed, not evidence_changed,
        not historic_changed, all(chunks_unchanged.values()), output["render_passed"], content_valid,
        output["branch"] == "codex/inference-dev", output["head"] == baseline["head"], output["index_empty"],
        output["check_status_OK"], output["check_passed"], resources_passed, output["archive_identity"],
        output["previous_archive_identity"], archive_protected_same, previous_complete["audit"]["passed"],
        not complete["installed_tests_rerun"]))
    return output


if __name__ == "__main__":
    try:
        result = verify()
        print(json.dumps(result, indent=2))
        sys.exit(0 if result["passed"] else 1)
    except (OSError, ValueError, KeyError, tarfile.TarError) as error:
        print(json.dumps({"passed": False, "error": str(error)}))
        sys.exit(1)
