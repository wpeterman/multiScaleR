"""Read-only final source/evidence audit. Never starts R or a numerical producer."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = HERE.parent
CLOSED = ["prediction_contracts_20261028", "functional_precision_20261028",
          "spatial_reference_20261028", "probability_precision_20261029",
          "response_sd_profile_20261029", "spatial_source_contract_20261029",
          "inference_contract_guide_20261029", "annual_sd_profile_20261029"]


def read(path):
    return json.loads(path.read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).strip()


def verify():
    latest = read(BASE / "annual_sd_profile_20261029/final-source-evidence.json")
    source = latest["source_sha256"]
    changed = [name for name, expected in source.items() if digest(ROOT / name) != expected]
    files = set(git("ls-files").splitlines())
    for folder in ("R", "man", "tests"):
        files.update(path.relative_to(ROOT).as_posix() for path in (ROOT / folder).rglob("*") if path.is_file())
    files.update(git("ls-files", "--others", "--exclude-standard", "vignettes", "inst", "src").splitlines())
    files = {name for name in files if (ROOT / name).is_file() and not name.startswith(("tools/", ".git/"))
             and Path(name).suffix.lower() not in (".o", ".dll", ".so", ".a")}
    additional = sorted(files - set(source))
    evidence_changed = []
    manifests = {}; evidence_files = evidence_bytes = 0
    for folder in CLOSED:
        directory = BASE / folder
        manifest = directory / "ARTIFACT_MANIFEST.json"
        manifests[folder] = digest(manifest)
        for row in read(manifest):
            evidence_files += 1; evidence_bytes += row["bytes"]
            if digest(directory / row["file"]) != row["sha256"]:
                evidence_changed.append(folder + "/" + row["file"])
    unified = read(BASE / "unified_development_20261028/checkpoint.json")
    historic_changed = [name for name, expected in unified["pinned_sha256"].items()
                       if name.startswith("tools/") and digest(ROOT / name) != expected]
    probability = read(BASE / "probability_precision_20261029/final-source-evidence.json")
    archive = read(BASE / "annual_sd_profile_20261029/build01.json")
    description = (ROOT / "DESCRIPTION").read_text()
    checkpoint_changed = []
    final_manifest = HERE / "FINAL_MANIFEST.json"
    if final_manifest.exists():
        for row in read(final_manifest):
            if digest(HERE / row["file"]) != row["sha256"]:
                checkpoint_changed.append(row["file"])
    static_changed = []; manifest_changed = []
    if (HERE / "final-state.json").exists():
        state = read(HERE / "final-state.json")
        static_changed = [name for name, expected in state["static_documents_sha256"].items()
                          if digest(ROOT / name) != expected]
        manifest_changed = [name for name, expected in state["checkpoint_closed_manifests_sha256"].items()
                            if manifests.get(name) != expected]
    output = {"source_files": len(source), "source_changed": changed, "additional_source": additional,
              "closed_studies": len(CLOSED), "closed_evidence_changed": evidence_changed,
              "historic_evidence_changed": historic_changed, "closed_manifest_sha256": manifests,
              "evidence_files_excluding_manifests": evidence_files, "evidence_MiB": evidence_bytes / 2**20,
              "branch": git("branch", "--show-current"), "head": git("rev-parse", "HEAD"),
              "index_empty": not git("diff", "--cached", "--name-only"),
              "version_0_7_4": "Version: 0.7.4\n" in description,
              "latest_complete_suite": latest["audit"]["installed_tests"],
              "latest_check_status_OK": latest["audit"]["check_status_OK"],
              "archive_SHA256": archive["SHA256"],
              "latest_archive_identity": digest(Path(archive["archive"])).upper() == archive["SHA256"],
              "smooth_reference_qualified_only_under_declared_numerical_gates": probability["selected_smooth_reference_passed_declared_gates"],
              "divergent_null_reference_remains_unqualified": probability["divergent_null_reference_remains_unqualified"],
              "SD_traces_global_fit_eligible": latest["global_fit_eligible"],
              "SD_intervals_available": latest["interval_available"],
              "interval_design_sha256": digest(BASE / "likelihood_uncertainty_design_20261029/DESIGN.md"),
              "final_checkpoint_files_changed": checkpoint_changed,
              "final_static_documents_changed": static_changed,
              "closed_manifests_changed_since_final_checkpoint": manifest_changed,
              "scope": "Read-only source/evidence integrity; no tests, fitting, sampling or calibration rerun."}
    output["passed"] = all((not changed, not additional, not evidence_changed, not historic_changed,
        output["branch"] == "codex/inference-dev", output["head"] == "87bf8779fe571921f2fd5d3e45c0f9e478530022",
        output["index_empty"], output["version_0_7_4"], output["latest_check_status_OK"],
        output["latest_complete_suite"]["failures"] == 0, output["latest_archive_identity"],
        output["divergent_null_reference_remains_unqualified"], not output["SD_traces_global_fit_eligible"],
        not output["SD_intervals_available"], not checkpoint_changed, not static_changed, not manifest_changed))
    return output


if __name__ == "__main__":
    result = verify()
    print(json.dumps(result, indent=2))
    sys.exit(0 if result["passed"] else 1)
