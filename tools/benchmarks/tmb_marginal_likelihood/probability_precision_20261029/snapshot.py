"""Freeze a dependency-complete small source snapshot without Git staging."""
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]


def git(*args):
    return subprocess.check_output(["git", "-C", str(ROOT), *args], text=True).rstrip("\r\n")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.is_file() else None


def inventory():
    files = set(git("ls-files").splitlines())
    for folder in ("R", "man", "tests"):
        files.update(p.relative_to(ROOT).as_posix() for p in (ROOT / folder).rglob("*") if p.is_file())
    files.update(git("ls-files", "--others", "--exclude-standard", "vignettes", "inst", "src").splitlines())
    return sorted(name for name in files if (ROOT / name).is_file()
                  and not name.startswith(("tools/", ".git/"))
                  and Path(name).suffix.lower() not in (".o", ".dll", ".so", ".a"))


if __name__ == "__main__":
    destination = HERE / "snapshot01-input.json"
    if destination.exists():
        raise RuntimeError("Retain existing snapshot attempt")
    if git("branch", "--show-current") != "codex/inference-dev" or git("diff", "--cached", "--name-only"):
        raise RuntimeError("Branch/index changed")
    baseline = json.loads((HERE / "baseline.json").read_text(encoding="utf-8-sig"))["source_sha256"]
    changed = sorted(name for name, expected in baseline.items() if digest(ROOT / name) != expected)
    allowed = {"R/bayes_scale_probabilities.R", "man/diagnose_bayes_probabilities.Rd", "NEWS",
               "tests/testthat/test-bayesian-probabilities.R", "vignettes/bayesian-computation.Rmd"}
    if set(changed) != allowed:
        raise RuntimeError(f"Unexpected existing-source changes: {changed}")
    files = inventory()
    new = sorted(set(files) - set(baseline))
    if new:
        raise RuntimeError(f"Unexpected source additions: {new}")
    for folder in ("prediction_contracts_20261028", "functional_precision_20261028", "spatial_reference_20261028"):
        previous = HERE.parent / folder
        closed = json.loads((previous / "ARTIFACT_MANIFEST.json").read_text(encoding="utf-8-sig"))
        if any(digest(previous / row["file"]) != row["sha256"] for row in closed):
            raise RuntimeError("Closed evidence changed: " + folder)
    workspace = Path(tempfile.mkdtemp(prefix="msr-probability-availability-"))
    package = workspace / "multiScaleR"
    source = {}
    for name in files:
        source[name] = digest(ROOT / name)
        target = package / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(ROOT / name, target)
        if digest(target) != source[name]:
            raise RuntimeError("Snapshot copy mismatch")
    record = {"branch": git("branch", "--show-current"), "head": git("rev-parse", "HEAD"),
              "workspace": str(workspace), "package": str(package), "source_sha256": source,
              "source_files": len(source), "source_bytes": sum((package / n).stat().st_size for n in files),
              "changed_from_baseline": changed, "new_source": new,
              "closed_contract_evidence_unchanged": True}
    temporary = destination.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(record, indent=2), encoding="utf-8")
    temporary.replace(destination)
    print(json.dumps({key: value for key, value in record.items() if key != "source_sha256"}))
