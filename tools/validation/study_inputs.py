"""Current source/dependency identities for one declared inference refinement study."""
import json
from pathlib import Path


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def study_inputs(repo, study, candidate_files=None):
    repo, study = Path(repo).resolve(), Path(study).resolve()
    baseline = read(study / "baseline01.json")
    feature = study.parent / "feature_compatibility_20261029"
    normal = read(feature / "gate74/normal-check-complete01.json")
    pin = read(feature / "gate72/dependency-pin02.json")
    identities = [dict(root=str(repo), files=baseline["root_source_sha256"]),
                  dict(root=normal["installed_library"], files=normal["installed_code_sha256"]),
                  dict(root=normal["exact_archive_tests"], files=normal["test_sha256"]),
                  dict(root=str(Path(pin["overlay"]) / "MASS"), files=pin["files"])]
    candidate = read(study / "candidate192-initial01.json")["package"]
    if candidate_files is not None:
        identities.append(dict(root=candidate, files=candidate_files))
    libraries = [str(Path(normal["installed_library"]).parent), pin["overlay"],
                 str(repo / "tools/benchmarks/joint_scale_assessment/phase31/library"),
                 "C:/Users/peterman.73/.codex/scratch/multiScaleR-bayes-library",
                 "C:/Users/peterman.73/AppData/Local/R/win-library/4.6",
                 "C:/Program Files/R/R-4.6.1/library"]
    metadata = [study / "baseline01.json", feature / "gate74/normal-check-complete01.json",
                feature / "gate72/dependency-pin02.json"]
    return dict(identities=identities, library_paths=libraries, candidate=candidate,
                sources=metadata, dependency_overlay=pin["overlay"])
