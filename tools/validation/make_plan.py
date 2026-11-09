"""Freeze reusable runtime bytes and a structured, single-use validation plan."""
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess

RUNTIME_FILES = ("phase_runner.py", "run_phase.ps1", "resource_control.ps1", "seal_phase.ps1",
                 "document.R", "focus.R", "installed.R", "installed_focus.R", "test_metrics.R", "make_plan.py", "study_inputs.py")

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def make_plan(path, *, repo, output, control, identities, phases, library_paths,
              candidate=None, generated_files=(), sources=()):
    path, repo, output = map(lambda x: Path(x).resolve(), (path, repo, output))
    if path.exists() or not path.is_relative_to(output):
        raise ValueError("Use a new plan inside its output directory")
    code = Path(__file__).resolve().parent
    names = RUNTIME_FILES
    tag = hashlib.sha256(b"".join((code / name).read_bytes() for name in names)).hexdigest()
    runtime = output / "runtime" / tag
    runtime.mkdir(parents=True, exist_ok=True)
    for name in names:
        target = runtime / name
        data = (code / name).read_bytes()
        if target.exists():
            if target.read_bytes() != data:
                raise ValueError("Immutable runtime snapshot changed")
        else:
            with target.open("xb") as stream:
                stream.write(data)
    snapshot_spec = importlib.util.spec_from_file_location("frozen_phase_runtime", runtime / "phase_runner.py")
    snapshot = importlib.util.module_from_spec(snapshot_spec)
    snapshot_spec.loader.exec_module(snapshot)
    base = [dict(group) for group in identities]
    extra = [runtime / name for name in names] + [Path(control).resolve()]
    extra += [Path(name).resolve() for name in sources]
    prepared_phases = []
    for phase in phases:
        phase = dict(phase)
        phase["executable"] = str(Path(phase["executable"]).resolve())
        phase["working_directory"] = str(Path(phase["working_directory"]).resolve())
        arguments = []
        for value in phase["arguments"]:
            if value.startswith("@runtime/"):
                value = runtime / value[len("@runtime/"):]
                if not value.resolve().is_relative_to(runtime) or not value.is_file():
                    raise ValueError("Unknown runtime resource")
                value = str(value)
            arguments.append(value)
        phase["arguments"] = arguments
        phase["argument_line"] = subprocess.list2cmdline(phase["arguments"])
        prepared_phases.append(phase)
        extra.append(Path(phase["executable"]))
    base += [dict(root=str(q.parent), files={q.name: sha(q)}) for q in extra]
    plan = dict(schema_version=1, repo=str(repo), output=str(output), candidate=str(candidate or repo),
        branch="inference-dev", head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo).decode().strip(),
        control_file=str(Path(control).resolve()), identities=base, phases=prepared_phases,
        library_paths=list(library_paths), generated_files=list(map(str, generated_files)),
        frozen_runtime=str(runtime), runtime_sha256=tag)
    snapshot.validate_plan(plan)
    snapshot.verify_git(plan)
    snapshot.verify_inputs(plan)
    snapshot.atomic_json(path, plan)
    return runtime / "phase_runner.py"
