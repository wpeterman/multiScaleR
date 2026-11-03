"""Run one declared validation phase against immutable inputs, without replay."""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    if path.exists() or temporary.exists():
        raise ValueError(f"Retain existing record: {path}")
    with temporary.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2)
    os.replace(temporary, path)


def inside(path, root):
    return Path(path).resolve().is_relative_to(Path(root).resolve())


def identity_map(group):
    root = Path(group["root"]).resolve()
    for name, expected in group["files"].items():
        path = root / name
        if not inside(path, root) or not re.fullmatch(r"[0-9a-f]{64}", expected):
            raise ValueError(f"Invalid identity entry: {name}")
        yield path, expected


def verify_inputs(plan, after=False):
    generated = {str(Path(name).resolve()) for name in plan.get("generated_files", [])} if after else set()
    for group in plan["identities"]:
        for path, expected in identity_map(group):
            if str(path.resolve()) in generated:
                continue
            if not path.is_file() or digest(path) != expected:
                raise ValueError(f"Input identity changed: {path}")


def validate_plan(plan):
    if plan.get("schema_version") != 1:
        raise ValueError("Unsupported phase-plan schema")
    repo = Path(plan["repo"]).resolve()
    output = Path(plan["output"]).resolve()
    if not inside(output, repo / "tools/benchmarks") or not output.is_dir():
        raise ValueError("Phase outputs must use an existing benchmark directory")
    if plan["branch"] != "inference-dev" or not re.fullmatch(r"[0-9a-f]{40}", plan["head"]):
        raise ValueError("Invalid frozen branch or commit")
    control = read_json(plan["control_file"])
    names = ("long_launch_until", "launch_until", "worker_stop_at", "window_end")
    dates = [dt.datetime.fromisoformat(control[name].replace("Z", "+00:00")) for name in names]
    if any(value.tzinfo is None for value in dates) or dates != sorted(dates):
        raise ValueError("Controller deadlines must be explicit and ordered")
    required = {"threads": 1, "maximum_owned_interpreters": 3, "maximum_owned_fit_workers": 2,
                "RAM_reserve_bytes": 8 * 2**30, "disk_runtime_reserve_bytes": 5 * 2**30,
                "disk_launch_reserve_bytes": 20 * 2**30}
    if any(control.get(name) != value for name, value in required.items()):
        raise ValueError("Resource limits disagree with the reviewed controller")
    if not control.get("foreign_jobs_untouched") or control.get("timing_comparison_eligible") is not False:
        raise ValueError("Preserve foreign work and exclude timing claims")
    phases = plan["phases"]
    if not phases or len({phase["name"] for phase in phases}) != len(phases):
        raise ValueError("Declare unique, nonempty phases")
    for phase in phases:
        if not re.fullmatch(r"[a-z][a-z0-9_-]*", phase["name"]):
            raise ValueError("Unsafe phase name")
        if Path(phase["executable"]).name.lower() not in ("rscript.exe", "rcmd.exe"):
            raise ValueError("Only declared installed R executables are admitted")
        if not isinstance(phase["arguments"], list) or not all(isinstance(x, str) for x in phase["arguments"]):
            raise ValueError("Use structured argument strings")
        if phase["argument_line"] != subprocess.list2cmdline(phase["arguments"]):
            raise ValueError("Windows argument quoting differs from the structured arguments")
        if not Path(phase["working_directory"]).is_dir():
            raise ValueError("Phase working directory is unavailable")
    paths = {str(path.resolve()) for group in plan["identities"] for path, _ in identity_map(group)}
    expected = {str(Path(plan["control_file"]).resolve()), str(Path(__file__).resolve()),
                str(Path(__file__).with_name("run_phase.ps1").resolve()),
                str(Path(__file__).with_name("resource_control.ps1").resolve())}
    expected.update(str(Path(phase["executable"]).resolve()) for phase in phases)
    if not expected.issubset(paths):
        raise ValueError("Controller, runner and executable identities must be frozen")
    generated = {str(Path(x).resolve()) for x in plan.get("generated_files", [])}
    if not generated.issubset(paths):
        raise ValueError("Generated documentation must have a declared starting identity")
    for name in generated:
        if Path(name).suffix.lower() != ".rd" or not inside(name, Path(plan["candidate"]) / "man"):
            raise ValueError("Only isolated candidate Rd files may be regenerated")
    return plan


def verify_git(plan):
    cwd = plan["repo"]
    call = lambda *args: subprocess.check_output(["git", *args], cwd=cwd).decode().strip()
    if call("branch", "--show-current") != plan["branch"] or call("rev-parse", "HEAD") != plan["head"]:
        raise ValueError("Live Git identity changed")
    if call("diff", "--cached", "--name-only"):
        raise ValueError("Index changed; preserve unrelated staging")


def run(plan_path, phase_name):
    plan_path = Path(plan_path).resolve()
    plan = validate_plan(read_json(plan_path))
    phase = next((x for x in plan["phases"] if x["name"] == phase_name), None)
    if phase is None:
        raise ValueError("Phase was not declared")
    output = Path(plan["output"])
    # Reserve before launching; failed and interrupted attempts remain closed.
    request = output / (phase_name + "-request.json")
    markers = [phase_name + "-" + name + ".json" for name in
               ("request", "admission", "admission-blocked", "process", "exit", "verified")]
    markers += [phase_name + "-request.json.tmp", phase_name + "-resources.jsonl",
                phase_name + ".log", phase_name + ".err"]
    if any((output / name).exists() for name in markers):
        raise ValueError("Retain prior attempt; declare a new phase for a continuation")
    verify_git(plan)
    verify_inputs(plan)
    atomic_json(request, dict(plan=str(plan_path), plan_sha256=digest(plan_path), phase=phase_name,
                             requested_at=dt.datetime.now().astimezone().isoformat()))
    wrapper = Path(__file__).with_name("run_phase.ps1")
    result = subprocess.run(["pwsh", "-NoProfile", "-File", str(wrapper), str(plan_path), phase_name], cwd=plan["repo"])
    if result.returncode:
        raise RuntimeError(f"Retain failed phase {phase_name}, exit {result.returncode}")
    verify_git(plan)
    verify_inputs(plan, after=True)
    generated = {name: digest(name) for name in plan.get("generated_files", [])}
    atomic_json(output / (phase_name + "-verified.json"),
                dict(plan_sha256=digest(plan_path), phase=phase_name, inputs_verified=True,
                     generated_files=generated, source_edit_admission=False,
                     verified_at=dt.datetime.now().astimezone().isoformat()))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("plan", type=Path)
    parser.add_argument("phase")
    args = parser.parse_args()
    try:
        run(args.plan, args.phase)
    except (ValueError, RuntimeError, OSError, subprocess.SubprocessError) as error:
        print(str(error), file=sys.stderr)
        sys.exit(1)
