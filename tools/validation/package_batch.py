"""Run one frozen candidate through normal build/check and its exact installed tests."""
import argparse
import datetime as dt
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tarfile
import time

from make_plan import make_plan, sha


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def write(path, value):
    path = Path(path)
    if path.exists() or path.with_suffix(path.suffix + ".tmp").exists():
        raise ValueError("Retain the existing batch record")
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
    temporary.replace(path)


def archive_tests(archive, candidate, sources, destination):
    """Extract regular test members only, after validating every archive path."""
    members, normalized = {}, {}
    destination = Path(destination).resolve()
    if destination.exists():
        raise ValueError("Exact-archive test destination must be new")
    with tarfile.open(archive, "r:gz") as stream:
        checked = []
        for item in stream.getmembers():
            parts = item.name.split("/")
            if parts[0] != "multiScaleR" or any(x in (".", "..") or ":" in x or "\\" in x for x in parts):
                raise ValueError("Unexpected archive path")
            if item.isdir():
                continue
            if not item.isfile():
                raise ValueError("Unexpected non-file archive member")
            relative = "/".join(parts[1:])
            if not relative or relative in members:
                raise ValueError("Duplicate or empty archive member")
            data = stream.extractfile(item).read()
            members[relative] = hashlib.sha256(data).hexdigest()
            if parts[1] == "tests":
                checked.append((Path(*parts[2:]), data))
        for name, expected in sources.items():
            path = Path(name)
            if path.parts[0] in ("R", "man", "tests", "src") or (path.parts[0] == "vignettes" and path.suffix == ".Rmd"):
                if members.get(name) != expected:
                    native = hashlib.sha256((Path(candidate) / name).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
                    if path.parts[0] != "src" or members.get(name) != native:
                        raise ValueError("Archive source mismatch: " + name)
                    normalized[name] = dict(source_sha256=expected, archive_sha256=native, native_LF_only=True)
        for name in ("build/vignette.rds", "inst/doc/feature-compatibility.html"):
            if name not in members:
                raise ValueError("Missing normal-build vignette artifact: " + name)
        destination.mkdir()
        for name, data in checked:
            target = (destination / name).resolve()
            if not target.is_relative_to(destination):
                raise ValueError("Escaping test member")
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(data)
    if not (destination / "testthat.R").is_file():
        raise ValueError("Missing exact-archive test entry point")
    return dict(archive=str(archive), archive_sha256=sha(archive), tests=str(destination), members=members,
                source_identity_passed=True, normal_vignette_index_present=True, native_line_endings=normalized)


def run(input_path):
    input_path = Path(input_path).resolve()
    config = read(input_path)
    repo, output, package = [Path(config[name]).resolve() for name in ("repo", "output", "candidate")]
    workspace = Path(config["workspace"]).resolve()
    if (not output.is_relative_to(repo / "tools/benchmarks") or output == repo / "tools/benchmarks"
            or workspace.parent != Path(os.environ["TEMP"]).resolve()
            or not workspace.name.startswith("msr-") or workspace.exists() or not package.is_dir()):
        raise ValueError("Require a new owned temporary workspace and contained benchmark output")
    output.mkdir(exist_ok=False)
    sources = read(config["source_manifest"])["source_sha256"]
    if any(sha(package / name) != value for name, value in sources.items()):
        raise ValueError("Candidate changed before full validation")
    workspace.mkdir(exist_ok=False)
    copied = workspace / "multiScaleR"
    for name, value in sources.items():
        target = (copied / name).resolve()
        if not target.is_relative_to(copied):
            raise ValueError("Escaping candidate source")
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(package / name, target)
        if sha(target) != value:
            raise ValueError("Copied source mismatch")
    identities = config["identities"] + [dict(root=str(copied), files=sources)]
    libraries = config["library_paths"]
    metadata = config["sources"] + [str(input_path), str(Path(__file__).resolve())]
    baseline = dict(head=config["head"], workspace=str(workspace), candidate_package=str(package),
                    source_sha256=sources, source_files=len(sources), root_source_sha256=config["root_source_sha256"])
    write(output / "baseline01.json", baseline)

    def phase(name, executable, arguments, extra=(), installed=False):
        plan = output / (name + "-plan01.json")
        runtime = make_plan(plan, repo=repo, output=output, control=config["control_file"],
            identities=identities + list(extra), phases=[dict(name=name, executable=executable,
                arguments=arguments, working_directory=str(workspace), long=True, installed_suite=installed)],
            library_paths=libraries, candidate=package, sources=metadata)
        subprocess.run([sys.executable, "-B", str(runtime), str(plan), name], cwd=repo, check=True)
        subprocess.run(["pwsh", "-NoProfile", "-File", str(runtime.with_name("seal_phase.ps1")),
            "-Output", str(output), "-Phase", name, "-TargetName", name + "-owned-seal01.json"], cwd=repo, check=True)

    samples = []
    for i in range(4):
        free = shutil.disk_usage(workspace).free
        samples.append(dict(time=dt.datetime.now().astimezone().isoformat(), disk_bytes=free))
        write(output / ("build-stable-sample%d.json" % (i + 1)), samples[-1])
        if free <= 30 * 2**30:
            raise ValueError("Require stable 30 GiB free disk before build")
        if i < 3:
            time.sleep(15)
    rcmd = config["rcmd"]
    phase("build01", rcmd, ["build", "multiScaleR"])
    version = re.search(r"(?m)^Version: ([^\r\n]+)$", (package / "DESCRIPTION").read_text())[1]
    archive = workspace / ("multiScaleR_" + version + ".tar.gz")
    built = archive_tests(archive, package, sources, workspace / "exact-archive-tests")
    write(output / "archive01.json", built)
    # Include the archive identity in both the normal check and complete suite.
    identities.append(dict(root=str(workspace), files={archive.name: built["archive_sha256"]}))
    phase("check01", rcmd, ["check", archive.name, "--no-manual", "--no-tests"])
    log = (workspace / "multiScaleR.Rcheck/00check.log").read_text(encoding="utf-8-sig")
    (output / "R-CMD-check.log").write_text(log, encoding="utf-8")
    status = re.findall(r"(?m)^Status: (.*)$", log)[-1].strip()
    notes = [line for line in log.splitlines() if line.startswith("* checking") and line.endswith(" NOTE")]
    accepted = status == "OK" or (re.fullmatch(r"[12] NOTEs?", status) and notes and all(line.startswith("* checking examples") for line in notes))
    write(output / "check-assessment01.json", dict(status=status, notes=notes, accepted=bool(accepted), timing_comparison_eligible=False))
    if not accepted:
        raise ValueError("Retain normal-check findings before any continuation")
    installed = workspace / "multiScaleR.Rcheck/multiScaleR"
    installed_files = {q.relative_to(installed).as_posix(): sha(q) for group in ("R", "libs", "Meta") for q in (installed / group).rglob("*") if q.is_file()}
    installed_files.update({name: sha(installed / name) for name in ("DESCRIPTION", "NAMESPACE")})
    tests = Path(built["tests"])
    test_files = {q.relative_to(tests).as_posix(): sha(q) for q in tests.rglob("*") if q.is_file()}
    normal = dict(normal_check_completed=True, tests_complete=False, check_status=status, normal_vignette_enabled=True,
        installed_library=str(installed), installed_code_sha256=installed_files, exact_archive_tests=str(tests),
        test_sha256=test_files, archive=str(archive), archive_sha256=built["archive_sha256"])
    write(output / "normal-check-complete01.json", normal)
    extra = [dict(root=str(installed), files=installed_files), dict(root=str(tests), files=test_files)]
    inp = output / "installed-input01.json"
    write(inp, dict(installed_library=str(installed), exact_archive_tests=str(tests),
        dependency_overlay=config["dependency_overlay"], result=str(output / "test-counts01.json")))
    metadata.append(str(inp))
    libraries = [str(installed.parent)] + libraries
    phase("installed01", config["rscript"], ["--vanilla", "@runtime/installed.R", str(inp)], extra=extra, installed=True)
    write(output / "complete01.json", dict(passed=True, tests_complete=True, check_accepted=True, check_status=status,
        normal_vignette_enabled=True, tests_from_exact_archive=True, archive_sha256=built["archive_sha256"],
        source_files=len(sources), completed_at=dt.datetime.now().astimezone().isoformat(), timing_comparison_eligible=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    run(parser.parse_args().input)
