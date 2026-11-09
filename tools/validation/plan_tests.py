"""Inventory all tests and plan conservative focused feedback, never release acceptance."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import subprocess


def read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def relative(name):
    name = str(name).replace("\\", "/")
    path = PurePosixPath(name)
    if not name or path.is_absolute() or any(p in (".", "..") for p in name.split("/")) or ":" in name:
        raise ValueError("Changed files must use contained repository-relative paths")
    return name


def inventory(repo, profiles):
    repo = Path(repo).resolve()
    known = set(profiles["focused"])
    if any(set(rule["profiles"]) - known for rule in profiles.get("source_profiles", ())):
        raise ValueError("Source map references an unknown profile")
    files = sorted(p.name for p in (repo / "tests/testthat").glob("test-*.R") if p.is_file())
    filters = {name: re.compile(pattern) for name, pattern in profiles["focused"].items()}
    members = {name: [f for f in files if pattern.search(f[5:-2])] for name, pattern in filters.items()}
    missing_smoke = set(profiles.get("smoke_files", ())) - set(files)
    if missing_smoke:
        raise ValueError("Smoke files no longer exist: " + repr(sorted(missing_smoke)))
    mapped = {name: [group for group, chosen in members.items() if name in chosen] for name in files}
    # Source references are an audit starting point, not runtime coverage or
    # proof that a meaningful assertion exercises every exported behavior.
    namespace = (repo / "NAMESPACE").read_text(encoding="utf-8-sig")
    exports = re.findall(r"(?m)^export\(([^)]+)\)$", namespace)
    tests = {name: (repo / "tests/testthat" / name).read_text(encoding="utf-8-sig") for name in files}
    references = {name: [f for f, text in tests.items() if re.search(
        r"(?<![\w.])" + re.escape(name.strip('"')) + r"\s*\(", text) or re.search(
        r"\bdo[.]call\s*\(\s*" + re.escape(name.strip('"')) + r"\s*,", text)] for name in exports}
    methods = []
    for generic, cls in re.findall(r"(?m)^S3method\(([^,]+),([^)]+)\)$", namespace):
        generic, cls = generic.strip(), cls.strip()
        name = generic.split("::")[-1]
        methods.append(dict(generic=generic, cls=cls, method=name + "." + cls,
            generic_source_references=[f for f, text in tests.items() if re.search(
                r"(?<![\w.])" + re.escape(name) + r"\s*\(", text)],
            runtime_dispatch_coverage_established=False))
    return dict(repository=str(repo), profiles_key=hashlib.sha256(json.dumps(profiles, sort_keys=True).encode()).hexdigest(),
        test_files=files, profile_members=members, memberships=mapped,
        unclassified=[name for name, groups in mapped.items() if not groups],
        test_sha256={name: sha(repo / "tests/testthat" / name) for name in files},
        exported_source_references=references,
        registered_S3_methods=methods,
        exports_without_direct_source_reference=[name for name, refs in references.items() if not refs],
        coverage_scope="Static direct/do.call source references, including possible comments; not branch, runtime or assertion coverage.")


def plan(repo, profiles, requested=(), changed=(), full=False, *, catalog=None):
    repo = Path(repo).resolve()
    if catalog is None:
        catalog = inventory(repo, profiles)
    elif catalog.get("repository") != str(repo) or catalog.get("profiles_key") != hashlib.sha256(
            json.dumps(profiles, sort_keys=True).encode()).hexdigest():
        raise ValueError("Reuse only this process's frozen inventory for the same repository and profiles")
    unknown = set(requested) - (set(profiles["focused"]) | {"smoke"})
    if unknown:
        raise ValueError("Unknown profiles: " + repr(sorted(unknown)))
    chosen, reasons, unmapped = set(), [], []
    groups = set(requested) - {"smoke"}
    if "smoke" in requested:
        chosen.update(profiles["smoke_files"])
    changed = sorted({relative(name) for name in changed})
    for name in changed:
        if name.startswith("tests/testthat/test-") and name.endswith(".R"):
            base = PurePosixPath(name).name
            if base in catalog["memberships"]:
                chosen.add(base)
                groups.update(catalog["memberships"][base])
                continue
        matches = [rule for rule in profiles.get("source_profiles", ()) if re.search(rule["pattern"], name)]
        if not matches:
            unmapped.append(name)
        for match in matches:
            groups.update(match["profiles"])
            if match.get("full"):
                full = True
        reasons.append(dict(file=name, profiles=sorted({p for rule in matches for p in rule["profiles"]})))
    if unmapped:
        # Unknown or new helpers/core dependencies never get a narrow automatic
        # plan. A developer can request an explicit profile as partial feedback.
        full = True
    if not set(groups) <= set(profiles["focused"]):
        raise ValueError("Source map references an unknown profile")
    for name in groups:
        chosen.update(catalog["profile_members"][name])
    if full:
        chosen.update(catalog["test_files"])
    if not chosen:
        raise ValueError("Select a profile, changed files or the full inventory")
    names = sorted(chosen)
    # Matches both testthat's stripped names and its test-*.R file names.
    pattern = "(^|test-)(" + "|".join(re.escape(name[5:-2]) for name in names) + ")($|[.]R$)"
    actual = [name for name in catalog["test_files"] if re.search(pattern, name)]
    if actual != names:
        raise ValueError("Focused regex does not select the exact declared files")
    return dict(schema_version=1, selected_files=names, excluded_files=sorted(set(catalog["test_files"]) - chosen),
        selected_profiles=sorted(requested), inferred_profiles=sorted(groups), changed_files=changed,
        unmapped_changed_files=unmapped, reasons=reasons, filter=pattern,
        all_files_selected=names == catalog["test_files"],
        full_release_acceptance_required=True, release_eligible=False,
        scientific_gates_removed=False, test_files_modified=False,
        scope="Developer feedback only; full normal check and exact-archive suite remain the publication gate.",
        catalog=catalog)


def write_new(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    if path.exists() or temporary.exists():
        raise ValueError("Retain the existing plan")
    with temporary.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, indent=2)
    temporary.replace(path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--profile", action="append", default=[])
    parser.add_argument("--changed-files", type=Path, help="JSON array of repository-relative paths")
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    config_path = args.repo / "tools/validation/profiles.json"
    changed = read(args.changed_files) if args.changed_files else []
    if not isinstance(changed, list) or not all(isinstance(x, str) for x in changed):
        raise ValueError("Changed files must be a JSON array of paths")
    value = plan(args.repo, read(config_path), args.profile, changed, args.full)
    value.update(head=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=args.repo).decode().strip(),
        planner_sha256=sha(__file__), profiles_sha256=sha(config_path))
    write_new(args.output, value)
    print(json.dumps(dict(selected=len(value["selected_files"]), total=len(value["catalog"]["test_files"]),
        unclassified=len(value["catalog"]["unclassified"]),
        unreviewed_exports=value["catalog"]["exports_without_direct_source_reference"],
        release_eligible=False, output=str(args.output))))
