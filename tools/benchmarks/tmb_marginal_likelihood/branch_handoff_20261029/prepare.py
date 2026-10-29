"""Stage an explicit dependency-complete source/evidence set and snapshot Git blobs."""
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

ROOT = Path.cwd()
BASE = Path('tools/benchmarks/tmb_marginal_likelihood')
STUDY = BASE / 'branch_handoff_20261029'


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def save(name, value):
    target = STUDY / name
    temporary = target.with_suffix(target.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
    temporary.rename(target)


if (STUDY / 'snapshot01.json').exists():
    raise RuntimeError('Retain existing candidate; use a new numbered attempt.')
assert git('branch', '--show-current').decode().strip() == 'inference-dev'
old = json.loads((BASE / 'annual_sd_profile_20261029/final-source-evidence.json').read_text())
sources = sorted(old['source_sha256'])
allowed = {'.Rbuildignore', 'DESCRIPTION', 'NEWS'}
changed = [p for p in sources if sha(Path(p).read_bytes()) != old['source_sha256'][p]]
assert set(changed) == allowed, changed
paths = set(sources) | {'.gitattributes', (BASE / 'CONTINUE.md').as_posix()}
frozen = {}
for name in ['prediction_contracts_20261028', 'functional_precision_20261028',
             'spatial_reference_20261028', 'probability_precision_20261029',
             'response_sd_profile_20261029', 'spatial_source_contract_20261029',
             'inference_contract_guide_20261029', 'annual_sd_profile_20261029',
             'overnight_checkpoint_20261029']:
    folder = BASE / name
    manifest = folder / ('FINAL_MANIFEST.json' if name.startswith('overnight') else 'ARTIFACT_MANIFEST.json')
    for row in json.loads(manifest.read_text()):
        path = folder / row['file']
        assert sha(path.read_bytes()) == row['sha256'].lower(), str(path)
        frozen[path.as_posix()] = row['sha256'].lower()
        paths.add(path.as_posix())
    paths.add(manifest.as_posix())
unified = BASE / 'unified_development_20261028'
for path in unified.iterdir():
    if path.is_file():
        paths.add(path.as_posix())
for path, expected in json.loads((unified / 'checkpoint.json').read_text())['pinned_sha256'].items():
    if path.startswith('tools/'):
        assert sha(Path(path).read_bytes()) == expected.lower(), path
        paths.add(path)
for name in ['DESIGN.md', 'RESUME.md']:
    paths.add((BASE / 'likelihood_uncertainty_design_20261029' / name).as_posix())
for path in STUDY.iterdir():
    if path.is_file() and path.suffix not in {'.tmp'}:
        paths.add(path.as_posix())
assert all(Path(p).is_file() for p in paths)
assert max(Path(p).stat().st_size for p in paths) < 100 * 1024 ** 2
save('staging-plan02.json', {'paths': sorted(paths), 'package_sources': sources,
                          'frozen_evidence': frozen, 'metadata_changed': changed})
paths.add((STUDY / 'staging-plan02.json').as_posix())
assert set(git('diff', '--cached', '--name-only').decode().splitlines()) <= paths, 'Unexpected index changes'
plan = STUDY / 'explicit-paths.txt'
plan.write_text('\n'.join(sorted(paths)) + '\n', encoding='utf-8')
subprocess.run(['git', '-c', 'core.safecrlf=false', 'add', '--force', '--pathspec-from-file=' + str(plan)], check=True)
whitespace = subprocess.run(['git', 'diff', '--cached', '--check'], capture_output=True)
(STUDY / 'staged-whitespace-review.txt').write_bytes(whitespace.stdout + whitespace.stderr)
tree = git('write-tree').decode().strip()
workspace = Path(tempfile.mkdtemp(prefix='msr-inference-handoff-'))
snapshot = workspace / 'multiScaleR'
hashes = {}
normalized = []
for path in sources + ['.gitattributes']:
    data = git('show', tree + ':' + path)
    target = snapshot / path
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(data)
    hashes[path] = sha(data)
    live = Path(path).read_bytes()
    if live != data:
        assert live.replace(b'\r\n', b'\n') == data.replace(b'\r\n', b'\n'), path
        normalized.append(path)
for path, expected in frozen.items():
    assert sha(git('show', tree + ':' + path)) == expected, path
save('snapshot01.json', {'workspace': workspace.as_posix(), 'snapshot': snapshot.as_posix(),
                         'staged_tree': tree, 'source_sha256': hashes,
                         'source_files': len(hashes), 'line_endings_normalized_only': normalized,
                         'frozen_evidence_git_blobs_verified': len(frozen),
                         'allowed_metadata_changes': changed})
print(json.dumps({'workspace': str(workspace), 'source_files': len(hashes),
                  'frozen_evidence_verified': len(frozen), 'staged_paths': len(paths)}))
