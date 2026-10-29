"""Verify the actual staged candidate, archive, complete tests and compact evidence."""
import hashlib
import json
from pathlib import Path
import re
import subprocess
import tarfile

STUDY = Path(__file__).resolve().parent
ROOT = Path.cwd()
snapshot = json.loads((STUDY / 'snapshot01.json').read_text())
workspace = Path(snapshot['workspace'])


def sha(data):
    return hashlib.sha256(data).hexdigest()


def canonical(path, data):
    text = Path(path).suffix.lower() in {'.r', '.rd', '.cpp', '.h', '.c', '.rmd',
        '.md', '.yaml', '.yml', '.win', '.rproj'} or Path(path).name in {
        '.Rbuildignore', '.gitignore', '.gitattributes', 'Makevars', 'DESCRIPTION',
        'NEWS', 'NAMESPACE', 'CRAN-SUBMISSION'}
    return data.replace(b'\r\n', b'\n') if text else data


def blobs(tree, paths):
    request = ''.join(tree + ':' + p + '\n' for p in paths).encode()
    output = subprocess.check_output(['git', 'cat-file', '--batch'], input=request)
    result = {}; at = 0
    for path in paths:
        end = output.index(b'\n', at)
        header = output[at:end].decode().split()
        assert header[1] == 'blob', header
        size = int(header[2]); start = end + 1
        result[path] = output[start:start + size]
        at = start + size + 1
    assert at == len(output)
    return result


def dcf(data):
    fields = {}; key = None
    for line in data.decode('utf-8').splitlines():
        if not line.strip():
            continue
        if line[0].isspace():
            fields[key] += ' ' + line.strip()
        else:
            key, value = line.split(':', 1)
            fields[key] = value.strip()
    return {k: re.sub(r'\s+', ' ', v.strip()) for k, v in fields.items()}


assert json.loads((STUDY / 'validation01.json').read_text())['passed']
assert 'Status: OK' in (STUDY / 'R-CMD-check02.log').read_text()
match = re.findall(r'\[ FAIL (\d+) \| WARN (\d+) \| SKIP (\d+) \| PASS (\d+) \]',
                   (STUDY / 'installed01.log').read_text(encoding='utf-8'))[-1]
counts = dict(zip(['failures', 'warnings', 'skips', 'passed'], map(int, match)))
assert counts == {'failures': 0, 'warnings': 12, 'skips': 10, 'passed': 5560}, counts
completed_phases = ['documentation01', 'build01', 'check02', 'installed01']
for phase in completed_phases:
    record = json.loads((STUDY / (phase + '-exit.json')).read_text())
    assert record['operational_passed'] and record['exit_code'] == 0, phase
tree = subprocess.check_output(['git', 'write-tree']).decode().strip()
sources = blobs(tree, snapshot['source_sha256'])
for path, expected in snapshot['source_sha256'].items():
    assert sha(sources[path]) == expected, 'Staged source changed: ' + path
    assert canonical(path, Path(path).read_bytes()) == canonical(path, sources[path]), path
    assert canonical(path, (Path(snapshot['snapshot']) / path).read_bytes()) == canonical(path, sources[path]), path
plan = json.loads((STUDY / 'staging-plan02.json').read_text())
evidence = blobs(tree, plan['frozen_evidence'])
for path, expected in plan['frozen_evidence'].items():
    assert sha(evidence[path]) == expected, path
manifests = [p for p in plan['paths'] if Path(p).name in {'ARTIFACT_MANIFEST.json', 'FINAL_MANIFEST.json'}]
old_manifests = blobs(snapshot['staged_tree'], manifests)
current_manifests = blobs(tree, manifests)
for path in manifests:
    assert old_manifests[path] == current_manifests[path] == Path(path).read_bytes(), path
archive = workspace / 'multiScaleR_0.7.5.tar.gz'
archive_sources = 0
with tarfile.open(archive) as stream:
    for member in stream.getmembers():
        if not member.isfile():
            continue
        path = member.name.removeprefix('multiScaleR/')
        assert Path(path).suffix not in {'.o', '.dll', '.so'}, path
        data = stream.extractfile(member).read()
        if path in sources:
            if path == 'DESCRIPTION':
                original, built = dcf(sources[path]), dcf(data)
                assert set(built) - set(original) == {'NeedsCompilation', 'Packaged'}
                assert built['NeedsCompilation'] == 'yes'
                assert all(built[k] == v for k, v in original.items()), path
            else:
                assert canonical(path, data) == canonical(path, sources[path]), path
            archive_sources += 1
        if path.startswith('tests/'):
            assert data == (workspace / 'exact-archive-tests' / path.removeprefix('tests/')).read_bytes(), path
resources = [json.loads(line) for phase in completed_phases
             for p in [STUDY / (phase + '-resources.jsonl')]
             for line in p.read_text().splitlines()]
assert resources and not any(r['unknown'] for r in resources)
assert max(r['owned_interpreters'] for r in resources) <= 4
record = {'passed': True, 'branch': 'inference-dev', 'version': '0.7.5',
    'source_files': len(sources), 'package_source_staged_tree': snapshot['staged_tree'],
    'package_source_unchanged_after_stage': True, 'archive_matched_sources': archive_sources,
    'archive_SHA256': sha(archive.read_bytes()), 'exact_archive_tests_unchanged': True,
    'check_status_OK': True, 'installed_tests': counts,
    'frozen_evidence_git_blobs_verified': len(evidence), 'resource_samples': len(resources),
    'frozen_manifests_git_blobs_verified': len(manifests),
    'maximum_owned_interpreters': max(r['owned_interpreters'] for r in resources),
    'maximum_actual_interpreters': max(r['actual_interpreters'] for r in resources),
    'maximum_foreign_interpreters': max(r['foreign'] for r in resources),
    'minimum_RAM_GiB': min(r['RAM_bytes'] for r in resources) / 1024**3,
    'minimum_disk_GiB': min(r['disk_bytes'] for r in resources) / 1024**3,
    'user_concurrency_exception': True, 'timing_comparison_eligible': False,
    'retained_failed_attempt': json.loads((STUDY / 'check01-resource-stop.json').read_text()),
    'whitespace_review': 'Nine inherited source EOF blanks retained; frozen -text CRLF evidence retained byte-exact. Large local verbose whitespace output is not staged.'}
target = STUDY / 'final-validation01.json'
assert not target.exists(), 'Retain previous final validation.'
temporary = target.with_suffix('.json.tmp')
temporary.write_text(json.dumps(record, indent=2) + '\n', encoding='utf-8')
temporary.rename(target)
print(json.dumps(record, indent=2))
