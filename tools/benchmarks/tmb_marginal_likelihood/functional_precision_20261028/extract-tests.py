"""Extract only the exact built archive's tests into a new owned working folder."""
import hashlib
import json
from pathlib import Path
import sys
import tarfile

archive, destination = map(Path, sys.argv[1:])
destination = destination.resolve()
destination.mkdir(exist_ok=False)
rows = []
with tarfile.open(archive, "r:gz") as stream:
    for member in stream.getmembers():
        if not member.name.startswith("multiScaleR/tests/") or member.isdir():
            continue
        if not member.isfile():
            raise ValueError("Unexpected non-file archive test member")
        relative = Path(member.name.removeprefix("multiScaleR/tests/"))
        target = (destination / relative).resolve()
        if not target.is_relative_to(destination):
            raise ValueError("Archive member escapes the test workspace")
        data = stream.extractfile(member).read()
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(data)
        rows.append({"file": relative.as_posix(), "sha256": hashlib.sha256(data).hexdigest()})
if not (destination / "testthat.R").is_file():
    raise ValueError("Archive test entry point is missing")
(destination / "archive-test-manifest.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")
print(json.dumps({"tests_extracted": len(rows), "archive": str(archive)}))
