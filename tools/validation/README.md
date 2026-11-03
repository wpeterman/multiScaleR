# Validation phases

New validation batches use `make_plan.py` and `phase_runner.py`. A plan declares
the source commit and hashes, candidate, installed dependencies, arguments,
output directory, libraries, and operational deadlines. The builder copies the
small runtime into a content-addressed directory. Later tooling edits cannot
change a closed or active batch's executable bytes.

Run one declared phase with the frozen runner:

```text
python <frozen-runtime>/phase_runner.py <plan.json> <phase-name>
```

Each attempt reserves a new record before launch. Existing requests, partial
reservations, admissions, logs, processes and results block replay. Inputs and
Git identity are checked before and after the phase. Only explicitly listed Rd
files in an isolated candidate can change during documentation generation.
Declare a new continuation when a phase fails; retain its original records.

The reviewed Windows resource engine inventories actual R interpreters, including
RStudio and Rserve, using executable/module and process-birth identities. Unknown
inventory blocks work. Owned work uses one numerical thread, at most two fit
workers and three R interpreters including the driver. Runtime disk reserve is
5 GiB, launch reserve 20 GiB, and RAM reserve 8 GiB. Normal scheduling applies
only to verified owned interpreters. Foreign jobs are never stopped or reprioritized.

The saved controller may record the user's existing permission to work beside
known foreign R jobs. Without that permission, more than two foreign interpreters
block admission. Unknown inventory, CPU admission, RAM/disk reserves, owned limits
and deadlines remain enforced in either case. Timing comparisons stay ineligible.
Admission refusals are recorded before any R process starts.

Before a fresh package build, require stable 30 GiB free disk. Generate help,
run focused compatibility tests, then build vignettes normally. Pair the normal
`R CMD check --no-manual --no-tests` with the complete installed suite extracted
from that exact archive. Freeze and independently verify source/archive/library/
test identities, scientific markers, warning/skip changes and per-file results
before applying an isolated package candidate.

`seal_phase.ps1` checks only the named phase's recorded process identities. It
does not mistake another phase or foreign worker for unfinished owned work.
Keep compact records rather than full models, chains, tapes or prediction arrays.

Run `python -B tools/validation/test_phase_runner.py` for non-numerical runner
contract tests. Historical benchmark scripts and failure records remain intact;
the reusable runner replaces future script cloning.

`package_batch.py <input.json>` prepares a new independent workspace from a
source manifest, requires four stable disk samples above 30 GiB, and runs the
normal vignette build, normal check and complete exact-archive installed suite.
It uses the same frozen phase runtime and retains each failed attempt. The
archive reader rejects escaping paths, links, source drift and missing vignette
artifacts. Run `python -B tools/validation/test_package_batch.py` for its
non-numerical archive and full-profile tests.

`profiles.json` lists focused feedback filters and the full acceptance policy.
Focused tests never replace the complete suite. Assertion totals describe
coverage; they are not a development target. Remove duplication by sharing
unchanged fixtures and helpers. Keep scientific gates, edge cases and regression
assertions unless a documented replacement preserves their coverage. Current
notes link to archived history, so old validation results retain their exact
source and operational context.
