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

`python -B tools/validation/test_resource_control.py` exercises completed,
running and unknown-inventory process lifecycles with PowerShell process mocks.
It launches no R process. The controller records at least one actual inventory
even when a short phase finishes before the first poll. Existing unknown-state,
owned-process, resource and deadline guards remain in force. This change does
not fill gaps in historical resource logs with invented observations.

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

## Faster development feedback

Use a focused plan while changing a feature. Reserve the complete normal check
and exact-archive suite for a coherent package batch before publication. Repeating
the entire suite after each local edit adds little useful feedback. A focused
pass always lists its selected and excluded files and remains ineligible for
release acceptance.

```text
python -B tools/validation/plan_tests.py --profile smoke --output <new-plan.json>
python -B tools/validation/plan_tests.py --profile gdr_tmb --output <new-plan.json>
python -B tools/validation/plan_tests.py --changed-files <paths.json> --output <new-plan.json>
```

`paths.json` is a JSON array of repository-relative changed paths. Known source
edges combine relevant feature groups. Unknown dependencies, changed shared
helpers and native source select the complete inventory for automatic plans.
Explicit profiles remain available for narrower partial feedback. Inspect the
plan and supply its `filter` to a newly declared `focus.R` phase. The planner
does not start R or bypass the resource controller. Existing focused profile
names continue to work.
Batch callers can build one immutable inventory and pass it to multiple
`plan()` calls in the same process. Do not reread every source for each profile
or reuse an inventory after source changes.

For native PSOCK integration feedback, `installed_focus.R` can run selected
tests from the same frozen candidate after a genuine private installation.
Pin the installation, test files, providers and selected file list in the phase
plan. This avoids R's source-only package-environment serialization warnings;
it does not suppress model warnings or replace normal archive validation.
Its results remain ineligible for release acceptance. The full installed suite
still runs every exact-archive test after the normal vignette-enabled build.

The three-file smoke profile checks coordinate maps, kernel mathematics and
RAM discovery. It is deliberately small and does not qualify native model
fits, predictions or uncertainty. The full profile includes every test file,
including any new file not yet classified. Inventory reports expose unmatched
files, exported-function source references and registered S3 methods. Static
references can include comments and do not prove runtime dispatch, branch
coverage or meaningful assertions. Review those aspects when adding methods.

Future focused and installed phases write a compact `.timings.json` sidecar
beside the unchanged count ledger. It retains elapsed time and outcomes for
each test block and summarizes blocks by file. Setup and unassigned elapsed
time remain separate. Missing timing fields are reported as unavailable.
These are descriptive development measurements, not controlled model speed
benchmarks, child CPU totals or process-tree peak memory estimates. Active and
historical runtime snapshots are not changed or retrospectively profiled.

Use measured slow blocks to find safe fixture reuse. Share immutable numeric
inputs and unchanged read-only fixture preparation where appropriate. Do not
share mutable TMB tapes, replace independent/cold fits with saved likelihoods,
search for passing seeds or hide failed fits. Consolidate duplicate assertions
only after identifying the contract and regression coverage retained by the
replacement. Assertion totals are evidence summaries, not growth targets.

Run `python -B tools/validation/test_plan_tests.py` for planner, inventory,
selection and record-retention unit tests. The existing archive/controller
contract suites remain applicable. Tooling-only changes excluded from the
package archive can be verified separately when all package payload bytes
remain identical; report that scope rather than replaying unchanged numerical
qualification or claiming new mathematical coverage.

The explicit `unmarked_prediction` profile selects the five native unmarked/GDR prediction, trained-recipe and raster guard files. Use it for fast feedback in that area; the broader `native_prediction` profile remains available across model classes. The small profile has the exact file selection exercised by the installed GDR prediction qualification. It does not replace the complete release gate or narrow automatic source-dependency selection.
