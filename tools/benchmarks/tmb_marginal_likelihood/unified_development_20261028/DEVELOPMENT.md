# Inference development

Current source: version 0.7.14 explicit-fold GLM cross-validation and compatibility updates.
Read `../publication_cv_20261103/RESUME.md`. Earlier receipts below
describe their own source epochs.

Keep native refitting, analytic likelihood methods, marginal-likelihood research
and Bayesian development on inference-dev. Version 0.7.13 contains the accepted
BAM family, likelihood-bookkeeping, reference-memory, documentation and validation
refinements. [Shared contracts](INFERENCE_DEV_CONTRACTS.md) define the interfaces,
parameter counts and uncertainty distinctions that subsequent work must retain.

## Validation evidence

The accepted full75 archive passed a normal vignette-enabled check and 28,463
installed assertions, with zero failures/warnings and ten unchanged development
helper skips. Archive SHA256:
`cc2cedd01f8535d8f46a94a20403672002662d623ebe567eab69592bc29ecc49`.
All 14 scientific markers and old per-file results were retained, apart from the
reviewed reference/CRS/warning changes and new family coverage. The reusable
validation tools passed 17 non-R contract tests.

The 0.7.13 publication changes only DESCRIPTION version and NEWS release metadata
after that full suite. Its normal vignette-enabled archive check is fresh. R,
generated help, native source, tests and vignette code retain the validated
source identity; exact-archive test files are unchanged. The full-suite result
belongs to the 0.7.12 code epoch and is qualified by this source-equivalence proof.

The working tree additionally contains an accepted, uncommitted NB compatibility
batch. Fresh full78 passed its normal vignette-enabled check and 29,935 installed
assertions, with zero failures/warnings and the same ten skips. All fourteen
scientific markers and every old per-file result were preserved; the new
weighted-NB file adds 1,472 assertions from eight cached software fits. Accepted
archive SHA256:
`1d3e95c268fa666e7995ad4464d723c7550ed4054661a6bf5161b76bd758a776`.
At that checkpoint, the 440 curated source files matched its accepted candidate, and all
fourteen build-ignored legacy example files remain untouched. The reusable
validation tools pass nineteen non-R tests, including three controller lifecycle
scenarios. Records and the current resume instructions are in
`../native_nb_compatibility_20261103`.

The subsequent uncommitted reference-comparison batch passed fresh full79:
Status OK, 30,020 installed assertions, zero failures/warnings and the same ten
skips. All old per-file ledgers and all fourteen scientific markers were
preserved. Its 85 new checks cover exact/tolerance boundaries, nonfinite inputs,
custom comparison dispatch, complete native design verification and fitted-model
preservation. The 441 curated source files match the applied candidate; fourteen
ignored legacy files remain untouched. Archive SHA256:
`bd25d903769ba98acc08abf6772afcc5f6ef103bcbcbb32e12705447cdb56e49`.
That resume records the full79 checkpoint; full80 is the later current source.

Use [tools/validation](../../../validation/README.md) for new frozen, single-use
plans and fresh normal archive/installed-suite validation after behavioral edits.
Preserve every scientific gate and failed/interrupted attempt. Assertion totals
record coverage rather than a growth target. Shared unchanged helpers reduce
duplication without removing regression or failure checks.

The subsequent uncommitted logical-response batch passed fresh full80: normal
vignette check Status OK, 30,726 installed assertions, zero failures/warnings
and the same ten skips. All old per-file results and fourteen scientific
markers remain. Installed224 passes the same 1,184 focused cases without
the two retained source-loading warnings. The 442 curated files match the
applied candidate; fourteen ignored legacy files remain untouched. Archive
SHA256: `d5c7cce4cdba52d81738269c53a54c86f8cbde697d48e767c8a4e0d6905caa68`.

Fresh full81 then passed the binomial selection-copy normalization: normal
vignette check Status OK,31,045 installed assertions, zero failures/warnings
and the same ten skips. All old ledgers and fourteen scientific markers
remain. Focused233 passes1475 checks; the new file adds319. The443 curated
files match the applied candidate. Native fits, covariance, count penalties,
priors, diagnostics and private guards remain unchanged. Archive SHA256:
`3e5a12eb856409b3257c0d79b9dd4602176a1d4b75d4383e95d7d59ec7d1742b`.

Fresh full82 passed the pcount actual-visit comparison correction: normal
vignette check Status OK,31,207 installed assertions, zero failures/warnings
and ten unchanged skips. All old ledgers/fourteen scientific markers remain.
Focused240 passes1637 checks; the new file adds162. Native marginal targets,
variance/mixture counts, covariance, diagnostics, priors and private guards
remain. The444 curated files match the applied candidate. Archive SHA256:
`c7fcd892ffa89c2dc219ae1771baa2f686e2fdb95f121c1467bef4a4b678ca5a`.

Fresh full83 passed the default unmarked frame repair: normal vignette-enabled
check Status OK, 31,283 installed assertions, zero failures/warnings and ten
unchanged skips. Focus244 passes 1,713; the new file adds 76. Native retained
site counts and original IDs now support siteCovs=NULL, including null and
detection-only mixed pcount models. Parameter penalties, priors, covariance,
diagnostics, uncertainty and guards remain unchanged; all old file ledgers
and fourteen scientific markers match. The 445 curated files match archive
SHA256 `99e84388b40705a9115822ff0bd1e35f9e32548d1d50bc01366d068347f69474`.

Fresh full86 passed the first CV stage: normal vignette-enabled check Status
OK,31,545 installed assertions, zero failures/warnings and ten unchanged
skips. The new CV file adds218; the same installed-focused filter passes
2566 checks without warnings. Original source serialization warnings remain
qualified by that installed context. All old per-file ledgers and fourteen
scientific markers match. cross_validate() uses explicit IDs, training-only
scale/coefficient fits, frozen training bins/transforms and retained failed
folds. Existing methods, priors, counts, diagnostics, uncertainty and private
guards are unchanged.449 curated files match archive SHA256
`1ea3e2d7de06cde028a5c2c9262f5f4ead73f448a7a09507ae652a343bc4bc75`.

## Next compatibility work

Known-scale legacy NB and zero/fractional case weights now have point, snapshot,
map, raw-reference curve, worker and selection regression coverage. Large-weight
NB bookkeeping verifies canonical native AIC arithmetic with explicit theta and
rescaled density at the original tolerance, then reports the original weighted
target and true nuisance count. Changed theta state or AIC metadata still fails.
The legacy family's default estimated working scale remains a separate native
convention and is not newly qualified or silently replaced.

Large exact finite comparisons now avoid temporary comparison copies, keeping
complete native prediction batches and the original tolerance fallback for every
difference. The identical-matrix allocation audit used about 82 percent fewer
temporary bytes; it does not measure peak package RAM or uncontended timing.
Further blockwise basis construction remains deferred pending native batch and
equivalence proofs. Logical FALSE/TRUE preparation is now supported with original trial totals,
frames and priors. Equivalent logical/numeric selection signatures are now covered. Mixed
pcount AIC suitability and default null/detection-only frame compatibility
are reviewed in `../unmarked_mixed_aic_20261103/REPORT.md`. The direct CV-stage request is now complete. Next refine explicit-fold
target/comparison contracts and external fold adapters before additional
model-class scores. Other unmarked used-visit checks remain a later
separately declared zero-fit audit.
The later trial/weight audit remains in
`../reference_memory_20261103/GLM_COMPARABILITY.md` before extending GLM compatibility. Numeric 0/1 Bernoulli is already supported using binomial/logit.
Do not reinterpret frequency weights or broaden links without a declared native
audit and derivative/target checks. Private weighted-SQ/spatial fitting admission, normalized
priors and targets, starts/cold/zero/limiting/curvature checks, and unresolved
calibration/reference requirements remain unchanged. Numerical software success
does not promote an experimental backend or establish interval calibration.

Inspect actual global R, source and resources before launches. Use one thread,
at most two owned fit workers and three owned R interpreters including drivers.
Unknown inventory blocks launch. Leave foreign work untouched; known-foreign-R
permission does not supply uncontended timing evidence. Use compact atomic
records and preserve source, results, failures and pinned libraries. Delete only
individually verified unnecessary files.

Local publication/refinement records are retained under
`tools/benchmarks/tmb_marginal_likelihood/publication_20261103` and
`refinement_20261103` when available. Historical protocols and deadlines describe
their own source epochs; they do not override the latest direct user request or
a new declared operational window.

## Version 0.7.14 publication

The accepted full86 code and tests passed 31,545 installed checks, zero failures
or warnings and ten unchanged helper skips. Version 0.7.14 changes only release
metadata after that full suite. Its normal vignette-enabled check is fresh;
its exact installed-focused filter passes 2,566 without warnings. Code, generated
help, native source, tests and vignettes retain the full86 source identity.
Nineteen validation-tool contract tests passed without launching R.

The next authorized batch is grouped/spatial CV folds, guarded CV score tables
and model-selection vignette integration. Preserve parameter counts, priors,
diagnostics, uncertainty and every private guard. Continue curated-only source
copies and compact checkpoints. Keep scientific failures and unique inputs.
Only the present accepted batch is covered by the current commit/push request.
