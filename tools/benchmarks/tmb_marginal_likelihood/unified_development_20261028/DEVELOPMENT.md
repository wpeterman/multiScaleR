# Inference development

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

Use [tools/validation](../../../validation/README.md) for new frozen, single-use
plans and fresh normal archive/installed-suite validation after behavioral edits.
Preserve every scientific gate and failed/interrupted attempt. Assertion totals
record coverage rather than a growth target. Shared unchanged helpers reduce
duplication without removing regression or failure checks.

## Next compatibility work

Inspect the legacy fixed-theta native NB family and case-weighted likelihood
bookkeeping before declaring a small software audit. Verify native family,
scale, covariance and parameter-count conventions before changing an adapter.
Keep all existing uncertainty distinctions and public/native dispatch rules.

Blockwise reference-basis comparison remains deferred pending exact equivalence
and failure tests. Private weighted-SQ/spatial fitting admission, normalized
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
