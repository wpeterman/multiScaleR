# Assertion coverage during integration

Assertion totals describe a particular run. They are not a development target,
a measure of scientific validity, or a reason to retain duplicate checks.
Coverage should grow when a feature introduces another failure mode. It need
not grow when an existing native equality already verifies the same behavior.

Keep a small set of distinct checks for each applicable contract:

| Contract | Evidence that must remain distinct |
| --- | --- |
| Fitted target | Unchanged coefficients, covariance, likelihood, priors, reference transforms and diagnostics |
| Native output | Values, dimensions, labels, ordering, component spans and family-specific quantities |
| Data context | Full, reordered, unequal and single-row inputs; trained factor/logical coding; fixed learned transformations and offsets |
| Public route | Point, numerical snapshot, checked raster, effect curve and model-selection behavior |
| Uncertainty | Native coefficient covariance, optional native smoothing covariance, ecological-scale uncertainty, posterior uncertainty and observation noise |
| Model selection | Native likelihood/EDF metadata plus estimated ecological parameters; fixed radii and unestimated nuisance terms must not be counted again |
| Scientific admission | Original starts, cold/reference, zero/limiting, normalized-density, curvature and diagnostic gates; retained failures remain evidence |
| Provenance | Actual native provider, source/archive identity, fitted rows, physical scales, normalization and installed code |

A native equality can usually check a whole returned object, including its
attributes, instead of separate near-duplicate assertions for every field.
Keep a separate assertion when it has a different oracle or interpretation.
For example, a native value comparison does not replace an independent
coefficient derivative check, and a conditional SE does not validate joint
scale uncertainty. A class name or finite result does not establish eligibility.

Before pruning, map each proposed deletion to the surviving check and the
failure it would detect. Verify that the surviving oracle is independent and
that all relevant routes, boundary cases and uncertainty distinctions remain.
Remove duplicate expectations together in a reviewable maintenance batch,
with an explicit coverage explanation and a fresh full validation. Do not
delete a failing test, increase a tolerance, remove an original numerical
gate or discard a failed fit to reduce a count. No pruning is part of the
current integration batch.

Use focused tests for each edit and the complete exact-archive suite for
coherent batches. The current per-file ledgers compare an unchanged baseline
to catch accidental coverage loss; new expectations are derived from the
actual new test file. They are not quotas. Explain intentional count changes
when a later coverage review consolidates equivalent checks. A baseline from
another source, provider version or harness must be qualified before comparing
totals or timings.

Keep interrupted validation separate from completed validation. A clean normal
vignette check is useful evidence, but an interrupted archive suite has no
complete assertion ledger. Preserve its outputs and resume only unfinished
validation under a new phase identity. External CRAN incoming failures and
resource interruptions do not become successful checks or scientific failures.

The active source/status is recorded in
`../continuation_20261030/RESUME.md` and `current-state.json`. Those records
distinguish applied, isolated, validated, incomplete and unlaunched work.
