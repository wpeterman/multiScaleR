# Next bounded batch: functional probability precision

Do not edit package source while check04/installed04 is running. Read RESUME.md
for the current external validation phases and retained earlier attempts. First finish the exact
prediction-contract candidate, inspect its installed tests/check and audit live
source identity. Then freeze a separate small development batch on the same
branch. The user authorized development only through04:45Eastern October29.

The next scientific priority is a diagnostic for probabilities such as
P(kernel distance <=250map units), rather than eligibility inferred from aggregate
Rhat/ESS. No new sampler study is needed to implement and test the diagnostic.
The existing diagnose_bayes_scales interface and all original gates remain intact.

Suggested addition: a separate diagnostic function accepting an existing Bayesian
fit, a named list of kernel-specific finite thresholds, their distance/log axis,
and an optional declared probability-MCSE target. Validate requests before
extracting any chains. Reuse original unpermuted post-warmup chains and the
original parameter map. Never use permuted/thinned draws or importance candidates
as chains. Compact and approximate objects must report explicit unavailability.

For each requested predictor, threshold, and distinct active fitted scale row,
reconstruct log distance using the frozen scale design and original grouping/SD
indices. Compare on the log axis with a documented <=event. Compute the event
probability, per-chain probabilities, indicator mean-ESS, MCSE and Rhat with the
installed posterior implementation. Keep only compact summaries. Do not serialize
full arrays, draw populations, tapes or samplers. Preserve caller RNG and fit.

Constant indicators are nonestimable for precision, even if their empirical
probability is0or1. Nonfinite values, missing/duplicate parameter columns, stuck
chains and single-chain limitations must stay explicit. A user-selected precision
target is a separate report, not a replacement for raw or derived fit diagnostics.
Do not return a new general fit-passed flag or silently change thresholds.
An empirical zero event rate is not proof of an absent tail or unvisited mode.

Use existing synthetic unpermuted-chain fixtures for deterministic tests of exact
probabilities, chain ordering, distance conversion, fixed/random/crossed designs,
scalar/bracketed parameter names, constant/stuck/nonfinite cases, saved compact
unavailability, request validation and RNG/object immutability. Add an independent
ordinary-array reconstruction and mean-MCSE reference; comparisons must not merely
mirror helper code. No NUTS fitting, tuning or new calibration panel in this batch.

Read tests/testthat/test-bayesian-derived-diagnostics.R and
R/bayes_scale_diagnostics.R first. The closed consolidation study's long-reference.R
contains retained functional-CDF diagnostics; inspect it read-only and reuse the
numerical lesson without rerunning its producer. Correlated private spatial blocks
need their own projection-aware adapter; do not coerce them into the public class
or send them into scalar factorization. Before any later spatial reference fit,
declare a fresh bounded protocol with fixed target/priors/seeds/controls/gates.

Generate documentation through roxygen, explain Monte Carlo versus posterior
uncertainty with a fast example, and run focused tests before further integration.
A new export requires a deliberate strict export-inventory test update. Preserve
the fullcheck evidence for the previous candidate and qualify subsequent drift.
