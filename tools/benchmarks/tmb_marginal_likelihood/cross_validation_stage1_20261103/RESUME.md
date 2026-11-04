# Explicit-fold GLM cross-validation: accepted and applied

Updated 2026-11-03T23:03:24.568524-05:00. Branch inference-dev, HEAD c5dab380, version0.7.13.
Eleven independently accepted files are applied without staging, commit, push,
merge or version change. All earlier dirty work is preserved. The463 recorded
root identities match449 curated source files plus fourteen unchanged
build-ignored legacy files. Read applied86.json and full86-verified-state01.json.

## Resulting behavior

cross_validate(formula, data, kernel_inputs, folds, ...) refits scales, shapes
where applicable and GLM coefficients on each explicit training split. Gaussian
identity, Poisson log and binomial logit use the unchanged profile_glm or refit
backends. Native formula handling retains factors, interactions, polynomial
terms and offsets subject to the existing backend's guards. Frequency weights
and binomial trial counts remain separate. No custom family/link, mixed, GAM,
unmarked or Bayesian scoring adapter is implied by this first stage.

Each fold supplies disjoint analysis/assessment character row IDs. Rows may be
unassessed for buffers or assessed repeatedly; assessment_counts records every
occurrence. The user defines ecological sampling units and dependence-aware
folds. Automatic random, grouped, spatial or temporal fold generation is a later
stage. Use the same planned targets/folds/weights across candidate comparisons.

Raw per-point buffers are reused read-only. Initial spatial normalization and
pooled representative bin distances are learned only from training rows, over
the predeclared physical buffer domain. Assessment values use those frozen bins
and the fitted training transformations. Simply slicing existing full-data bins
would not satisfy this contract because their representatives pool all sites.
Stored cell-level data are therefore required; lean-only and grouped preparation
are outside this stage's qualified input scope. At least one optimized spatial
predictor is required. No full-data model is fitted or used as a warm start.

Gaussian scores are RMSE/MAE; Poisson adds negative log score and deviance;
binomial uses negative log score and Brier score. Negative log scores retain
native count constants and frequency replication. Binomial Brier loss averages
over trials, with frequency times trial count as its aggregation weight. Zero
weight/trial rows remain recorded and contribute no score. An observed outcome
assigned zero probability has an infinite log score, without artificial clipping.
RMSE pools squared errors before taking the square root.

The result keeps every planned assessment occurrence, scores, scales, training
normalization, bin references and original diagnostics. Full fold fits are opt-in
in memory; compact records are the default. Warnings/errors are retained per
fold. Any failed fold makes pooled estimates incomplete/NA, while successful
partial results remain inspectable. No failed fold is dropped or retried.
complete describes calculation completion, not diagnostic eligibility or
calibration. These are plug-in predictions/scores with no additional parameter
integration or observation intervals. Existing priors, covariance/SE targets,
parameter penalties, starts/cold/convergence/curvature checks and every private
admission guard remain unchanged. No new dependency was introduced.

## Exact validation

Fresh full86 normal vignette-enabled package check: Status OK. The full
exact-archive installed suite passes31,545 assertions, with zero failures or
warnings and ten unchanged helper skips. Every old per-file result and all
fourteen scientific markers match full83. The new CV file adds218 assertions.
The same installed-focused filter37 passes 2,566 checks with zero failures,
errors, warnings or skips. It has identical per-file pass counts to source
focus36. Two existing source/pkgload serialization warnings at
test-profile-weighted-surface.R:113 are retained; the real installed namespace
resolves that context with zero warnings. No warning or numerical gate was
suppressed or relaxed. Roxygen document08 generated and explicitly verified
the two namespace entries and new help file. Later source/test corrections
leave its documented interface and generated help bytes unchanged.

The current root is byte-equivalent to the independently accepted archive,
SHA256 1ea3e2d7de06cde028a5c2c9262f5f4ead73f448a7a09507ae652a343bc4bc75.

The one seed9117 software fixture plans18 fold fits within a ceiling24 per
lifecycle and caches fits/failures once. Native fixed-scale GLMs independently
check predictions, coefficients and score arithmetic. Cases cover all three
families, logical/grouped responses, frequency weights, factors, interactions,
poly/offset terms, Gaussian/exponential/expow kernels, weighted SQ, native refit,
two spatial predictors, held-out perturbations, bin references, row IDs and
intentional prediction failure from an unseen level. Changing assessment
outcomes/predictors/buffers or full-preparation summaries leaves every training
mathematical quantity and diagnostic unchanged. Timing columns remain preserved
but are not equality targets or evidence of efficiency. These are software
regressions, not parameter-recovery or predictive-calibration claims.

## Retained failures and qualifications

document03 was refused before R launch by the generic existing-Rd output rule;
the separate document08 producer pins its own code, checks the old namespace
hash and exactly two generated entries, and leaves all other inputs frozen.
document05 was refused before R launch for disk reserve. The user freed space;
the subsequent inventory verified recovery. No shared resource/source guard was
changed. focus09 found a single-column vapply shape bug in preparation before
any CV scale fit. R12 corrects only the output matrix shape. focus13 completed
the CV fitting cases but exposed two test-reference errors: formula-environment
weight lookup and comparison of elapsed seconds. Test15 passes oracle weights
literally and compares all mathematical fields while retaining timing records.
All source, failed logs/results and original epochs remain. No passing-seed
search, optimizer retuning or diagnostic relaxation occurred.

Full84 completed the installed suite with31,500 passes, one failure, zero
warnings and ten unchanged skips. Its sole failure was the existing
expected-public-exports list omitting the new cross_validate entry. The
separate candidate19 adds that one expected export and keeps the exact
set-equality assertion. Production source, generated help, priors and
every numerical criterion are unchanged. Full84 and its original candidate
are retained; installed17 never launched. Fresh focus36/full86/installed37
qualify the corrected test contract.

Scores must match the intended prediction task. Model/kernel/metric selection
using CV needs outer validation or an independent final test set to assess the
selected procedure. Site/animal/group independence and spatial transfer are not
inferred from row IDs. Proper scoring of mixed, unmarked or Bayesian predictive
distributions needs separate target and integration qualification; setting
random effects to zero is not marginalization. Private fitting admission,
spatial recovery limitations and all earlier numerical failures remain intact.
Historical as-cran incoming metadata timeouts from full67 remain unresolved
external checks; this normal vignette-enabled package validation is complete.

## RAM discovery correction and disk use

Focus20 retained2110 passes, zero assertion failures, two errors and the same
two source/pkgload warnings. CV and the export test passed. Two existing
Bayesian tests stopped at their unchanged unknown-RAM guard before tapes:
a transient Windows CIM discovery failure had been cached as NA. Candidate35
uses a direct Windows .NET memory query, retains the CIM fallback and caches
only finite positive readings. Unavailable readings are retried later. Both
independent Windows queries and the actual R discovery measured68355309568
bytes. Document35 verified that total with no fit or tape. No guessed RAM
value, threshold relaxation, new dependency or prior/diagnostic change was
introduced. Forty-four new tests check direct discovery, fallback, invalid
output, cache recovery and the original unknown/insufficient-memory stops.
The initial source failure and original candidate remain retained; fresh
focus36/full86/installed37 qualify the corrected source.

Earlier user-authorized duplicate cleanup removed2254 exact legacy-fit
copies from131 old temporary trees (16.64GiB logical content). All fourteen
canonical originals and unique scientific/source/log evidence remain. Read
disk_audit24/REPORT.md, its deletion receipt and restoration manifest. Old
temporary manifests need those redundant inputs restored before reuse;
closed producers stay closed. Current source and failed CV candidates are
byte-verified. New copies use curated manifests and exclude those fits, .git
and caches. Size ledgers cover explicitly owned paths; no whole-drive polling
or full numerical arrays are checkpointed. DISK_POLICY.md records the
continuing limits. Logical file sizes are not net physical reclaim.
After all owned phases closed and acceptance/application passed, retirement41
removed13 generated source-tree compilation objects
(24.6MiB), excluded from every frozen source
identity. All449 source identities and the actual installed library remained
byte-intact. Source and installed DLLs are distinct; the accepted installed DLL
is retained. Failed candidates are untouched. The final explicit-path ledger
totals170.0MiB logical and
170.0MiB allocated across current/retained
CV candidates, the failed/current check workspaces and compact outputs.
Cross-root hardlinks are not deduplicated in that total. This ledger is
scoped to this batch, not a whole-drive inventory.

## Resource checkpoint and next work

No owned R remains. Known foreign jobs are untouched; unknown inventory stayed
zero. Peak owned R was3; minimum free RAM
25.50GiB and disk 74.39GiB.
One numerical thread and frozen resource controls were enforced. Timing
comparisons are ineligible. Only verified duplicates and expendable compiled source objects were removed; no full fit/chain/tape arrays
were checkpointed and no recurring automation was created.

The first CV stage is complete. Next supported refinements are explicit-fold
target/comparison contracts, adapters for external fold definitions, and then
qualified additional model-class scoring. Pointwise immutable bin sufficient
statistics could reduce repeated preparation, but require training-only
representative-distance and memory/equivalence proofs before optimization.
Do not silently broaden family/link/weight or private-model scope. Preserve all
old priors, guards, diagnostics, counts and uncertainty distinctions. Any later
source edit qualifies this historical coverage and needs fresh validation.
Do not replay closed producers or run the pre-application verifier against
the applied root. No staging/commit/push/version bump is authorized.
