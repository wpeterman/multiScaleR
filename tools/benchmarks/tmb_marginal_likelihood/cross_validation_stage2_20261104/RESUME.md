# Publication checkpoint: version 0.7.15

The user authorized commit and push of this validated CV batch on November 4.
Only DESCRIPTION/NEWS release metadata changed after full32; executable code,
generated help, tests, native source and vignettes retain their accepted identity.
No new numerical work was launched for logoff. The 31,917-assertion full suite and
2,991 compatibility subset cover identical code/tests at the 0.7.14 source epoch;
a fresh version-only 0.7.15 check was not rerun. Read
[NEXT_CV_COMPATIBILITY.md](NEXT_CV_COMPATIBILITY.md) when resuming, and use git log
for the publication commit. A compact local publication receipt records its SHA.
The earlier checkpoint below is retained as historical pre-publication context.

# Completed CV integration on inference-dev

Version 0.7.14 was committed and pushed as
f34f10daaf5d5b83ebeb2b59b2b3dea6ae966252. The remote branch still matches that
commit. The subsequent CV development is validated and applied to the working
checkout, with no additional staging, commit, version bump, merge or push.

## Current behavior

- grouped_cv_folds() keeps whole sampling groups together. spatial_cv_folds()
  keeps whole projected grid blocks, joins blocks connected by shared groups,
  and optionally excludes whole training units within a point-distance buffer.
  Both are deterministic, preserve the RNG and return explicit row-ID splits.
- cross_validate() refits scales and coefficients using training data only.
  cross_validate_glm() supplies canonical GLM null/reference models. Current
  scoring covers Gaussian identity, Poisson log and binomial logit, with native
  equations, factors, interactions, offsets and supported learned transforms.
- cv_tab() ranks compatible pooled or per-fold held-out losses. It matches
  targets, splits, outcomes, trials and frequency weights, retains infinite
  losses and rejects incomplete candidates. It adds no AIC penalty or weights.
- select_cv() evaluates a predeclared training-only AIC/AICc/BIC rule, using
  retained fold fits, native parameter counts and comparability guards. Held-out
  scores never choose that fold's winner. Criterion failures/ties remain failed
  folds. Candidate controls and selected row/map records remain explicit.
- Native grouped raster maps are validated before fitting, subset with their
  buffers and recorded for both parts of each split, including failed folds.
  These fixed covariate maps add no annual/random effects or fitted parameters.
- Both preparation constructors now record fixed x/y cell_size metadata.
  CV retains intersecting cells whose centers exceed the nominal extraction
  radius, using radius plus half the cell diagonal as its fixed geometric bin
  envelope. No observed maximum, clipping or dropped cell is used. The original
  1e-9 gate remains, and geometry is frozen with training bins. Legacy inputs
  lacking geometry retain their strict old gate and request fresh preparation
  when needed. Original extraction, non-CV covariates/bins, parameter bounds,
  priors, diagnostics, coefficient SEs and parameter penalties are unchanged.

Help was generated from roxygen. The new cross-validation vignette contains
worked spatial/grouped folds, ordinary GLM comparisons and training selection;
the main Guide, inference guide and compatibility guide link and explain it.
Model-class scoring beyond the canonical GLMs remains separately qualified work.
Spatial blocking is distinct from the guarded spatially varying scale prototype.

## Exact current validation

Document29: generated help. Focus30: 664 passes, zero failures/errors/warnings/
skips. Vignette31: worked guide rendered, 43,373 bytes. Full32: normal vignette-
enabled build and R CMD check --no-manual --no-tests, Status OK, paired with
its COMPLETE installed suite extracted from that exact archive: 31,917 passes,
zero failures/warnings and ten unchanged external-benchmark helper skips.
Installed33: 2,991 compatibility passes, zero failures/errors/warnings/skips.
The two phases together supply full check/test coverage; tests were not omitted
from acceptance merely because the normal check delegates them to its installed
phase. This is not an as-cran incoming-metadata claim.

Independent verification retains every old per-file ledger, unchanged skip
identity and all fourteen scientific markers. The original starts, cold, zero,
limiting, curvature, memory and private admission guards remain. Numerical
completion is not scientific eligibility, interval calibration or recovery.
Timing comparisons are ineligible because authorized foreign R jobs were active.

Archive SHA256: 78abd17959096e03a67645efbfa4c05ca255b9dfab79b85e7a41875686d5975f.
Read full32-verified-state01.json, installed33-result01.json and applied35.json.
Current 460 curated source files match the accepted archive/source map exactly;
fourteen unchanged legacy inputs are also protected, for 474 root identities.
Applied35 adds fifteen adapter/geometry diffs after the twenty files previously
accepted/applied by full20/applied22. The entire new package batch changes 25
files relative to the pushed commit. No subsequent package-source edit occurred.

## Failures and previous epochs retained

Full20 remains valid for its earlier source: normal check OK, 31,786 complete
installed passes, 2,807 installed-focused passes and the same ten skips/14 markers.
Subsequent edits ended its blanket coverage; full32 supplies the fresh coverage.
Focus09 retained 375 passes, two failures/two errors from test-fixture scope,
sf row IDs and a correctly empty buffered training set. Focus14 retained 479
passes and one incorrect factor-error expectation, then focus17 passed 137.
Vignette15 correctly rejected unreferenced simulated points. Their example now
explicitly declares a local Cartesian meter grid, without weakening CRS guards.
Focus27 retained 556 passes and a real native extraction-envelope error before
fitting. CELL_GEOMETRY_ADDENDUM28.md and focus27-review01.json explain its repair.
All original failure logs, reviews and source copies remain. The pre-geometry
460-file source can be reconstructed byte-for-byte without R using the witnesses
in before-geometry28 and historical26-reconstruction34.json. Earlier changed
fixtures/vignettes are retained in their named before16/before18/before26 copies.
Do not replay any closed producer.

The shared fixture reuses the original seed 9117 data, with exactly 24 scale fits
within its existing ceiling. Native square/rectangular preparation checks add
no optimized fits. No recovery study, sampler/profile panel or passing-seed
search was launched. Priors and tolerances were never relaxed.

## Disk and stopped-work checkpoint

All owned numerical work is complete and sealed. The final inventory found
zero owned R interpreters, one foreign RStudio session and no unknown inventory.
Foreign jobs were never changed. Fresh candidates copied curated source only,
excluding .git, legacy fits and caches. No full fits, chains, tapes or candidate
arrays were saved as checkpoints; results contain compact scores/identities.

Retire36 removed exactly 48 verified regenerable source .o files
(57,167,524 logical bytes, approximately 54.5 MiB).
Source, tests, generated help, DLLs, archives, failures and installed-library
hashes remain unchanged. No directories or unrelated temporary files were
removed. Do not attribute concurrent volume-space changes to this cleanup.
C: currently has approximately 73.2 GiB free. See retire36-plan01.json,
retired36.jsonl and retire36-complete01.json. Keep the stage1 DISK_POLICY.md.

## Next authorized decision / development boundaries

No work is running and no question is pending. Review the current source diff
and guide. A future explicit commit request should bump to 0.7.15 once, complete
the concise NEWS entry (including native cell geometry), validate that release
metadata and stage only the reviewed source/docs/tests and compact notes. Never
stage bulk benchmark evidence or archives. No source changes mean the current
full numerical evidence remains applicable; deliberate edits require qualified
coverage and validation proportional to the change.

Further scoring adapters should prioritize matched-set survival/amt targets and
GAM/mixed/unmarked predictive targets only after native oracles, training-only
refitting and correct observation/group units are declared. Whole-group folds
alone do not supply these likelihoods. Compact native criterion snapshots could
reduce the opt-in memory needed by training selection, while retaining exact
likelihood/K/row-target provenance. Keep all private spatial/ML guards closed
until their original prerequisites and numerical gates pass. Assertion totals
are coverage records, not a growth target; preserve distinct regressions and
scientific gates while sharing fixtures and removing only proven redundancy.
