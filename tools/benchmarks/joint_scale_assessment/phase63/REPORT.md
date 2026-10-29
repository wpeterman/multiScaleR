# Spatial prediction checkpoint, October 28, 2026

Development is authorized through 15:30 America/New_York today. No long work
after 14:45; final checkpoint starts at 15:10. Today's one-day automation is
`multiscaler-spatial-workflows-until-3-30pm`; delete it before logoff. All changes
remain uncommitted. Defaults, independent starts, Bayesian/owner changes, native
audits, explicit failures, Hessians and uncertainty restrictions are preserved.

## Integrated spatial workflows

- `prepare_unmarked_data(..., include_fixed = TRUE)` and the corresponding
  `prepare_prediction()` option reconstruct modeled fixed-radius predictors
  only after verifying their source/metric/radius definitions, original
  normalization and exact agreement with retained native training inputs.
  Explicit site keys preserve prediction and repeated-scenario row order.
  Default `FALSE` keeps fixed predictors as explicit native inputs.
- `kernel_scale.raster(..., include_fixed = TRUE, scale_center = TRUE)` uses
  the same verified normalization for unmarked, GLM and coxph fits. GLM/coxph
  require retained native model frames, numeric original predictor columns and
  original prepared row IDs. Opaque/transformed frames raise errors rather than
  guessing. Only a local object copy receives the additional scaling parameters.
- Native unmarked predictor discovery now includes modeled fixed metrics.
  The prior omission caused raster prediction to fail with `rough` missing.
  Clamping considers spatial numeric columns without coercion by site keys or
  unrelated factors. Nonzero `pct_mx` behavior was not changed.
- `kernel_scale.raster(..., match_buffer = TRUE)` matches the original prepared
  polygon-intersection footprint for kernel means and smooth weighted surface
  metrics at raster cell centers. Partially intersecting cells receive the
  original full-cell kernel weight. Finite mapped support is renormalized;
  missing central source cells remain masked. The default `FALSE` retains the
  historical `pct_wt` window and its numerical behavior.
- Matching requires retained preparation, projected square cells at the
  prepared resolution, `fft = TRUE` and `na.rm = TRUE`. Hard metric radii above
  the prepared buffer are refused because a capped center-distance circle does
  not reproduce the truncated intersection footprint. Weighted SQ remains a
  population weighted spread; hard-radius SQ remains sample SD.
- The guide documents explicit native unmarked prediction and conditional
  logistic raster scores, required nonspatial scenarios, and both opt-ins.
  Native Cox risk scores with `reference = "zero"` are relative selection scores.
  Choice probabilities require explicit alternatives and within-set
  normalization; these raster scores are not occurrence probabilities or a UD.

## Accepted measured gains

The previously pending phase61 pair completed on the same saved phase29
N-mixture fit (15 reordered sites, 100 state and 100 detection calls). Both
contention audits were empty; both modes retained the original fit, with no
refitting or changed starts. Fit/native equality passed at 1e-12 and within-mode
repeated predictions were identical. The pair's production source hashes match;
`default-prediction-path-audit.csv` additionally verifies that the five relevant
prediction/evaluation files still match the live checkout.

| Scope | Live spatial/native | Compact snapshots |
| --- | ---: | ---: |
| 200 repeated predictions | 3.21 s | 0.03 s |
| Snapshot preparation | none | 0.19 s |
| Rprofmem repeated-workload allocations | 7,393,264 bytes | 49,088 bytes |
| Sampled process-tree peak RSS | 1,246,011,392 bytes | 870,002,688 bytes |
| Sampled process-tree peak private commit | 1,548,632,064 bytes | 1,171,742,720 bytes |
| Cold whole-process wall time | 14.687 s | 16.078 s |

Including snapshot preparation, repeated-workload time was 0.22 s versus 3.21 s,
93.1% lower in this single pair. Snapshots retained 17,744 bytes. Process-tree
RSS was sampled every 0.2 s and was 30.2% lower; private commit includes native
allocations that Rprofmem excludes. Peaks occurred 13.984 and 14.250 s into the
processes. Cold process wall time was 9.5% slower for snapshots, so there is no
cold end-to-end speed gain. This is one small saved analysis and one pair, not a
broad or repeated-median memory result. No speed/peak-memory gain is claimed
for the new whole-raster matching option, which can increase work with a large
prepared buffer. Frozen evidence is in phase61/reference-1, candidate-1,
paired-audit.log and spatial-audit.rds.

## Validation and failures retained

- Phase62 final focused validation: 351 assertions, no failures, errors,
  warnings or skips. Covers fixed point/snapshot/native raster contracts,
  unchanged caller inputs, native GLM/Cox equivalence, row keys, normalization
  refusal and clamping. Source manifests are retained in phase62.
- Phase63 focused run 3: 75 assertions, no failures, errors, warnings or skips.
  Existing fixtures cover Gaussian, exponential, exponential-power and fixed
  kernels; independent radii/shapes; means, SQ, SA and higher moments; boundaries,
  central missing masks, constant support and invalid requests. Strict matching
  uses exact, unbinned evaluation at existing raster cell centers, tolerance
  1e-10. Earlier failed runs remain available: missing fixture CRS first, then
  the incorrect center-circle assumption. The original extraction footprint was
  reproduced instead of relaxing a tolerance or changing fitting behavior.
- Phase63 broad focused run: 1,454 assertions, zero failures/errors/skips and
  three existing binomial surface-fixture warnings (`glm.fit` fitted
  probabilities numerically 0 or 1). Those warnings remain visible; no new
  matching/fixed-prediction test warning occurred.
- An owner combined check finished with one export-inventory test failure
  because this batch added `prepare_unmarked_data` without updating the expected
  public export list. Its 4,605 passing assertions and vignette/example success
  are evidence from that frozen candidate, not a clean-check claim. The inventory
  now includes the new export; the focused correction passed 21 assertions,
  zero failures/errors/warnings. Our known-stale first combined attempt was
  stopped early using verified owned process creation times; its logs and reason
  remain in `check-early-stop.json` and attempt-1 files.
- Fresh full attempt 2 and its exact live-source audit completed successfully.
  The combined result, fixture warnings and development-only skips are detailed
  below. Existing ignored compiled objects and old rendered/tangled vignette
  files are listed separately and left intact.

## Recommended use and limits

Use `prepare_prediction()` for repeated outputs at fixed inputs and native
controls. Rebuild the snapshot after changing any rows, covariates, component,
random-effect controls or other native options; it stores fixed outputs rather
than a live refitting engine. Opt into `include_fixed = TRUE` only when retained
training information satisfies the verification contract. For cell-center map
projection intended to follow the fit's neighborhood, also opt into
`match_buffer = TRUE` when its explicit restrictions are suitable.

These changes do not restore missing source coverage, fingerprint raster
contents/CRS, undo distance binning or arbitrary-point versus cell-center
differences. Nonspatial factors, movement predictors, detection layouts, offsets
and repeated-site scenarios remain explicit. Native intervals condition on the
selected scales; scale uncertainty, temporal dependence, individual variation,
and out-of-sample calibration remain separate scientific priorities.

## 14:05 combined package validation completed

Fresh attempt 2 completed successfully in 16 minutes 31.3 seconds:
**R CMD check: 0 errors, 0 warnings, 0 notes**, including examples, donttest
examples, compiled-code inspection, package tests and rebuilt vignettes.
The installed test suite passed 4,763 assertions with zero failures. Its 12
retained test-level warnings are nine missing-CRS warnings in existing
kernel-binning fixtures and three existing binomial surface-fixture warnings.
Its 10 skips are development benchmark-policy/helper tests excluded from the
installed package, not skipped prediction contracts. These test warnings/skips
are distinct from the clean top-level package-check status and are not hidden.

`check-release-2/check-summary.txt`, `R-CMD-check.log` and `installed-tests.log`
retain the complete result. The post-check audit reports all 305 frozen source
files unchanged, no additional package source files, and separately lists old
ignored build products. Thus this combined result covers the current integrated
package sources, including the preserved Bayesian/owner work at the audit time.
No package source changed after that freeze. The subsequent same-grid projection
pair completed serially with empty contention records and no resource stops.

## 14:08 accepted same-grid projection comparison

The exact existing 50x50, resolution-10, seed-99 raster fixture was reused at all
2,500 cell centers, without fitting a model or changing starts. Reference uses
unbinned point preparation; candidate uses `match_buffer = TRUE` raster
transforms. Both evaluate the same Gaussian mean and weighted SQ with independent
scales 60/100 and buffer 150. Source rasters, scales, output columns, cell order,
missing-value masks and production hashes match. Across all cells, maximum
absolute differences were 3.89e-16 for the mean and 8.33e-17 for SQ, passing the
unchanged 1e-10 equivalence gate.

| Scope | Point preparation | Matched raster |
| --- | ---: | ---: |
| Transform elapsed time | 9.95 s | 0.79 s |
| Rprofmem allocations | 524,439,128 bytes | 8,260,024 bytes |
| Sampled process-tree peak RSS | 951,840,768 bytes | 789,594,112 bytes |
| Sampled process-tree peak private commit | 1,250,070,528 bytes | 1,087,864,832 bytes |
| Cold whole-process wall time | 21.188 s | 13.375 s |

Transform time was 12.6 times faster (92.1% lower); allocations were 98.4% lower,
sampled peak RSS 17.0% lower, and cold wall time 36.9% lower in this single pair.
Both fresh processes use one computation thread, identical small metadata setup,
and 0.2-second process-tree memory sampling. The launcher inspected R surfaces
before and during both modes; `competing` records are empty and both exits are
zero. Source/output proof is retained in `projection-audit.rds` and
`projection-audit-2.log`; process data are in `projection-reference-1/` and
`projection-candidate-1/`. No pair is repeated or overwritten.

The reference retains about 25.2 MB of reusable point buffers; the raster's R
object reports 1,304 bytes but contains an external pointer. Those object sizes
are not comparable total-memory measurements. Output matrices also differ in
row-name storage (200,616 versus 40,568 bytes). The sampled process peaks above
are the relevant memory evidence. This one small-grid result does not establish
large-map performance or remove the cost of wide buffers.

The first post-run audit correctly stopped on differing matrix row-name
attributes. The audit now verifies explicit point row IDs against cell order,
then compares mask/value arrays without row labels. Both original failed and
successful audit logs are retained. No package source, output values, numerical
definition or tolerance changed. The first audit shell also inherited an invalid
Windows locale; its warnings remain in that log and the successful audit used C.

## Next priorities

1. Preserve the accepted combined result and both performance pairs. Recheck the
   source manifest after any later owner/package edits before reusing coverage.
2. Retain this working spatial workflow and document remaining source-identity,
   CRS/unit and quantitative coverage diagnostics without silently filtering
   sites or changing existing preparation semantics.
3. Measure large-map runtime and peak memory before promising matched-raster
   efficiency. Use existing analyses and serial uncontended pairs with a useful
   target workload; do not repeat the already accepted small snapshot pair.
4. Evaluate out-of-sample and scale/temporal uncertainty only under a separately
   defined statistical validation design. Fast prediction and native equivalence
   alone do not establish calibration.

Exact check/resume commands and stop gates are in RESUME.md. No package source
edits should overlap the active frozen combined check absent a concrete failure.

## Final operational checkpoint, 14:10 Eastern

This bounded batch is complete ahead of the 15:30 logoff. The one-day automation
`multiscaler-spatial-workflows-until-3-30pm` was successfully deleted, with no
further launch scheduled. No verified task-owned worker remains active.
Unrelated owner `latent_gaussian_20261028/pilot-null.R` processes began after the
projection pair ended and were left untouched. The final audit still matches
all 305 checked source files, with no additional package sources. Diff whitespace
checks for this batch passed. All prior atomic outputs and failed development
audit logs remain available. No staging, commits, version bump, merge or push was
performed by this batch. See the current final handoff at the top of RESUME.md;
earlier operational statements in this report describe their checkpoint times.
