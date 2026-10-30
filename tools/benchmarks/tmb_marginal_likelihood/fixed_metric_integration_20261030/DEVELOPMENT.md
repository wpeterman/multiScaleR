# Spatial-metric integration, version 0.7.6

Continue on **inference-dev**. This batch integrates the current validated package
candidate, including the shared prediction/metric development on which the fixed
adapter depends. Benchmark outputs, failed attempts, and unrelated local work are
preserved separately. The evidence below establishes software behavior, not
posterior approximation accuracy, performance, or interval calibration.

## Shared contracts

- Declared fixed-radius surface and landscape metrics are ordinary numeric
  predictors in supported analytic GLM/NB2, joint mixed, conditional-logistic,
  and Bayesian workflows. Original prepared values and means/sample SDs are
  verified. Fixed radii add no estimated scale/shape or bandwidth prior. Their
  response coefficients still count as ordinary fitted coefficients for AIC.
- `include_fixed = TRUE` opts into checked reconstruction for supported point
  predictions, snapshots, and Bayesian prediction/simulation. The default keeps
  fixed columns explicit. Source, metric, radius, kernel, resolution, original
  buffer extent (`max_D`), physical units, and row/site identities remain checked.
  Direct `predict_scales()` continues to use explicitly supplied covariates.
- Bayesian fixed metrics can enter response and scale formulas. Their original
  preparation normalization remains distinct from the fixed prior-median
  reference for optimized predictors. No per-draw/year/subset restandardization
  or prior change is introduced. Older fixed metadata reuses its retained fitted
  support bound. Legacy mean-only inputs without `scale_vars` retain their
  empty fixed contract and original prediction-data return type.
- Point/raster/snapshot prediction, factor/offset scenarios, native effect plots,
  FFT/focal metric smoothing, and model-selection safeguards share their supported
  definitions. Native SEs condition on selected spatial scales. Bayesian mean
  intervals, posterior predictive noise, zero random effects, and integration
  over random effects remain different quantities. Metadata never grants fit
  eligibility or supplies missing uncertainty.

Start with `vignettes/feature-compatibility.Rmd` and the generated help. The shared
fixed adapter is `R/scale_fixed_adapter.R`; native preparation proof remains in
`R/prepare_unmarked_data.R`; common inference preparation is in
`R/scale_model_prepare.R`. Focused regressions are in
`tests/testthat/test-fixed-inference-compatibility.R` and the companion prediction,
smoothing, selection, and weighted-SQ test files.

## Validation and retained failures

Fresh package check: **Status: 1 NOTE** (zero errors/warnings).
The single NOTE reports only `aic_tab` elapsed time of 5.53 seconds (2.49 CPU
seconds) under concurrent foreign R work. Its original five-second threshold and
complete log are retained. This is an explained development-branch result,
not a clean CRAN timing result; repeat timing checks without competing jobs
before submission. The normal build rebuilt the vignettes and
their index. The exact checked archive's full suite passed **6799
assertions**, with zero failures, 12 existing warnings, and
10 development-only benchmark skips. The broader targeted
repair passed 2,042 assertions with no failures/errors/warnings/skips. Help
generation was idempotent with roxygen2 8.0.0. R was 4.6.1 on Windows.

The focused matrix includes all six surface metrics and 25 landscape metrics,
four response families, Gaussian/exponential kernels, objective/score/Hessian/
prior equivalence to ordinary numeric columns, native starts/cold agreement,
point/snapshot/raster parity, fixed-unit plots, parameter counts, repeated IDs,
changed distance units, and rejected inconsistent inputs. IJI's undefined rows
are explicitly selected (20/36 finite rows); Bayesian rejection is asserted
before that selection, rather than silently dropping observations.

All failed attempts remain in local development records. Gate05's clean package
check was followed by 157 full-suite failures exposing missing legacy-schema
guards. Those guards and native return types were repaired. A bounded buffer
audit found a 1.2152-training-SD discrepancy when a fixed radius of 25 was
recomputed with buffer extent 15 instead of 40; automatic reconstruction now
rejects that change. Pre-override RAM interruptions are recorded as incomplete,
not statistical successes or failures. A separate installed R CMD dispatcher
classification stop is also retained; verified dispatcher hashes now distinguish
those frontends from actual R.dll-bearing interpreters. The final successful gate supersedes
their source coverage without deleting them. Concurrent foreign jobs were left
untouched under explicit user authorization; timings are not comparable.
The first final suite was interrupted when automatic allocation created three
PSOCK children plus its driver. The successful retry limited the runner's cached
hardware budget to four cores, leaving two automatic fit workers; explicit
parallel-parity tests remained intact. Package source and numerical gates were
unchanged. The interrupted run is retained locally.

`validation.json` records source hashes and archive identity. Hashes bind the
Windows validation snapshot; Git may normalize text line endings. The built
archive required exact R/help/test/Rmd identity and recorded only proven native
C++ CRLF-to-LF conversion. Source integrity is separate from scientific validity.

## Next priorities and assertions

Keep private unpenalized SQ limiting/rank-change gates intact. Other optimized
surface moments, discontinuous radius inference, unsupported kernel/method
combinations, public prediction from coefficient-integrated likelihoods, and posterior raster
integration retain their explicit capability limits. Do not substitute native
conditional maps or point scales for unimplemented posterior integration.

Next bounded refinements are native zero-random-effect row requirements
(recognize formula syntax independently of environment) and provider-specific
GAM snapshot controls, with native parity and uncertainty tests. Preserve all
existing starts, cold/Hessian/tail checks, priors, failed fits, and uncertainty
distinctions before advancing those features.

For another code batch, regenerate help, run focused tests, then a fresh complete
package check. On this machine set the RStudio Pandoc path and clear inherited
invalid locale settings. Typical commands are:

```r
devtools::document()
devtools::test()
devtools::check(args = "--no-manual")
```
