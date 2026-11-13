# Bayesian CV with refitted common and annual bandwidths

This study compares `cross_validate_bayes()` with an independent prediction
calculation after actual fold refits. All simulation and fitting choices were
fixed before the production runs. Short 40-draw-per-chain smoke fits verified
the data and parameter mappings and were excluded. The smoke fits also
motivated the composite quadrature used in production.

## Design

The 16 independent datasets cross Gaussian and NB2 responses with a common
or annual bandwidth generator. Four datasets were generated in each cell.
The common bandwidth is 260 m. In the annual generator, year effects have
SD 0.55 on log bandwidth. A site random intercept affects the response in
both generators. Each site-year pair has two visits with unique observation
IDs and a shared prepared buffer. The raw buffers contain twelve distances
with one cell per distance bin, and site-year specific cover values. This
geometry permits an exact independent check of bin representatives. It is
not a test of geometry for arbitrarily irregular raster extraction.

Each dataset has two predeclared splits. Each split has one joint block of
eight visits.

| Fold | Analysis set | Assessment set | Effect integrated at prediction |
|---|---|---|---|
| Whole year | Years 1 through 4 | All visits in year 5, at four known sites | Annual candidate samples a new bandwidth deviation. The common candidate uses the fitted common bandwidth. Both retain fitted site effects. |
| Seen years, new site | All visits except those at site 7 | All visits to site 7, spanning years 1 through 4 | Both candidates sample one new site intercept shared by its eight visits. The annual candidate uses its fitted bandwidth deviations for those years. |

Both candidates use `y ~ habitat + offset(offset) + (1 | site)` and the same
priors. The annual candidate uses `list(habitat = ~1 + (1 | year))`.
The Gaussian and NB2 dispersion priors are centered at their respective
generating values. Four NUTS chains run 3,000 iterations each with 1,500
warmup iterations and one core, leaving 6,000 posterior vectors. Package
fit and derived scale diagnostics are applied unchanged. Unsuccessful fits,
failed diagnostics, warnings, missing references, and pending attempts
remain in the fold inventory. The two models and two folds reuse each
dataset. Only the independent datasets supply the variability estimates.

## Reference and scoring

- **Conditional on the actual refit.** An independent implementation uses
  the retained draws to rebuild the Gaussian kernel, training transformation,
  offsets, year bandwidths, and known-site intercepts. It checks the
  observed responses, prepared point IDs, bins and normalization against
  the analysis rows. It reconstructs the package's default augmentation
  at the identical score seed. The reference integrates a new Gaussian
  site intercept analytically for Gaussian observations and uses composite
  Gauss-Legendre quadrature over a standard normal on [-12, 12] otherwise.
  The two quadrature passes use eight and sixteen nodes per subinterval.
  Score differences against this reference measure new-group integration
  error at fixed draws. They do not estimate NUTS uncertainty.
- **Oracle with known generating parameters.** The oracle knows the fixed
  generating parameters and the realized bandwidths of years observed in
  training. For each assessed known site, it learns the intercept from
  training responses only. It integrates new years and sites under the
  generating distributions. Gaussian site conditionals use conjugacy.
  NB2 site conditionals use mode-centered quadrature with proposal-density
  adjustment. The oracle evaluates 65- and 129-node Gauss-Hermite
  quadratures. The first production oracle used unshifted Gaussian kernel
  exponents, which underflowed for very small bandwidths at some quadrature
  nodes. `repair-oracle.R` recomputes only the oracle after subtracting the
  nearest squared distance, which preserves the kernel ratio and keeps its
  denominator positive. The original oracle outputs remain in the per-fit
  records and in the `*-frozen-oracle.csv` exports. The repair has its own
  source fingerprint and compares every formerly finite score and marginal
  probability against its repaired value. The oracle has access to latent
  and global information a fitted model must learn. Realized score
  differences are descriptive and may favor a fitted model on a particular
  block.
- **Rare joint outcome stress.** The same fitted vectors are rescored for
  the actual assessment responses and for a prespecified upward shift in
  all eight responses. The Gaussian shift is four true residual SDs. The
  NB2 outcome is twice the original count plus eight. Four independent
  seeds evaluate nested budgets of 1, 8 and 32 new-effect vectors per
  posterior vector. This vectorized independent experiment uses the same
  conditional estimator and a different random-number ordering from the
  package when the augmentation budget exceeds one. The shifted data are
  never used for fitting, score
  comparisons between candidate models, or calibration summaries. The
  density contribution count is a concentration diagnostic on the
  integrated density across finite posterior vectors. It is not a
  sampling ESS and does not include within-vector integration error.
- **Predictive intervals.** The reference evaluates marginal predictive
  cumulative probabilities for each visit, randomized uniformly within
  the discrete mass for NB2 outcomes. The fraction between .05 and .95
  describes central 90% predictive coverage on the probability scale.
  These visits share group effects. Summaries first aggregate within each
  dataset and use the dataset as the independent unit. Four datasets per
  cell support descriptive comparisons only.

## Reproduce

Run from the package root against integrated source
`8f7a26a8ead71c29c89e725e40e61829fa11f0e7`. The archived configuration
records the exact scripts, source content, binary SHA-256, dependencies
and settings. Compile the package native library before running the study.
Set `OMP_NUM_THREADS=1` and `OPENBLAS_NUM_THREADS=1`.
The experiment used R 4.6.1 and verified RTMB 2.0 and tmbstan 1.2.1.
The optional `statmod` package supplies quadrature nodes for the
benchmark only.

```sh
Rscript tools/benchmarks/bayesian_scales/cv_refits_20261113/run.R gaussian
Rscript tools/benchmarks/bayesian_scales/cv_refits_20261113/run.R nbinom2
Rscript tools/benchmarks/bayesian_scales/cv_refits_20261113/run.R assess
Rscript tools/benchmarks/bayesian_scales/cv_refits_20261113/repair-oracle.R
Rscript tools/benchmarks/bayesian_scales/cv_refits_20261113/summarize.R
Rscript tools/benchmarks/bayesian_scales/cv_refits_20261113/verify.R frozen
```

The two family runs can run concurrently with one single-thread process
each. They write one atomic compact record per dataset and candidate.
The assessment, repair, and summary commands make flat inventories and
summary tables. The verification command checks the exports against the
original records. Use `verify.R records` to check the archive after later
package-source edits. The assessment command requires the original source
identity. Reproduce in a clean copy of that revision with the archived
scripts, native DLL and dependencies. Keep all outcomes when calculating
availability and use only matched scored pairs for model contrasts.
