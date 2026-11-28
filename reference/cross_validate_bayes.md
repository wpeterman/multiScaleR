# Cross-Validate Bayesian Scale Models

Refit a
[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
model in each analysis split, then score the held-out observations with
the posterior predictive density. Observations that share a predicted
group, such as a held-out site or year, also share its random effect
when they are scored.

## Usage

``` r
cross_validate_bayes(
  formula,
  data,
  kernel_inputs,
  folds,
  scales = NULL,
  family = c("nbinom2", "poisson", "binomial", "gaussian"),
  row_id = NULL,
  priors = list(),
  method = c("slice", "sample", "laplace", "marginal"),
  fit_args = list(),
  assessment_group = NULL,
  new_scale_levels = c("sample", "error", "zero"),
  new_response_levels = c("sample", "error", "zero"),
  new_group_draws = 1L,
  bin_reference = c("training", "fixed"),
  seed = NULL,
  keep_fits = FALSE
)
```

## Arguments

- formula:

  Two-sided response formula in
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  syntax. Use a data column as the response, or
  `cbind(successes, failures)` for grouped binomial data. An
  [`I()`](https://rdrr.io/r/base/AsIs.html) wrapper is also accepted.

- data:

  Data frame of raw observations and formula variables. Its unique row
  names identify observations in `folds`, even when several observations
  share one prepared point.

- kernel_inputs:

  Prepared landscape inputs from
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  or
  [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md).
  For the default training-only rebinning, prepare kernel means with
  `bin = TRUE` and `store_cell_data = TRUE` (both defaults of
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)).
  Weighted SQ and fixed-radius metrics also use the raw cells.

- folds:

  Uniquely named list of splits. Each split holds nonempty, disjoint
  character vectors `analysis` and `assessment` of row names of `data`,
  as returned by
  [`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md)
  or
  [`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md).
  Every analysis set needs at least two observations.

- scales:

  Named list of scale formulas, as in
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md).
  The default `NULL` estimates one shared scale per landscape predictor.

- family:

  Response family, one of `"nbinom2"` (default) or `"poisson"` for
  counts, `"binomial"` for successes and trials, or `"gaussian"` for
  continuous outcomes.

- row_id:

  Optional name of the column of `data` that holds prepared point IDs.
  The default `NULL` uses the row names of `data` as point IDs. Point
  IDs may repeat. They do not replace the unique row names that identify
  observations.

- priors:

  Named list of prior settings passed to
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md).
  Set them before you look at any assessment outcome. The default
  [`list()`](https://rdrr.io/r/base/list.html) uses the documented
  defaults.

- method:

  Computation for every fold. The default `"slice"` needs no optional
  sampling packages. `"sample"` runs NUTS. `"laplace"` and `"marginal"`
  use the approximations, and an importance-corrected marginal fit is
  scored directly with its candidate weights.

- fit_args:

  Named list of further fitting arguments for
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
  chosen from `chains`, `iter`, `warmup`, `control`, `laplace_draws`,
  and `cores`. Settings you leave out take the defaults of
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md).
  Folds always run one after another.

- assessment_group:

  Optional name of a column of `data` that defines joint scoring blocks
  within each assessment split. The default `NULL` scores each whole
  assessment split as one block. Choose the blocks to match your
  prediction question before you see any assessment outcome.

- new_scale_levels, new_response_levels:

  What to do with a group level that is absent from an analysis split,
  set separately for scale and response random effects. The default
  `"sample"` averages over new-group variation by simulation. `"error"`
  refuses an unseen level, and `"zero"` predicts at a zero deviation.
  Groups present in the analysis split keep their fitted effects.

- new_group_draws:

  Positive integer number of new-effect simulations per posterior vector
  when unseen levels are sampled; default 1. Raise it to reduce the
  simulation error of the average over new groups. See "Compare models".

- bin_reference:

  `"training"` (default) computes the representative bin distances from
  the unique analysis points of each fold and freezes them for its
  assessment points. It needs the raw cells, so it refuses inputs
  prepared with `store_cell_data = FALSE`, such as lean inputs from
  [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md).
  `"fixed"` uses the prepared bin distances as given, so lean inputs of
  kernel means work with it. Use `"fixed"` only when the supplied bins
  were fixed independently of all analysis and assessment data.

- seed:

  Optional positive integer seed. It generates distinct fitting and
  scoring seeds for every fold, and your own random-number state is
  restored afterward.

- keep_fits:

  Logical. The default `FALSE` stores compact fold records only. `TRUE`
  also keeps every fold's complete fit in `fits`, for investigation, at
  a cost in memory.

## Value

An object of class `multiScaleR_bayes_cv`, which inherits from
`multiScaleR_cv` (see
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)).
It is a list with these elements.

- `summary`: a one-row data frame with the pooled score. `estimate` is
  the sum of joint negative log scores divided by the number of
  assessment observations (`unit` is
  `"joint_block_density_per_observation"`), so lower is better.
  `n_assessment`, `n_scored`, `weight_sum`, and `complete` describe its
  basis. `estimate` is `NA` unless every fold scored.

- `metrics`: the same score for each `fold`.

- `predictions`: one row per assessed observation, with `fold`,
  `row_id`, `point_id`, `observed`, `trials`, `assessment_group`, the
  posterior mean response `prediction` (a per-trial probability for
  binomial models), `marginal_neg_log_score` for that observation alone,
  `new_scale_level` and `new_response_level` flags, and `eligible`,
  which is `TRUE` only in a scored fold.

- `joint_scores`: one row per fold and assessment block, with its size
  `n` and joint `neg_log_score`.

- `fold_results`: one row per fold, with `n_analysis`, `n_assessment`,
  and `status`, which is `"scored"`, `"diagnostic_failed"`, or
  `"failed"` (an error at any stage, recorded in `errors`).

- `diagnostics`: for each fold, the fit's stored checks (`raw`), the
  [`diagnose_bayes_scales()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_scales.md)
  report for sampled fits (`derived`), the overall `passed` flag, and
  the importance records of corrected fits.

- `errors`, `warnings`: for each fold, the failing `stage` and `message`
  (or `NULL`), and the warnings captured during the fold.

- `training`: for each fold, the analysis rows and point IDs, the fit
  and score seeds, the frozen normalization (`center`, `spread`, and
  those of fixed-radius metrics), the bin distances and their reference,
  and the priors.

- `prediction_contracts`: for each fold, the unseen groups, the number
  of posterior vectors, `new_group_draws`, and whether importance
  weights were used.

- `targets`, `folds`, `settings`, `assessment_counts`, `complete`,
  `call`, `integrity_sha256`: the scored observations, the validated
  splits, the settings, how often each row was assessed, whether every
  fold scored, the call, and a checksum that
  [`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
  uses to detect altered records.

- `fits`: the complete fold fits when `keep_fits = TRUE`, and `NULL`
  otherwise.

## Details

This function is experimental. Its arguments and results may change in a
later release without a deprecation period.

## What the score measures

The joint score of a block is the negative log of the posterior average
of the product of the observation densities in that block. Lower is
better. New effects are drawn once per group in each joint parameter
vector and shared by every assessment observation in that group. For
grouped binomial data, each observation contributes one normalized
probability of its count of successes. The sum of the marginal row
scores is a different target. The pooled score weights each block by its
number of observations. It is the only cross-validation score in the
package that averages over uncertainty in the coefficients, the scales,
and the group effects.

## Choose folds for your prediction question

Each fit learns the normalization, formula bases, coefficients, and
scales from its analysis rows only. Repeated visits can share one point
buffer. To predict new sites or years, hold out whole sites or years,
for example with
[`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md)
on the site labels. Where groups cross, make each assessment block cover
the observations whose joint response you want to predict.
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md)
discusses fold design.

## Handle failed folds

A fold that fails its native or derived diagnostics is marked
`"diagnostic_failed"` and leaves the result incomplete, though its
computed scores stay available for investigation. Slice sampling without
the posterior package can score folds with diagnostics recorded as
unavailable, so install posterior to check mixing and simulation
precision. Laplace and marginal scores inherit the limits of their
approximation and its checks. Importance-corrected fits are scored with
their normalized candidate weights. Passing numerical checks does not
establish predictive calibration. Every fold runs the full fitting
method, so a run takes about as long as one fit per fold. When any fold
is not scored, the function returns the result and gives one warning
that counts the failed folds and quotes their first distinct reasons.
Read `errors` for the stage and message of every failed fold.

## Compare models

Use
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
with `metric = "neg_log_score"` to rank complete Bayesian results
computed on the same targets, folds, blocks, and new-level policies.
Poisson and NB2 results can be compared, because both score probability
mass on the same counts. Training AIC, AICc, and BIC selection does not
apply to posterior fits, so use an outer split to evaluate a
model-selection rule. The average over unseen group effects uses a
finite number of simulated effects. When two candidates score very
closely, rerun both with a larger `new_group_draws` and the same `seed`
and fitting settings, and treat the comparison as unresolved if the
ranking changes.

## See also

[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
[`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md),
[`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md),
[`diagnose_bayes_scales()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_scales.md),
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md),
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md);
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md).

## Examples

``` r
# \donttest{
data("landscape_counts")
data("surv_pts")
land <- terra::rast(system.file("extdata", "landscape.tif",
                                package = "multiScaleR"))
# Keep the raw cells (the default) so each fold rebins from its own
# analysis points
inputs <- kernel_prep(pts = terra::vect(surv_pts), raster_stack = land[["land1"]],
                      max_D = 1700, kernel = "gaussian", verbose = FALSE)

# Two folds of alternating rows
folds <- grouped_cv_folds(rep(1:2, length.out = nrow(landscape_counts)),
                          k = 2, ids = rownames(landscape_counts))

# Very short slice chains keep this example quick. With the posterior
# package installed, they usually fail the diagnostics, so the folds are
# marked "diagnostic_failed", the pooled score is NA, and a warning says
# so. Use the default chains and iter in an analysis.
cv <- cross_validate_bayes(counts ~ land1, landscape_counts, inputs, folds,
                           family = "poisson",
                           priors = list(scale_median = 300),
                           method = "slice",
                           fit_args = list(chains = 2, iter = 200), seed = 1)
#> Warning: cross_validate_bayes() returned an incomplete result: 2 of 2 folds failed, so the pooled scores are NA. Reason: "The fold did not pass its posterior diagnostics." (fold1, fold2). Inspect `$errors` of the result for the stage and message of every failed fold.
cv$fold_results    # status of each fold
#>    fold n_analysis n_assessment            status
#> 1 fold1         50           50 diagnostic_failed
#> 2 fold2         50           50 diagnostic_failed
cv$joint_scores    # joint negative log score of each assessment block
#>    fold assessment_group  n neg_log_score
#> 1 fold1       whole_fold 50      105.5412
#> 2 fold2       whole_fold 50       90.7778
cv$summary         # pooled score per observation, lower is better
#>          metric estimate                                unit n_assessment
#> 1 neg_log_score       NA joint_block_density_per_observation          100
#>   n_scored weight_sum complete
#> 1        0          0    FALSE
# }
```
