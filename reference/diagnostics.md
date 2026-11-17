# Retrieve diagnostics from multiScaleR objects

Returns structured warning/diagnostic information stored on fitted
`multiScaleR` objects. Diagnostics are populated automatically during
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
and flag potential issues with the optimization result.

## Usage

``` r
diagnostics(object, ...)

# S3 method for class 'multiScaleR'
diagnostics(object, ...)
```

## Arguments

- object:

  An object to inspect. Must be of class `multiScaleR`.

- ...:

  Additional arguments passed to methods.

## Value

A named list with the following elements:

- `max_distance`:

  A list describing whether the estimated scale of effect approaches or
  exceeds `max_D`. Fields include: `triggered` (logical), `variables`
  (names of affected covariates), `effective_distance` (estimated 90%
  one-dimensional kernel distance per covariate; see
  [`kernel_dist`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md)),
  `max_D` (the limit used), `ratio` (`max_D / effective_distance`), and
  `suggested_max_D` (a recommended minimum `max_D` value). Triggered
  when `max_D / effective_distance < 2`.

- `sigma_precision`:

  A list describing whether sigma estimates are imprecise. Fields
  include: `triggered` (logical), `variables` (names of affected
  covariates), `estimate` (sigma means), `se` (standard errors), and
  `se_to_mean` (ratio of SE to mean). Triggered when `SE / Mean >= 0.5`
  for any covariate. `NULL` when no precision concern was flagged.

- `shape_precision`:

  Identical structure to `sigma_precision` but for the shape parameter
  of the exponential power kernel. `NULL` unless `kernel = "expow"` was
  used and a precision concern arose.

- `sample_size`:

  A list with `fitted_n`, the number of observations in the final model,
  and `prepared_n`, the number of observations in the preparation
  object, or model rows for repeated-site `unmarked` data. For
  `unmarked` fits, `prepared_sites` is the number of distinct spatial
  sites in the preparation object. `triggered` is `TRUE` when the fitted
  model uses fewer observations. Candidate scales that change the
  initial model's sample are rejected during optimization.

- `start_stall`:

  For the default refit backend, present only when the optimizer
  returned its starting values, so the reported scales are those starts.
  A search from the default start is flagged whenever it returns that
  start. A user-supplied or screened start is flagged only when the
  search stopped after one evaluation, because such a start may already
  be an optimum. Fields include `reason`, `screened` (whether the start
  came from `start_strategy = "screen"`), `start` (starting sigma values
  in map units), `function_evaluations`, and `message`, which names the
  remedy. `reason` is `"hard_radius"` when every optimized scale is a
  hard radius, as under `kernel = "fixed"`. Its likelihood changes in
  steps, so the gradient is zero almost everywhere. `reason` is
  `"refit_failed"` when no refitted model gave a finite log-likelihood
  at the start. A glmmTMB model with a non-positive-definite Hessian
  reports none. `reason` is `"no_progress"` otherwise. A warning gives
  the same message when the fit is made, and
  [`print()`](https://rdrr.io/r/base/print.html) and
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html)
  repeat it.

- `refit_hessian`:

  For the default refit backend, retains the outer curvature status and
  active smooth-parameter mask. A bound table records estimates,
  lower/upper bounds in reporting units, and bound flags. If any active
  parameter is within `1e-5` log units of a positive bound, joint outer
  Wald SEs are withheld. A zero lower bound is flagged only when
  reached. The original curvature status is retained separately. A
  status of `"hessian_includes_failed_refit"` means a refit failed at
  one of the finite-difference points used for the Hessian, so its
  curvature reflects the failure penalty; Wald SEs are withheld and a
  warning is given. This is a numerical support flag, not a test of
  identification or a calibrated confidence limit. Native coefficient
  covariance and parameter counts are unchanged.

- `radius_profile`:

  For a declared radius grid, one row per combination. `radius`, or one
  `radius:<predictor>` column per hard radius, uses map units. Any
  `sigma:<predictor>` and `shape:<predictor>` columns hold the
  reoptimized smooth bandwidths in map units and unitless shapes,
  respectively. The table also reports native `log_likelihood`, `status`
  (scored or failed), any failure `message`, and a logical `selected`
  flag. With smooth predictors, `smooth_starts` holds a table of
  optimizer and convergence diagnostics for each combination that
  reached the smooth search. An all-failed search carries the profile
  table on the error condition.

- `profile_glm`:

  Present for both opt-in analytic GLM backends. Includes
  `uncertainty_available`, explanatory `reasons`,
  `weak_identifiability`, `boundary`, `condition_number`, log-kernel
  standard errors, competing-mode flags, and attempted/eligible start
  counts. `poisson_limit` marks negligible fitted NB2 extra variance;
  `poisson_compatible` and `poisson_deviance` compare the NB2 fit with a
  Poisson fit at the selected kernel scales. See
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  for thresholds and limitations.

## Details

These diagnostics are evaluated automatically at the end of
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
Console warnings are printed when any diagnostic is triggered and
`verbose = TRUE`. `diagnostics()` provides programmatic access to the
same information without re-running the model. These stored diagnostics
retain the historical 1D 90% distance and threshold, independently of
the new 2D reporting default in plots and summaries. Changing
`dimension` when reporting a fit does not recompute its stored
diagnostics. Use `kernel_dist(object)` to inspect the ideal planar
radius separately.

When `max_distance` is triggered, consider re-running
[`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
with a larger `max_D` value. The suggested minimum is stored in
`suggested_max_D`.

When `sigma_precision` is triggered, the variable may not have a
meaningful scale of effect, or the data may be insufficient to estimate
it precisely. Interpret affected covariates with caution.

## See also

[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)

## Examples

``` r
# \donttest{
data('pts')
data('count_data')
hab <- terra::rast(system.file('extdata', 'hab.tif', package = 'multiScaleR'))

kernel_inputs <- kernel_prep(pts = pts,
                             raster_stack = hab,
                             max_D = 250,
                             kernel = 'gaussian')
#> 
#> No sigma values provided...
#> Creating necessary elements to optimize sigma
#> 
#> Extracting values from raster layers...

mod <- glm(y ~ hab, family = poisson, data = count_data)
opt <- multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs)
#> 
#> 
#> Optimization complete
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent.
#>  Consider increasing max_D in `kernel_prep` to >=250.56 to ensure accurate estimation of scale.
#> 
#> 

## Access diagnostics
diag <- diagnostics(opt)
diag$max_distance$triggered
#> [1] TRUE
diag$sigma_precision
#> NULL
# }
```
