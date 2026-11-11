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

A named list with four elements:

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

## Details

These diagnostics are evaluated automatically at the end of
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
Console warnings are printed when any diagnostic is triggered.
`diagnostics()` provides programmatic access to the same information
without re-running the model. These stored diagnostics retain the
historical 1D 90% distance and threshold, independently of the new 2D
reporting default in plots and summaries. Changing `dimension` when
reporting a fit does not recompute its stored diagnostics. Use
`kernel_dist(object)` to inspect the ideal planar radius separately.

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
