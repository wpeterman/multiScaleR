# Profile Model Fit Across Sigma Parameter Space

Evaluates the model log-likelihood and AICc at a series of sigma values
spanning the optimization range for each spatial covariate. This
provides a diagnostic view of the likelihood surface and helps assess
whether the optimized sigma is at a clear minimum or on a flat plateau.

## Usage

``` r
profile_sigma(
  x,
  n_pts = 10,
  metric = c("AICc", "LL"),
  verbose = TRUE,
  n_cores = NULL,
  spacing = c("log", "linear"),
  sigma_values = NULL,
  sigma_range = NULL
)
```

## Arguments

- x:

  A fitted `multiScaleR` object returned by
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

- n_pts:

  Positive integer (\>= 3). Number of sigma values to evaluate for each
  covariate along the profile grid. More points give a smoother profile
  at the cost of additional model refits. Default: `10`. Ignored when
  `sigma_values` is supplied.

- metric:

  Character. Metric to display on the y-axis of the profile. One of
  `"AICc"` (default) or `"LL"` (log-likelihood). The optimized sigma
  minimizes negative log-likelihood, so the `"LL"` profile should peak
  at the optimized value and `"AICc"` should show a minimum.

- verbose:

  Logical. Print per-covariate progress messages. Default: `TRUE`.

- n_cores:

  Optional positive integer. Number of CPU cores to use when profiling
  covariates in parallel. Parallel profiling is applied across
  covariates with a PSOCK cluster. Default: `NULL` (serial profiling).

- spacing:

  Character. Spacing of the automatically generated sigma grid. `"log"`
  (default) concentrates evaluation points at small sigma values where
  the likelihood surface typically changes more rapidly. `"linear"`
  provides equal spacing across the range.

- sigma_values:

  Optional positive numeric vector of sigma values (in the original
  projection units) to evaluate directly. When supplied, `spacing`,
  `sigma_range`, and `n_pts` are ignored. Duplicate values are removed
  and remaining values are sorted before profiling. Must contain at
  least 3 unique values.

- sigma_range:

  Optional positive numeric vector of length 2 giving the minimum and
  maximum sigma values (in projection units) for the generated profile
  grid. Defaults to the optimization range stored in `x` (`min_D` to
  `max_D`). Ignored when `sigma_values` is supplied.

## Value

A list of class `sigma_profile` containing:

- profiles:

  A data frame with columns `variable`, `sigma`, `LL`, and `AICc`.

- opt_sigma:

  A named numeric vector of the optimized sigma for each covariate.

- metric:

  The metric used for profiling.

- spacing:

  The profile grid type: `"log"`, `"linear"`, or `"custom"`.

- sigma_grid:

  The sigma values evaluated for each covariate.

## Details

For each spatial covariate, sigma is varied across a sequence of
candidate values, while all other sigma values are held at their
optimized values. By default this is a log-spaced sequence from the
minimum to maximum distance considered during optimization. Users can
instead request linear spacing with `spacing = "linear"` or supply exact
values with `sigma_values`. At each evaluation point the model is refit
and the log-likelihood extracted. AICc uses all free fitted likelihood
parameters (including ancillary and unmarked detection or variance
parameters), plus estimated sigma and shape parameters, and the existing
model observation count. Parameter counting follows
[`aic_tab`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md).
AICc is `NA` when `n <= K + 1`. If all refits return the optimizer
penalty, profiling stops with an error instead of returning an
artificial flat profile. Keep the model formula and data available to
the stored model call, or provide a working `refit_fn` to
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

Log-spacing concentrates evaluation points at smaller sigma values where
the likelihood surface often changes more rapidly, and spaces them out
at larger sigma values where the surface tends to be flatter.

Linear spacing can be useful when the profile needs equal resolution
across a specific sigma interval. User-supplied `sigma_values` are
sorted and duplicate values are removed before profiling to avoid
redundant refits.

## See also

[`plot.sigma_profile`](https://wpeterman.github.io/multiScaleR/reference/plot.sigma_profile.md),
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)

## Examples

``` r
# \donttest{
data('pts')
data('count_data')
hab <- terra::rast(system.file('extdata',
                   'hab.tif', package = 'multiScaleR'))

kernel_inputs <- kernel_prep(pts = pts,
                             raster_stack = hab,
                             max_D = 250,
                             kernel = 'gaussian')
#> 
#> No sigma values provided...
#> Creating necessary elements to optimize sigma
#> 
#> Extracting values from raster layers...

mod <- glm(y ~ hab,
           family = poisson,
           data = count_data)

opt <- multiScale_optim(fitted_mod = mod,
                        kernel_inputs = kernel_inputs)
#> 
#> 
#> Optimization complete
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent.
#>  Consider increasing max_D in `kernel_prep` to >=250.56 to ensure accurate estimation of scale.
#> 
#> 

## Profile sigma
prof <- profile_sigma(opt)
#> Profiling 'hab' (1 of 1)...
#> Done.
plot(prof)


## More evaluation points
prof <- profile_sigma(opt, n_pts = 20)
#> Profiling 'hab' (1 of 1)...
#> Done.
plot(prof)


## Linearly spaced values between explicit limits
prof <- profile_sigma(opt, n_pts = 8, spacing = "linear",
                      sigma_range = c(25, 250))
#> Profiling 'hab' (1 of 1)...
#> Done.
plot(prof)


## User-supplied sigma values
prof <- profile_sigma(opt, sigma_values = c(25, 50, 100, 150, 250))
#> Profiling 'hab' (1 of 1)...
#> Done.
plot(prof)


## Profile log-likelihood instead of AICc
prof <- profile_sigma(opt, metric = "LL")
#> Profiling 'hab' (1 of 1)...
#> Done.
plot(prof)

# }
```
