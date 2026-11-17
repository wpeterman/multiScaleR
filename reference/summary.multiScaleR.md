# Summarize a multiScaleR fit

Summarizes a fit from
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md):
the estimated scales of effect with standard errors and intervals, the
distances that contain most of each kernel's weight, and the fitted
model at those scales.

## Usage

``` r
# S3 method for class 'multiScaleR'
summary(object, profile = FALSE, dimension = c("2d", "1d"), ...)
```

## Arguments

- object:

  An object of class `multiScaleR` from
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

- profile:

  Logical. If `TRUE`, report profile-likelihood confidence limits for
  sigma where feasible instead of Wald limits. Default: `FALSE`, which
  keeps summaries fast. See "Profile-likelihood intervals" below.

- dimension:

  Character. `"2d"` (default) reports, as the selected distance, the
  radius that holds `prob` of the weight of an ideal kernel on a plane.
  `"1d"` restores the original distance along a line. Both definitions
  remain in the returned object and the printed output.

- ...:

  Optional arguments. `prob`, a number strictly between 0 and 1, sets
  the proportion of cumulative kernel weight used for the distances.
  Default: `0.9`.

## Value

An object of class `summary_multiScaleR`, a list that prints as a
report. Its elements:

- `opt_scale`: data frame with one row per optimized covariate and
  columns `Mean` (estimated sigma, map units), `SE`, `2.5%`, and
  `97.5%`. Its `interval_method` attribute records `"wald"` or
  `"profile"` for each row.

- `opt_shape`: the same table for the shape parameter under
  `kernel = "expow"`, otherwise `NULL`.

- `opt_dist_1d`: distances in map units, with columns `Mean`, `2.5%`,
  and `97.5%`, at which a kernel's cumulative weight along a line
  through the point reaches `prob` (about 1.65 sigma for a Gaussian
  kernel at 90%).

- `opt_dist_2d`: the radius, in map units and with the same columns, of
  the circle around the point that holds `prob` of the weight of an
  ideal, untruncated kernel on a plane (about 2.15 sigma for a Gaussian
  kernel at 90%).

- `opt_dist` and `opt_distance`: two names for the distance table
  selected by `dimension` (2D by default, rather than the 1D definition
  used through version 0.7.3).

- `dimension`: the selected definition, `"2d"` or `"1d"`.

- `prob`: the cumulative weight proportion used for the distances.

- `fitted_mod`: the model refit at the estimated scales (`opt_mod` of
  the fit).

- `kernel`: the kernel name.

- `hard_radius`: names of covariates whose scale is a hard radius.

- `diagnostics`, `warn_message`, `next_run`, `call`: copied from the
  fit; see
  [multiScaleR](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-class.md).

- `backend`: the backend name, present for the analytic and
  conditional-logistic backends.

## Read the summary

Sigma is the kernel's scale parameter, and it is not itself a distance
you can report. The distance tables translate it into the distance that
holds 90 percent (by default) of the kernel's weight, which is usually
the number to report as the scale of effect. Hard-radius covariates
report their radius in both tables. See
[`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md)
and
[multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md)
for when to use each definition and for the formulas.

Kernel-scale limits default to Wald intervals from the Hessian. For
hard-radius landscape and unweighted surface metrics, and every radius
estimated under `kernel = "fixed"`, Hessian standard errors and Wald
limits are shown as `NA`; use `profile = TRUE` for profile-likelihood
limits. An error during profiling is reported rather than returning
invalid limits.

Set `dimension = "1d"` to keep the original `opt_dist` and
`opt_distance` values in existing scripts. For extraction that does not
depend on the setting, use `opt_dist_1d` or `opt_dist_2d`. Neither
choice changes the fitted sigma, coefficients, or likelihood. The 2D
radius is for an ideal continuous kernel; `max_D`, raster edges, missing
cells, and cell geometry can change the realized weight in the prepared
data.

## Profile-likelihood intervals

The default refit backend caches profile results during the current R
session; the analytic backends recompute them. Set
`options(multiScaleR.profile_cache_max_bytes = 64 * 1024^2)` to bound
the retained cache entries (64 MiB by default), or set it to zero to
disable caching and clear retained entries at the next profile request.
Oversized entries and custom refits are computed without the cache. A
failure to store a cache key warns and computes the same intervals
directly, while numerical profile failures still raise errors. The
returned refit scale table records cache availability in its
`profile_cache` attribute. The option bounds retained entries, not
temporary key construction or the memory used by the numerical profile.

## Analytic backends

For `backend = "profile_glm"` or `"joint_nb2"`, default scale limits use
joint information on the log scale and are omitted when the
identifiability checks fail (see
[multiScaleR-fitting-methods](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-fitting-methods.md)).
A returned Wald interval describes the local curvature around the
selected optimum and can extend beyond the searched scale range.
Coverage can fall below the nominal level even for an interior fit, so
inspect the likelihood over the admissible range when scale inference
matters.

With `profile = TRUE`, the analytic summary follows the connected
likelihood region around the fitted scale, and an endpoint at a search
bound is labeled bound-limited. For a supported one-scale NB2 fit,
[`scale_inference()`](https://wpeterman.github.io/multiScaleR/reference/scale_inference.md)
checks the full range and reports the regions it detects. If a
full-range profile suggests a competing optimum, compare additional
starts, including both bounds. The printed coefficient summary still
conditions on the estimated kernels; use
[`vcov.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/vcov.multiScaleR.md)
for the joint covariance.

## See also

[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`print.summary_multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/print.summary_multiScaleR.md),
[`profile_sigma()`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md),
[`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md),
and
[`vignette("v09-estimates-and-intervals", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v09-estimates-and-intervals.md).

## Examples

``` r
data(pts)
data(count_data)
hab <- terra::rast(system.file("extdata", "hab.tif",
                               package = "multiScaleR"))
inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
                      kernel = "gaussian", verbose = FALSE)
fit <- multiScale_optim(glm(y ~ hab, family = poisson, data = count_data),
                        inputs, verbose = FALSE)
s <- summary(fit)
s$opt_scale     # sigma with Wald limits, in map units
#>         Mean       SE     2.5%    97.5%
#> hab 76.16637 9.567451 57.09847 95.23428
s$opt_dist_2d   # radius holding 90 percent of the kernel weight
#>       Mean   2.5%  97.5%
#> hab 163.45 122.53 204.37

# Profile-likelihood limits instead of Wald limits
summary(fit, profile = TRUE)$opt_scale
#>         Mean       SE     2.5%    97.5%
#> hab 76.16637 9.567451 59.86998 99.38553
```
