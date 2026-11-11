# Summarize multiScaleR objects

Summarizes output from `multiScale_optim`.

## Usage

``` r
# S3 method for class 'multiScaleR'
summary(object, profile = FALSE, dimension = c("2d", "1d"), ...)
```

## Arguments

- object:

  An object of class `multiScaleR`.

- profile:

  Logical. If `TRUE`, use profile-likelihood confidence limits for
  \`sigma\` when feasible. Defaults to `FALSE` so summaries remain fast;
  profile results are cached for repeated calls on the same fitted
  object during the current R session.

- dimension:

  Character. `"2d"` (default) reports the radius containing `prob` of
  ideal planar kernel weight as the selected distance. `"1d"` restores
  the original distance along a line. Both definitions remain available
  in the returned object and printed output.

- ...:

  Optional arguments passed to the method (e.g., `prob` for cumulative
  kernel weight threshold).

## Value

An object of class `summary_multiScaleR`. Kernel-scale limits default to
Wald intervals. Scale distances, in map units, come in two forms:
`opt_dist` (also `opt_distance`) holds the selected definition, which
defaults to 2D rather than the 1D definition used through version 0.7.3.
`dimension` records that selection. `opt_dist_1d` always holds the
one-dimensional distance, the distance at which the kernel's cumulative
density along a line through the point reaches `prob` (about 1.65 sigma
for a Gaussian kernel at 90%). `opt_dist_2d` holds the two-dimensional
distance, the radius of the circle around the point that holds `prob` of
the weight of an ideal, untruncated kernel on a plane (about 2.15 sigma
for a Gaussian kernel at 90%). Both tables have columns `Mean`, `2.5%`,
and `97.5%`, and hard-radius covariates report their radius in both. See
[`kernel_dist`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md)
for when to use each. For hard-radius landscape and unweighted surface
metrics, Hessian standard errors and Wald limits are shown as `NA`; use
`profile = TRUE` for profile-likelihood limits. Profiling errors are
reported instead of silently returning invalid limits.

## Details

Set `dimension = "1d"` to preserve the original `opt_dist` and
`opt_distance` values in existing scripts. For a dimension-independent
extraction, use `opt_dist_1d` or `opt_dist_2d`. Neither selection
changes the fitted sigma, coefficients, or likelihood. The 2D radius is
for an ideal continuous kernel; `max_D`, raster edges, missing cells,
and cell geometry can change the realized weight in the prepared data.
See
[`kernel_dist`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md)
for formulas and examples.
