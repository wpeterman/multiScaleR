# Scale Distance

Converts a kernel scale parameter (`sigma`) into a distance in map
units: the distance within which a specified proportion of the kernel
weight falls. Distances can be measured along a line through the point
(one-dimensional) or as the radius of a circle around it
(two-dimensional). Can be used with a fitted `multiScaleR` object to
report optimized scale distances, or with manually supplied kernel
parameters to explore kernel behavior.

## Usage

``` r
kernel_dist(model, prob = 0.9, dimension = c("2d", "1d"), ...)
```

## Arguments

- model:

  A fitted `multiScaleR` object from
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
  When provided, scale distances and 95% confidence intervals are
  returned for all optimized covariates. When omitted, `sigma`,
  `kernel`, and (for `expow`) `beta` must be supplied via `...`.

- prob:

  Numeric between 0 and 1 (exclusive). Cumulative proportion used to
  define the distance. Default: `0.9`. The proportion refers to weight
  along a line for `dimension = "1d"` and weight on a plane for
  `dimension = "2d"`.

- dimension:

  Character. How the distance is measured. `"2d"` (default) returns the
  radius of the circle around the point that holds `prob` of the total
  weight the ideal kernel gives to locations on an unbounded plane.
  Because raster cells lie across a plane, the `"2d"` radius answers
  "within what distance does the ideal kernel assign 90% of its total
  weight?". For distance-weighted kernels it is larger than the `"1d"`
  distance; hard-radius covariates report the same radius for both.
  `"1d"` restores the original line-based distance, enclosing `prob` of
  the area under the kernel curve along a transect through the point.
  See Details.

- ...:

  Additional parameters used when `model` is not supplied:

  `sigma`

  :   Numeric (positive). The kernel scale parameter in the same units
      as the projection of `pts` and `raster_stack` passed to
      [`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).
      For Gaussian kernels this is the standard deviation; for negative
      exponential kernels this is the decay length.

  `kernel`

  :   Character. The kernel function to use. One of `"gaussian"`,
      `"exp"` (negative exponential), `"fixed"` (fixed-radius buffer),
      or `"expow"` (exponential power). Required when `model` is not
      provided.

  `beta`

  :   Numeric (positive). Shape parameter for the exponential power
      kernel. Required when `kernel = "expow"` and `model` is not
      provided. Ignored for all other kernels.

## Value

When `model` is provided: a data frame with one row per optimized
covariate and three columns, all in map units: `Mean` (distance at the
estimated sigma), `2.5%` (distance at the lower 95% limit of sigma), and
`97.5%` (distance at the upper 95% limit of sigma). Values are rounded
to two decimal places. For hard-radius covariates (landscape metrics and
unweighted surface metrics) the row reports the radius itself, whichever
`dimension` is requested.

When kernel parameters are supplied directly via `...`: a single numeric
value giving the distance, in map units, at which the cumulative kernel
weight under the requested dimensional definition reaches `prob`.

## Details

**Default change.** Version 0.7.4 and later use `dimension = "2d"` for
distances, summaries, and kernel plots. Earlier versions used `"1d"`.
Set `dimension = "1d"` explicitly to reproduce earlier reported
distances. This changes reporting, not fitted sigma values, kernel
weights, or likelihoods. Existing simulation extents, projection
cutoffs, and optimization diagnostics retain their explicit 1D
definitions.

**One or two dimensions.** A kernel assigns each raster cell a weight
that declines with the cell's distance from the point. Two distances can
summarize how far that weight reaches:

- `dimension = "1d"` treats the kernel as a curve along a single line (a
  transect) through the point and returns the distance on either side
  that encloses `prob` of the area under that curve. For a Gaussian
  kernel the 90% distance is about 1.65 sigma.

- `dimension = "2d"` accounts for locations lying on a plane, so the
  number of locations at distance *d* grows with the circumference 2 pi
  *d*. It returns the radius of the circle that holds `prob` of the
  ideal kernel weight. For a Gaussian kernel the 90% radius is about
  2.15 sigma.

With sigma = 300 m, a Gaussian kernel gives a 90% distance of about 493
m in one dimension and 644 m in two. On an unbounded plane, locations
within 493 m carry only about 74% of the total weight. When you describe
the area associated with kernel weighting, or compare with values
defined on the plane (for example, a radius holding 90% of the kernel
weight), use `"2d"`. Always state which definition you report.

**Kernel formulas.** In one dimension:

- **Gaussian**: the inverse normal CDF, `qnorm((1 + prob) / 2) * sigma`.

- **Negative exponential**: `-sigma * log(1 - prob)`.

- **Fixed buffer**: `sigma * prob` (the fraction of the buffer radius).

- **Exponential power**:
  `sigma * qgamma(prob, shape = 1 / beta)^(1 / beta)`; both `sigma` and
  `beta` (shape) must be specified.

In two dimensions:

- **Gaussian**: `sigma * sqrt(-2 * log(1 - prob))`.

- **Negative exponential**: `sigma * qgamma(prob, shape = 2)`, about
  3.89 sigma at 90%.

- **Fixed buffer**: `sigma * sqrt(prob)`, the radius of the circle
  holding `prob` of the buffer's area.

- **Exponential power**:
  `sigma * qgamma(prob, shape = 2 / beta)^(1 / beta)`.

Both definitions describe the ideal continuous kernel. They do not
account for truncation at `max_D`, raster boundaries, missing cells, or
cell geometry. The realized share of weight among extracted raster cells
can therefore differ. If the two-dimensional radius approaches or
exceeds `max_D`, consider increasing `max_D` in
[`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
and refitting the model.

Confidence intervals for the fitted `model` case are derived from the
Hessian-based standard errors (or profile-likelihood intervals when
[`profile_sigma`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md)
has been run and stored on the object). Each limit of sigma is converted
to a distance with the same formula. For exponential-power kernels, beta
is held at its fitted estimate during this conversion; the limits do not
include uncertainty in beta or in site-specific retained weight.

[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
of a fitted model reports both distances; see
[`summary.multiScaleR`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md).

## See also

[`summary.multiScaleR`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md),
[`plot.multiScaleR`](https://wpeterman.github.io/multiScaleR/reference/plot.multiScaleR.md),
[`plot_kernel`](https://wpeterman.github.io/multiScaleR/reference/plot_kernel.md)

## Examples

``` r
# \donttest{
## Using package data
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

## Optimize scale
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

## Uses of `kernel_dist`
kernel_dist(model = opt)                     # radius of a circle (2D)
#>       Mean   2.5%  97.5%
#> hab 163.45 122.53 204.37
kernel_dist(model = opt, dimension = "1d")   # original distance along a line
#>       Mean  2.5%  97.5%
#> hab 125.28 93.92 156.65
kernel_dist(model = opt, prob = 0.95)
#>       Mean   2.5%  97.5%
#> hab 186.44 139.76 233.11
kernel_dist(sigma = 500, kernel = 'gaussian', prob = 0.95)
#> [1] 1223.87
kernel_dist(sigma = 100, prob = 0.975, kernel = "exp")
#> [1] 557.16
kernel_dist(sigma = 100, prob = 0.95, kernel = "expow", beta = 1.5)
#> [1] 235.51
kernel_dist(sigma = 100, kernel = "fixed")
#> [1] 94.87
# }

## The two definitions for a Gaussian kernel with sigma = 300
kernel_dist(sigma = 300, kernel = "gaussian")                     # about 644
#> [1] 643.79
kernel_dist(sigma = 300, kernel = "gaussian", dimension = "1d")   # about 493
#> [1] 493.46
```
