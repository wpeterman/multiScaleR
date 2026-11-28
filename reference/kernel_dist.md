# Convert Sigma into an Effective Distance

Converts a kernel scale parameter (sigma) into an *effective distance*
in map units, the distance that holds a stated proportion of the kernel
weight. Use it with a fitted `multiScaleR` object to report the
estimated scales of effect as distances, or with kernel parameters you
supply to explore how far a kernel reaches before fitting. The kernels,
the two distance definitions, and their formulas are described in
[multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md).

## Usage

``` r
kernel_dist(model, prob = 0.9, dimension = c("2d", "1d"), ...)
```

## Arguments

- model:

  A fitted `multiScaleR` object from
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
  When supplied, distances and 95% confidence limits are returned for
  every optimized covariate, and any `sigma`, `kernel`, or `shape`
  passed through `...` is ignored with a warning. When omitted, supply
  `sigma` and `kernel` (and `shape` for `"expow"`) through `...`.

- prob:

  Numeric between 0 and 1 (exclusive). Proportion of the kernel weight
  that the reported distance must hold. Default: `0.9`. Larger values
  give longer distances. Under `dimension = "1d"` the proportion is of
  the weight along a line through the point. Under `"2d"` it is of the
  weight on a plane.

- dimension:

  Character. How the distance is measured. `"2d"` (default) returns the
  radius of the circle around the point that holds `prob` of the ideal
  kernel weight on an unbounded plane. `"1d"` returns the original
  line-based distance, which encloses `prob` of the area under the
  kernel curve along a transect through the point. For distance-weighted
  kernels the `"2d"` radius is longer than the `"1d"` distance.
  Hard-radius covariates in a fitted model report the same radius under
  both. See Details.

- ...:

  Kernel parameters used when `model` is not supplied:

  - `sigma`: positive numeric scalar. The kernel scale parameter, in the
    map units of the points and rasters passed to
    [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).
    For the Gaussian kernel it is the standard deviation, for the
    negative exponential kernel the decay length, and for the fixed
    kernel the buffer radius.

  - `kernel`: character. One of `"gaussian"`, `"exp"` (negative
    exponential), `"expow"` (exponential power), or `"fixed"`
    (fixed-radius buffer). Required.

  - `shape`: positive numeric scalar. Shape of the exponential power
    kernel. Required when `kernel = "expow"` and ignored otherwise.

  - `beta`: deprecated former name of `shape`. It is still accepted,
    with a warning, and will be removed in a future release. Supplying
    both `shape` and `beta` with different values is an error.

## Value

When `model` is supplied, a data frame with one row per optimized
covariate (row names are the covariate names) and three columns, all in
map units and rounded to two decimal places:

- `Mean`: the effective distance at the estimated sigma.

- `2.5%`: the distance at the lower 95% limit of sigma.

- `97.5%`: the distance at the upper 95% limit of sigma.

For a hard-radius covariate (a landscape metric, an unweighted surface
metric, or any covariate fitted with `kernel = "fixed"`), the row
reports the radius itself under either `dimension`. Limits are `NA` when
no interval is available (see Details). A covariate whose sigma estimate
is not finite returns `NaN`.

When kernel parameters are supplied through `...`, a single numeric
value, the distance in map units that holds `prob` of the kernel weight
under the chosen `dimension`, rounded to two decimal places.

## Read the distances

Sigma is a kernel parameter, not a radius. An effective distance puts
the estimate in map units you can compare across kernels and studies.
With sigma = 300 m, a Gaussian kernel holds 90% of its planar weight
within about 644 m (`dimension = "2d"`), while the line-based 90%
distance is about 493 m (`dimension = "1d"`). Use `"2d"` when you
describe the area a kernel draws on, and always state the definition and
the value of `prob` you report. The formulas for each kernel are
tabulated in
[multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md).

Both definitions describe the ideal continuous kernel. They do not
account for truncation at `max_D`, raster edges, missing cells, or cell
geometry, so the share of weight among the cells actually extracted
around a point can differ. If the `"2d"` radius approaches or exceeds
`max_D`, increase `max_D` in
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
and refit the model.

## Reproduce distances from earlier versions

Versions 0.7.4 and later report `dimension = "2d"` distances by default
in `kernel_dist()`,
[`summary()`](https://rspatial.github.io/terra/reference/summary.html),
[`print()`](https://rdrr.io/r/base/print.html), and kernel plots.
Earlier versions reported the line-based distance. The switch changes
reporting only. It does not change fitted sigma values, kernel weights,
or likelihoods. To reproduce distances reported by earlier versions, set
the dimension explicitly:

    kernel_dist(fit, dimension = "1d")
    summary(fit, dimension = "1d")
    print(fit, dimension = "1d")
    plot(fit, dimension = "1d")
    plot_kernel(sigma = 300, kernel = "gaussian", dimension = "1d")

Simulation extents in
[`sim_dat()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat.md),
the `pct_wt` window in
[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md),
and the `max_D` warning of
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
keep their line-based definitions.

## Confidence limits from a fitted model

Each limit of sigma is converted to a distance with the same formula as
the estimate. By default the sigma limits are Wald limits from the
Hessian-based standard errors. When the fitted object stores its own
interval table in `profile_scale_est` (see the [fitted-object help
page](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-class.md)),
those limits are used instead. `summary(fit, profile = TRUE)` converts
profile-likelihood limits to distances in its own output, as described
in
[`summary.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md).
For exponential-power kernels, the shape is held at its estimate during
this conversion, so the limits do not include uncertainty in the shape.
These rules depend on the fitting route:

- With `backend = "profile_glm"` or `"joint_nb2"` and an estimated
  exponential-power shape, distance limits use joint log-scale
  delta-method uncertainty for both sigma and shape. They remain Wald
  limits even when sigma limits are profiled in
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html),
  and they are missing when the fit fails the joint identifiability
  checks. Point distances remain available.

- With `backend = "conditional_logit"`, point distances are available
  but scale uncertainty has not been calibrated, so distance limits are
  missing.

- With `backend = "joint_mixed"`, distance limits convert the log-scale
  Wald limits of sigma and are missing for a scale whose standard error
  is withheld. See
  [multiScaleR-fitting-methods](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-fitting-methods.md)
  for their coverage caution. In `summary(fit, profile = TRUE)`, they
  convert the profile limits instead, which are reported even when the
  standard error is withheld.

- With the default refit backend, point distances remain available when
  a sigma standard error is not finite, but its Wald limits are missing.
  When a smooth parameter sits at a search bound, the whole Wald block
  is withheld, while explicitly marked profile limits remain available.

- A hard-radius covariate reports its radius directly and has no kernel
  shape, including in an exponential-power fit. It has no Hessian-based
  interval, but stored profile limits are kept. See the section on
  smooth kernels and hard radii in
  [multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md).

Native coefficient standard errors remain conditional on the selected
scales. Manually supplied fixed-kernel parameters still return the
requested weight quantile under the chosen `dimension`.

## See also

[multiScaleR-kernels](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-kernels.md)
for the kernel formulas and both distance definitions.
[`summary.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md),
[`plot.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/plot.multiScaleR.md),
and
[`plot_kernel()`](https://wpeterman.github.io/multiScaleR/reference/plot_kernel.md)
report and draw the same distances.
[`vignette("v02-user-guide", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v02-user-guide.md)
converts sigma to effective distances on a worked example.

## Examples

``` r
## Distances from kernel parameters, without a fitted model
# A Gaussian kernel with sigma = 300 m: planar radius holding 90% of the weight
kernel_dist(sigma = 300, kernel = "gaussian")                     # about 644
#> [1] 643.79
# The original line-based definition gives a shorter distance
kernel_dist(sigma = 300, kernel = "gaussian", dimension = "1d")   # about 493
#> [1] 493.46
# Other kernels and proportions
kernel_dist(sigma = 500, kernel = "gaussian", prob = 0.95)
#> [1] 1223.87
kernel_dist(sigma = 100, kernel = "exp", prob = 0.975)
#> [1] 557.16
kernel_dist(sigma = 100, kernel = "expow", shape = 1.5, prob = 0.95)
#> [1] 235.51
kernel_dist(sigma = 100, kernel = "fixed")   # radius holding 90% of the buffer area
#> [1] 94.87

# \donttest{
## Distances from a fitted model
data("pts")
data("count_data")
hab <- terra::rast(system.file("extdata", "hab.tif",
                               package = "multiScaleR"))

kernel_inputs <- kernel_prep(pts = pts,
                             raster_stack = hab,
                             max_D = 250,
                             kernel = "gaussian",
                             verbose = FALSE)
mod <- glm(y ~ hab, family = poisson, data = count_data)
opt <- multiScale_optim(fitted_mod = mod,
                        kernel_inputs = kernel_inputs,
                        verbose = FALSE)

# Planar radius (m) holding 90% of the weight, with limits from the
# 95% interval of sigma
kernel_dist(model = opt)
#>       Mean   2.5%  97.5%
#> hab 163.45 122.53 204.37
# Line-based distances, as reported before version 0.7.4
kernel_dist(model = opt, dimension = "1d")
#>       Mean  2.5%  97.5%
#> hab 125.28 93.92 156.65
# Radius holding 95% of the weight
kernel_dist(model = opt, prob = 0.95)
#>       Mean   2.5%  97.5%
#> hab 186.44 139.76 233.11
# }
```
