# Kernels and scales of effect

`multiScaleR` summarizes the landscape around each survey point with a
*kernel*, a function that gives nearby raster cells more weight than
distant ones. The width of the kernel is its scale parameter, sigma,
measured in map units. This page defines the four kernel families,
explains which covariates use smooth kernel weights and which use a hard
radius, shows how sigma is converted into a distance you can report, and
explains how to choose `max_D`. For worked examples, see the section
"Understand the kernel and its parameters" in
[`vignette("v02-user-guide", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v02-user-guide.md).

## What a kernel does

The *scale of effect* is the distance over which a landscape variable
influences a response measured at a survey point. Instead of fitting one
model per buffer radius, `multiScaleR` weights every raster cell around
a point by its distance and estimates the kernel width from your data,
along with the regression coefficients.

For a survey point and a raster layer,
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
extracts the cells within `max_D` map units of the point and records the
distance `d` from the point to each cell center. A kernel-weighted
covariate is then the weighted mean of those cell values `v`:

\$\$x = \frac{\sum_i w(d_i) v_i}{\sum_i w(d_i)}\$\$

The sums run over cells with a finite value. A missing cell drops out of
both sums, so the remaining weights renormalize. Because the weight
function appears above and below the line, any normalizing constant
cancels and only the shape of `w(d)` matters. The model then uses this
covariate after centering and scaling, as described in
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).

Sigma is estimated by
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
(or
[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)).
It is always in the linear units of the projection, usually meters.

## Choose a kernel family

Set the kernel with the `kernel` argument of
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).
All covariates prepared together share one kernel, so comparing kernels
means one
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
call per kernel. The weight functions below come from the package code.

|  |  |  |  |
|----|----|----|----|
| Kernel | `kernel =` | Weight of a cell at distance `d` | Estimated parameters |
| Gaussian (default) | `"gaussian"` | `exp(-d^2 / (2 * sigma^2))` | sigma |
| Negative exponential | `"exp"` | `exp(-d / sigma)` | sigma |
| Exponential power | `"expow"` | `exp(-(d / sigma)^shape)` | sigma and shape |
| Fixed buffer | `"fixed"` | 1 when `d < sigma`, 0 otherwise | sigma, which is the buffer radius |

How to read the parameters:

- **Gaussian.** Sigma is the standard deviation of a normal curve. A
  cell at distance sigma gets about 61% of the weight of a cell at the
  point, and a cell at twice sigma gets about 14%.

- **Negative exponential.** Sigma is the decay length. A cell at
  distance sigma gets about 37% (1/e) of the weight at the point. The
  weight drops faster near the point than under a Gaussian kernel and
  has a longer tail.

- **Exponential power.** A cell at distance sigma always gets about 37%
  of the weight at the point, and `shape` sets the form of the curve.
  With `shape = 1` the kernel is the negative exponential with the same
  sigma. With `shape = 2` it has the Gaussian form, with a standard
  deviation of sigma / sqrt(2). Larger values flatten the top and
  steepen the edge, so the kernel approaches a fixed buffer of radius
  sigma. Values below 1 give a sharper peak and a heavier tail.
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  starts the shape at 2, and the default refit search of
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  keeps it between 0.75 and 50.

- **Fixed buffer.** Every cell whose center lies within sigma gets equal
  weight, so the covariate is the plain mean of those cells.

Start with the Gaussian kernel. Fit the other kernels when you have a
reason to test a different decay shape, and compare the fits by AICc
with
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md).
Not every fitting route accepts every kernel, as
[multiScaleR-fitting-methods](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-fitting-methods.md)
explains.

## Smooth kernels and hard radii

A covariate's scale is either a smooth kernel width or a hard radius,
and the two behave differently during estimation.

**Smooth kernel weights.** Kernel means from
[`kernel_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
and kernel-weighted surface metrics from
`surface_var(..., weighted = TRUE)` use the Gaussian, negative
exponential, or exponential power weights above. Each weight changes
continuously with sigma, so the likelihood is a smooth function of
sigma.
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
can then report a Hessian-based (Wald) standard error and confidence
interval for sigma.

**Hard radius.** Some covariates count each cell fully or not at all:

- Every landscape metric from
  [`landscape_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md).

- Every surface metric from
  [`surface_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  with `weighted = FALSE` (the default).

- Every optimized covariate prepared with `kernel = "fixed"`, including
  kernel means and weighted surface metrics.

The estimated scale of such a covariate is a radius in map units. The
covariate changes only when the radius passes the distance of a cell
center, so the likelihood is a step function of the radius. It is flat
between steps and jumps at them, which leaves no useful curvature for a
Hessian standard error. The default summary therefore shows `NA` for the
standard error and Wald interval of a hard radius. Examine the
likelihood profile instead, with
[`profile_sigma()`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md)
or `summary(fit, profile = TRUE)`, and report the estimate as a radius
in map units. An interval read from a profile cutoff does not by itself
establish coverage for a stepped likelihood. Because a flat objective
can stall a search at its start,
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
screens several starting values by default when the kernel is `"fixed"`.

Under `kernel = "expow"`, a hard-radius covariate has no shape
parameter. Its reported shape and shape standard error are `NA`, and
model selection leaves that unused position out of the parameter count.

A metric given a fixed `radius` in
[`landscape_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
or
[`surface_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
has no scale to estimate. It is computed once and enters the model as an
ordinary covariate.

## Report an effective distance

Sigma is a parameter, not a radius. The same sigma reaches a different
distance under each kernel, so report an *effective distance*, the
distance that holds a stated share `prob` of the kernel weight (90% by
default).
[`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md)
makes the conversion,
[`plot_kernel()`](https://wpeterman.github.io/multiScaleR/reference/plot_kernel.md)
draws it, and
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
of a fitted model reports it next to each sigma.

Two definitions are available through the `dimension` argument:

- **Planar, `dimension = "2d"` (the default since version 0.7.4).**
  Raster cells lie on a plane, and the number of cells at distance `d`
  grows with the circumference 2 pi `d`. The planar distance is the
  radius of the circle around the point that holds `prob` of the total
  kernel weight. Use it when you describe the area a kernel draws on.

- **Line-based, `dimension = "1d"`.** The kernel is treated as a curve
  along a single transect through the point. The distance encloses
  `prob` of the area under that curve. This was the default before
  version 0.7.4.

|  |  |  |  |
|----|----|----|----|
| Kernel | Line-based distance (`"1d"`) | Planar radius (`"2d"`) | At `prob = 0.9`, 1D and 2D |
| Gaussian | `qnorm((1 + prob) / 2) * sigma` | `sigma * sqrt(-2 * log(1 - prob))` | 1.645 and 2.146 times sigma |
| Negative exponential | `-sigma * log(1 - prob)` | `sigma * qgamma(prob, shape = 2)` | 2.303 and 3.890 times sigma |
| Exponential power | `sigma * qgamma(prob, shape = 1 / shape)^(1 / shape)` | `sigma * qgamma(prob, shape = 2 / shape)^(1 / shape)` | Depends on `shape` |
| Fixed buffer | `sigma * prob` | `sigma * sqrt(prob)` | 0.900 and 0.949 times sigma |

For example, a Gaussian kernel with sigma = 300 m holds 90% of its
planar weight within about 644 m, while the line-based 90% distance is
about 493 m. On an unbounded plane, the circle of radius 493 m holds
only about 74% of the weight. A hard-radius covariate reports its radius
under either definition. Always state the definition and the share of
weight you report.

Both definitions describe the ideal, continuous kernel on an unbounded
plane. They ignore truncation at `max_D`, raster edges, missing cells,
and cell geometry, so the share of weight held by the cells actually
extracted around a point can differ. The `pct_wt` smoothing window of
[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md),
the default `max_D` of
[`sim_dat()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat.md),
and the `max_D` warning of
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
keep the line-based definition.

## Choose max_D

`max_D`, set in
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md),
does three jobs.

- **It sets the neighborhood.** Each point gets every raster cell that
  intersects a circle of radius `max_D` around it, including partly
  covered cells at the edge, which keep their full kernel weight. Cells
  beyond `max_D` are never seen, so each kernel is cut off at `max_D`
  and its weights renormalize over the extracted cells.

- **It bounds the search.**
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  estimates sigma, or a hard radius, between `min_D` (the spacing of
  adjacent cell centers, about one cell width) and `max_D`.
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  places priors on sigma without these hard bounds, but its covariates
  still use only the cells extracted within `max_D`.

- **It sets the cost.** A larger `max_D` extracts more cells per point,
  which costs preparation time and memory. With `bin = TRUE` (the
  default) the cost of optimizing kernel means no longer depends on
  `max_D`. The raster must cover each point's full buffer.

Set `max_D` above the largest scale you consider plausible. After
fitting,
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
flags an estimate whose line-based 90% distance is more than half of
`max_D`. It prints a warning (when `verbose = TRUE`), records the
`max_distance` entry of
[`diagnostics()`](https://wpeterman.github.io/multiScaleR/reference/diagnostics.md),
and suggests a larger `max_D` in its `next_run` element. For a Gaussian
kernel the flag appears once sigma exceeds about `max_D / 3.3`, and for
a negative exponential kernel once it exceeds about `max_D / 4.6`. So
for a plausible upper sigma `s`, choose `max_D` of at least 3.3 `s`
(Gaussian) or 4.6 `s` (negative exponential). Also compare the planar
radius from `kernel_dist(sigma = s, kernel = ...)` with `max_D`. When
that radius approaches `max_D`, the extracted cells hold noticeably less
of the kernel weight than the ideal kernel describes. If an estimate
comes back near `max_D`, increase `max_D`, run
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
again, and refit. Do not shrink a scientifically chosen `max_D` only to
save memory.
[`estimate_multiscale_ram()`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md)
describes leaner preparation options.

## Functions that use kernels

|  |  |
|----|----|
| Function | What it does with the kernel |
| [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md), [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md) | Set the kernel and `max_D`, and extract cells and distances |
| [`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md) | Declare kernel means, kernel-weighted surface metrics, and hard-radius metrics |
| [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md) | Estimate sigma (and shape) with the model |
| [`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md) | Convert sigma into an effective distance |
| [`plot_kernel()`](https://wpeterman.github.io/multiScaleR/reference/plot_kernel.md) | Draw a kernel without a fitted model |
| [`profile_sigma()`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md) | Trace the likelihood across sigma or a hard radius |
| [`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md) | Apply fitted kernels to whole rasters for mapping |
| [`sim_dat()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat.md) | Simulate responses at known scales of effect |

## See also

- [`vignette("v02-user-guide", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v02-user-guide.md)
  for kernel plots, effective distances, and kernel comparison on a
  worked example.

- [`vignette("v04-landscape-metrics", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v04-landscape-metrics.md)
  and
  [`vignette("v05-surface-metrics", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v05-surface-metrics.md)
  for hard-radius and kernel-weighted metrics.

- [`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md),
  [`plot_kernel()`](https://wpeterman.github.io/multiScaleR/reference/plot_kernel.md),
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md),
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
  [`profile_sigma()`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md),
  [multiScaleR-fitting-methods](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-fitting-methods.md),
  [multiScaleR-glossary](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-glossary.md).
