# Simulate data for optimizing scales of effect

Generates simulated response data with known kernel-weighted landscape
covariates at a controlled scale of effect. Useful for testing and
demonstrating
[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
with data where the true parameters are known.

## Usage

``` r
sim_dat(
  alpha = 1,
  beta = NULL,
  kernel = c("gaussian", "exp", "expow", "fixed"),
  type = c("count", "count_nb", "occ", "gaussian"),
  StDev = 0.5,
  n_points = 50,
  min_D = NULL,
  raster_stack = NULL,
  sigma = NULL,
  shape = NULL,
  max_D = NULL,
  user_seed = NULL,
  ...
)
```

## Arguments

- alpha:

  Numeric scalar. Intercept term for the linear predictor. Default: `1`.

- beta:

  Numeric vector of slope coefficients, one per raster layer in
  `raster_stack`. Length must equal `nlyr(raster_stack)`.

- kernel:

  Character. Kernel function used to weight raster values by distance
  when generating the true covariate values. One of `"gaussian"`
  (default), `"exp"` (negative exponential), `"expow"` (exponential
  power), or `"fixed"` (fixed-radius buffer).

- type:

  Character. Distribution of the simulated response variable. One of:

  `"count"`

  :   (Default) Poisson-distributed counts.

  `"count_nb"`

  :   Negative binomial counts with dispersion `StDev`.

  `"occ"`

  :   Bernoulli (0/1) occupancy data via a logistic link.

  `"gaussian"`

  :   Normally distributed continuous response with standard deviation
      `StDev`.

- StDev:

  Positive numeric. Dispersion parameter for `"count_nb"` (the `size`
  argument of `rnbinom`) or the standard deviation for `"gaussian"`
  responses. Default: `0.5`.

- n_points:

  Positive integer or a `SpatVector`/`sf` point object. When an integer,
  that many points are placed on a hexagonal grid covering the raster
  extent (trimmed to 85% of the x/y range) and then randomly subsampled.
  When a spatial object is supplied, those points are used directly.
  Default: `50`.

- min_D:

  Positive numeric. Minimum inter-point spacing on the hexagonal grid.
  If `NULL` (default), set automatically to `1.55 * max(sigma)`.

- raster_stack:

  A `SpatRaster` object. Layer names become covariate names in the
  returned data frame.

- sigma:

  Positive numeric vector. True kernel scale parameters, one per raster
  layer. These are the values that
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  should recover. Must be in the same units as the raster projection.

- shape:

  Positive numeric vector. Shape parameters for the exponential power
  kernel (`kernel = "expow"`), one per raster layer. Values between 1
  and 50 are typical. Required when `kernel = "expow"`.

- max_D:

  Positive numeric. Maximum buffer radius for
  [`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  during simulation. If `NULL` (default), set automatically to 110% of
  the one-dimensional distance enclosing 99% of kernel weight along a
  line at `max(sigma)`. This historical default is smaller than the
  radius enclosing 99% of weight on a plane; supply a larger `max_D` if
  that planar cutoff is needed.

- user_seed:

  Optional integer seed for reproducibility. Default: `NULL`.

- ...:

  Additional arguments. Not currently used.

## Value

A named list with three elements:

- `obs`:

  Numeric vector of length `n_points` containing the simulated response
  values.

- `df`:

  Data frame with `n_points` rows and one column per raster layer plus a
  `y` column. The raster columns contain the true kernel-weighted mean
  values, scaled to zero mean and unit variance. This data frame can be
  used to fit a model for
  [`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

- `pts`:

  An `sf` POINT object with `n_points` rows. An `obs` column is appended
  containing the simulated response values. Pass this to
  [`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  as the `pts` argument.

## Details

Points are distributed on a hexagonal grid across the interior of the
raster extent (85% of the range in each direction to avoid edge
effects), then randomly subsampled to `n_points`. The `min_D` spacing
controls the grid resolution; if the grid produces fewer points than
`n_points`, `min_D` is reduced iteratively by 3% until enough points are
generated.

Kernel-weighted covariate values are computed using
[`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
at the specified `sigma` (and `shape`) values. These represent the
"true" covariate values that the optimization should recover.

## Examples

``` r
rs <- sim_rast()
rs <- terra::subset(rs, c(1,3))
s_dat <- sim_dat(alpha = 0.5,
                 beta = c(0.75,-0.75),
                 kernel = 'gaussian',
                 sigma = c(75, 150),
                 type = 'count',
                 raster_stack = rs,
                 max_D = 400)

plot(s_dat$df$y ~ s_dat$df$bin1)

plot(s_dat$df$y ~ s_dat$df$cont1)

```
