# Distance at Cumulative Kernel Proportion

Compute the distance at which a given cumulative proportion of kernel
weight is reached for several kernel types, measured either along a line
through the point (one-dimensional) or as the radius of a circle around
it (two-dimensional).

## Usage

``` r
k_dist(
  sigma,
  prob = 0.95,
  kernel = c("gaussian", "exp", "expow", "fixed"),
  beta = NULL,
  dimension = c("2d", "1d")
)
```

## Arguments

- sigma:

  Numeric. Scale parameter. For Gaussian and exponential, this is
  standard deviation or decay length. For expow, this is the kernel
  bandwidth.

- prob:

  Numeric. Desired cumulative proportion (e.g., 0.95).

- kernel:

  Character. One of "gaussian", "exp", "expow", or "fixed".

- beta:

  Numeric. Shape parameter for exponential power kernel. Ignored unless
  kernel = "expow".

- dimension:

  Character. `"2d"` (default) for the radius of the circle holding
  `prob` of ideal planar weight; `"1d"` for the original distance along
  a line through the point.

## Value

Numeric distance in map units at which the cumulative kernel weight
reaches `prob` under the selected dimensional definition.
