# Simulate Paired Distance-Removal Surveys with Retained Truth

Simulate single-visit point counts in which each bird is recorded by
distance band and by the time interval of its first detection. Abundance
is Poisson, detection by distance follows a half-normal curve, and
detection over time follows a removal (first-detection) process. The two
count tables describe the same birds. The observed data and the true
values used to generate them are returned separately, so you can fit a
model to the data and compare its estimates with the truth.

## Usage

``` r
sim_dat_gdr(
  density,
  distance_scale,
  removal,
  dist_breaks,
  interval_lengths,
  interval_units,
  units_in = c("m", "km"),
  site_covs = NULL,
  event_keys = NULL,
  seed,
  return_joint = FALSE,
  provenance = list()
)

# S3 method for class 'msr_gdr_simulation'
print(x, ...)
```

## Arguments

- density:

  Numeric vector of expected bird densities, in birds per square
  kilometer, one value per survey point. Zero is allowed. This is the
  response scale, not the log scale. Compute any ecological formula or
  kernel-weighted covariate effect first, for example
  `exp(log(400) + 0.5 * habitat)`.

- distance_scale:

  Positive number, or one per point. The standard deviation (sigma) of
  the half-normal detection function, in `units_in`. Detection falls to
  about 61% of its maximum at this distance. It is unrelated to the
  kernel scale of effect.

- removal:

  Probability, between 0 and 1, that a bird present is first detected
  during one unit of interval time. Supply one number, one value per
  point, or a matrix with points in rows and removal intervals in
  columns for interval-specific probabilities.

- dist_breaks:

  Increasing numeric distance boundaries starting at zero, in
  `units_in`. The last boundary is the survey radius.

- interval_lengths:

  Positive whole numbers giving the duration of each removal interval,
  in `interval_units`. Intervals may differ in length.

- interval_units:

  Character string naming the time unit of `interval_lengths`, for
  example `"minutes"`. It is recorded, not converted.

- units_in:

  Distance units, `"m"` (default) or `"km"`.

- site_covs:

  Optional data frame of site covariates, one row per point, such as
  kernel-weighted landscape covariates from
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).
  Columns are copied into the observed frame unchanged. As in
  `unmarked`, character columns are converted to factors with a warning.

- event_keys:

  Optional data frame with one row per point holding the original survey
  identifiers. It is returned unchanged, separate from the covariates,
  for use with
  [`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md).

- seed:

  Nonnegative whole number used as the random seed. Your random number
  state is restored afterward, so the simulation does not change later
  random draws in your session.

- return_joint:

  Logical, default `FALSE`. If `TRUE`, also return the simulated counts
  by point, distance band and interval as `joint`.

- provenance:

  Optional named list of notes about how the inputs were made, such as
  the generating formula, coefficients or the kernel scale used for a
  covariate. It is stored in `design$provenance` without being
  evaluated.

- x:

  An object returned by `sim_dat_gdr()`.

- ...:

  Not used.

## Value

A list of class `msr_gdr_simulation`:

- `observed`: the `unmarkedFrameGDR` to analyze, with the distance and
  removal tables, `site_covs`, distance breaks, units and interval
  lengths. It holds no true values.

- `truth`: a data frame with one row per point. See "Read the truth".

- `probabilities`: a list with `distance`, a points by distance-band
  matrix of the probability that a bird in the circle is detected in
  each band, and `removal`, a points by interval matrix of the
  probability that a bird is first detected in each interval.

- `joint`: the points by distance band by interval array of simulated
  counts when `return_joint = TRUE`, otherwise `NULL`.

- `event_keys`: as supplied.

- `design`: the settings used, including `density_units`
  (`"birds/km2"`), `area_km2` (the survey circle), `distance_scale` and
  `removal` per point, the breaks, interval lengths and units, the
  `seed` and `RNGkind`, the model assumptions (`mixture = "P"`,
  `keyfun = "halfnorm"`, one primary period, closed population), the
  `unmarked` version and `provenance`.

Printing shows a short summary of the number of points, bands and
intervals, the density range and the units, and returns the object
invisibly. The full frame and truth stay in the object.

## How the counts are generated

The survey area is `pi * radius^2`, in square kilometers, and expected
abundance is density times that area. Detection by distance and by time
are independent for each bird. The probability of a first detection in
an interval accounts for not being detected in earlier intervals, and an
interval of duration `t` has detection probability `1 - (1 - p)^t` for
per-unit probability `p`. The function draws independent Poisson counts
for every combination of distance band and interval, plus one for birds
never detected. Their sum is the Poisson abundance in the circle, and
both count tables are totals of the same simulated counts, so their row
totals always match.

## Read the truth

|  |  |  |
|----|----|----|
| Quantity | Where | How to read it |
| Expected density | `truth$expected_density` | Birds per square km, before imperfect detection |
| Expected abundance | `truth$expected_abundance` | Expected birds in the survey circle: density times `design$area_km2` |
| Realized abundance | `truth$realized_abundance` | Birds simulated in the circle, a Poisson draw around the expectation |
| Distance detection | `truth$distance_probability` | Chance that a bird in the circle is detected by distance, averaged over the circle |
| Removal detection | `truth$removal_probability` | Chance that a bird is first detected at some point during the count |
| Overall detection | `truth$observation_probability` | The product of the two, the share of birds an observer records |
| Expected observations | `truth$expected_observed` | Expected count after both detection processes |
| Realized observations | `truth$realized_observed` | Birds recorded, equal to the row totals of both count tables |
| Distance scale | `design$distance_scale` | Half-normal sigma, in `units_in` |
| Your notes | `design$provenance` | Whatever you passed in `provenance` |

The truth never enters the observed frame, so it cannot leak into fitted
covariates. `truth$site_row` gives the frame row of each point.

## What the simulator covers

The simulator covers single-visit Poisson surveys with a half-normal
distance function and equal or unequal whole-number interval lengths. It
assumes a closed population inside the radius. It does not simulate
repeat visits with availability, negative binomial or zero-inflated
abundance, movement across the radius, or random variation in detection.
It fits no model and leaves
[`sim_dat_unmarked()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_unmarked.md)
unchanged. See the "Simulate a survey with a known scale of effect"
section of
[`vignette("v12-distance-removal", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v12-distance-removal.md)
for a full example.

## See also

[`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md),
[`prepare_gdr_histories()`](https://wpeterman.github.io/multiScaleR/reference/prepare_gdr_histories.md),
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md)

## Examples

``` r
# Three points: two with birds and one with density zero
sim <- sim_dat_gdr(density = c(100, 200, 0), distance_scale = 40,
  removal = 0.3, dist_breaks = c(0, 25, 75),
  interval_lengths = c(1, 2, 1), interval_units = "minutes", seed = 24)
sim                        # compact summary
#> Single-primary Poisson distance-removal simulation
#> 3 survey events; 2 distance bins; 3 removal intervals
#> Expected density:   0 to 200 birds/km2 
#> Distance units: m ; removal time units: minutes 
#> Full observations: $observed; generating truth: $truth
# Truth behind the counts: about 1.8 and 3.5 birds expected in the circle
sim$truth[, c("expected_density", "expected_abundance",
              "realized_abundance", "realized_observed")]
#>   expected_density expected_abundance realized_abundance realized_observed
#> 1              100           1.767146                  1                 0
#> 2              200           3.534292                  5                 2
#> 3                0           0.000000                  0                 0
# The same birds by distance band and by time interval
sim$observed@yDistance
#>      [,1] [,2]
#> [1,]    0    0
#> [2,]    1    1
#> [3,]    0    0
sim$observed@yRemoval
#>      [,1] [,2] [,3]
#> [1,]    0    0    0
#> [2,]    0    2    0
#> [3,]    0    0    0
```
