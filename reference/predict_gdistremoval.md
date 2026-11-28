# Predict Whole-Survey Detection and Expected Counts for Distance-Removal Models

Combine the fitted components of an
[`unmarked::gdistremoval()`](https://ecoverseR.github.io/unmarked/reference/gdistremoval.html)
model into the quantities a survey produces: expected abundance and
density in the survey circle, the probability that an individual present
is detected at each visit, and the expected number of individuals
counted. Fitted coefficients and kernel scales are held at their
estimates.

## Usage

``` r
predict_gdistremoval(
  object,
  newdata = NULL,
  kernel_inputs = NULL,
  join_by = NULL,
  include_fixed = FALSE,
  re.form = NULL
)
```

## Arguments

- object:

  A fitted `unmarkedFitGDR` model, or a `multiScaleR` fit from
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  that contains one.

- newdata:

  An `unmarkedFrameGDR` describing the prediction sites, their visits,
  covariates, and removal interval durations. The default `NULL` uses
  the fitted frame. For new sites, supply response matrices of zeros
  with the same layout. Distance breaks and their units must match the
  fitted survey footprint. Removal durations must be positive integers
  and may differ from the fitted ones, so you can predict a longer or
  shorter count.

- kernel_inputs:

  Optional
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  output for the prediction sites. With a `multiScaleR` fit and explicit
  `newdata`, the spatial site covariates are rebuilt at the fitted
  scales with the training normalization, as in
  [`prepare_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/prepare_unmarked_data.md).
  Without it, the covariates in `newdata` are used as they are.

- join_by:

  Data frame of prediction-site keys, required with `kernel_inputs`. See
  [`prepare_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/prepare_unmarked_data.md).

- include_fixed:

  Logical. Default `FALSE`. With `kernel_inputs`, set `TRUE` to rebuild
  and check fixed-radius metrics in the model.

- re.form:

  Random-effect policy passed to
  [`unmarked::predict()`](https://rdrr.io/r/stats/predict.html). The
  default `NULL` uses the fitted group effects. `NA` sets group
  deviations to zero, for example at sites in unsampled groups. Neither
  averages over the random-effect distribution.

## Value

A data frame with one row per site and visit (primary period), in site
order and then visit order, with columns:

- `site`, `period`: row number in the frame and visit number.

- `lambda`: the native abundance or density parameter. For a
  zero-inflated (`"ZIP"`) model it is the Poisson-component mean, before
  structural zeros.

- `zero_probability`: the probability of a structural zero, 0 for `"P"`
  and `"NB"` models.

- `abundance`: expected number of individuals in the survey circle,
  including structural zeros. When phi is modeled it counts the
  individuals that use the site across visits, not those present at one
  visit.

- `density`: `abundance` divided by the survey-circle area, in
  individuals per hectare or square kilometer (the fit's `unitsOut`).

- `phi`: probability that an individual is available at this visit, 1
  for a single-visit model.

- `availability`: probability that an available individual is detected
  (gives a cue) in at least one removal interval.

- `perceptibility`: probability that an individual giving a cue anywhere
  in the circle is detected, from the distance detection function.

- `detection`: `phi * availability * perceptibility`, the probability
  that an individual present is counted at this visit.

- `expected_count`: `abundance * detection`, the expected number of
  individuals counted at this visit.

Attributes `density_units`, `survey_area` (in those units), and
`lambda_output` record the units and whether the native lambda is an
abundance or a density. The `prediction_contract` attribute records what
the predictions condition on (see
[prediction_contract](https://wpeterman.github.io/multiScaleR/reference/prediction_contract.md)).
These are point estimates without standard errors or intervals.

## Details

This function is experimental. Its arguments and results may change in a
later release without a deprecation period.

The components combine as they do in the native likelihood. The removal
probability refers to one time unit, so an interval lasting `d` units
has a first-detection probability of `1 - (1 - p)^d`. Distance and
removal detection are independent given the covariates. Each visit uses
its own distance parameter. The half-normal and exponential keys use
closed-form integrals over the distance bins, and the hazard-rate key
uses numerical integration.

[`unmarked::getP()`](https://ecoverseR.github.io/unmarked/reference/getP.html)
and functions built on it, such as
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`simulate()`](https://rdrr.io/r/stats/simulate.html), are not used. In
unmarked 1.5.2 they apply the first visit's distance parameter to every
visit, which matters when the distance formula contains a visit-level
(yearly site) covariate.

For a density map, `predict(fit, raster, type = "lambda")` remains the
route (see
[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md)).
For a zero-inflated model that map shows the Poisson-component mean. Use
this function on a frame of grid sites when you need the mean including
structural zeros, or the expected count under a planned survey protocol.

Zero-inflated models fitted with `engine = "TMB"`, models with active
offsets, and GDR response fields from
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md)
are not supported.

## See also

[`simulate_gdistremoval()`](https://wpeterman.github.io/multiScaleR/reference/simulate_gdistremoval.md),
[`prepare_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/prepare_unmarked_data.md),
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`sim_dat_gdr()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_gdr.md)

## Examples

``` r
y <- matrix(c(2, 1, 0, 1, 1, 0, 3, 1, 1, 2, 0, 0), ncol = 3, byrow = TRUE)
frame <- unmarked::unmarkedFrameGDR(y, y, dist.breaks = c(0, 25, 50, 100),
  unitsIn = "m", period.lengths = c(1, 2, 2))
model <- unmarked::gdistremoval(data = frame, K = 25, se = FALSE)
prediction <- predict_gdistremoval(model)
prediction[, c("density", "detection", "expected_count")]
#>    density detection expected_count
#> 1 6.043314 0.1166391       2.214467
#> 2 6.043314 0.1166391       2.214467
#> 3 6.043314 0.1166391       2.214467
#> 4 6.043314 0.1166391       2.214467
# density corrects for detection; expected_count is what a survey records
```
