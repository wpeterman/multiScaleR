# Simulate Distance-Removal Surveys from a Fitted Model

Draw new distance and removal count histories from a fitted
[`unmarked::gdistremoval()`](https://ecoverseR.github.io/unmarked/reference/gdistremoval.html)
model, holding its coefficients and kernel scales at their estimates.
Each visit uses its own detection probabilities, and the distance and
removal histories of a visit describe the same simulated individuals.
Use the replicates to study how the model behaves or to build
conditional predictive checks.

## Usage

``` r
simulate_gdistremoval(
  object,
  nsim = 1,
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

- nsim:

  Positive whole number of replicate frames. Default `1`.

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

A list of `nsim` `unmarkedFrameGDR` objects with the covariates and
protocol of the prediction frame and simulated counts. A visit with a
missing distance or removal history in the frame stays missing in both.
The `prediction_contract` attribute records that coefficients and scales
are fixed and that latent abundance and observation noise are simulated.

## Details

This function is experimental. Its arguments and results may change in a
later release without a deprecation period.

Latent abundance is drawn once per site and shared by its visits:
Poisson for `"P"`, negative binomial with the fitted size for `"NB"`,
and a structural zero with the fitted probability followed by a Poisson
draw for `"ZIP"`. The number counted at each visit is binomial with the
`detection` probability from
[`predict_gdistremoval()`](https://wpeterman.github.io/multiScaleR/reference/predict_gdistremoval.md).
That total is then split independently among the distance bins and among
the removal intervals with the fitted cell probabilities.

Call [`set.seed()`](https://rdrr.io/r/base/Random.html) first for
reproducible draws. Group effects follow `re.form`; new group effects
are not drawn. These conditional replicates alone are not a calibrated
bootstrap of random-effect variances, selected scales, or model
selection, which would need the whole fitting and sampling process
repeated.

[`unmarked::simulate()`](https://rdrr.io/r/stats/simulate.html) is not
used. In unmarked 1.5.2 it applies the first visit's distance parameter
to every visit and back-transforms the zero-inflation parameter with
[`exp()`](https://rdrr.io/r/base/Log.html), whereas the likelihood uses
[`plogis()`](https://rdrr.io/r/stats/Logistic.html).

## See also

[`predict_gdistremoval()`](https://wpeterman.github.io/multiScaleR/reference/predict_gdistremoval.md),
[`sim_dat_gdr()`](https://wpeterman.github.io/multiScaleR/reference/sim_dat_gdr.md)

## Examples

``` r
y <- matrix(c(2, 1, 0, 1, 1, 0, 3, 1, 1, 2, 0, 0), ncol = 3, byrow = TRUE)
frame <- unmarked::unmarkedFrameGDR(y, y, dist.breaks = c(0, 25, 50, 100),
  unitsIn = "m")
model <- unmarked::gdistremoval(data = frame, K = 25, se = FALSE)
set.seed(14)
replicate <- simulate_gdistremoval(model)[[1]]
# Both histories count the same simulated individuals
rowSums(replicate@yDistance) == rowSums(replicate@yRemoval)
#> [1] TRUE TRUE TRUE TRUE
```
