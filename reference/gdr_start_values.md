# List the Starting-Value Order of a Distance-Removal Model

Read the fitted coefficient and log SD values of an `unmarked`
distance-removal (GDR) model in the order that
[`unmarked::gdistremoval()`](https://ecoverseR.github.io/unmarked/reference/gdistremoval.html)
expects for its `starts` argument, with a map that names the model
component and position of each value. Use it to build extra starting
values for
[`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md)
and
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md).

## Usage

``` r
gdr_start_values(model)
```

## Arguments

- model:

  An `unmarkedFitGDR` fitted by
  [`unmarked::gdistremoval()`](https://ecoverseR.github.io/unmarked/reference/gdistremoval.html)
  with `unmarked` 1.5.2, using the C or TMB engine. Its stored estimates
  and optimizer values must be complete and agree with each other.
  Models with [`offset()`](https://rdrr.io/r/stats/offset.html) terms
  and zero-inflated TMB models are not supported. Finite values do not
  show that the model converged.

## Value

A plain list:

- `starts`: an unnamed numeric vector of the fitted values in starting
  order, on the link scales (log, logit, and log for SDs). Random-effect
  values are not included.

- `map`: a data frame with one row per element of `starts`:
  `constructor_position`, the TMB parameter `block` (such as
  `"beta_lambda"` or `"lsigma_lambda"`), `within_block` (position in the
  block), `fitted_position` (position in the fitted optimizer vector),
  `component` (`"lambda"`, `"dist"`, `"rem"` and so on), `kind`
  (`"regression"` or `"log_sd"`) and `label`, the coefficient or
  variance name stored by `unmarked`. For a TMB-engine model, paste
  `block` and `within_block` as `"beta_lambda[1]"` to name a start for
  [`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md).

- `fitted_values`: the fitted optimizer vector with its own names and
  order. `coordinate_count` is its length, not a count for AIC.
  `random_modes_omitted` gives the number of random-effect values left
  out in each component.

- `optimizer`: the stored objective `value` and `convergence` code, as
  reported by the fit.

- `provider`: the `unmarked` version, the engine and, for TMB, the TMB
  version. `format_version` identifies the layout of this list.

- `fingerprint`: a SHA-256 of this list. It does not identify the data,
  settings or likelihood.

- `fit_eligibility_established`, `likelihood_evaluated` and
  `uncertainty_evaluated`: always `FALSE`. The function only reads
  stored values.

## Why the order matters

With the TMB engine, the order of the `starts` argument can differ from
the order of the fitted optimizer vector. For example, `unmarked` may
place the negative binomial dispersion after the regression and variance
blocks in the starts, but next to abundance in the fitted vector.
Repeated coefficient labels cannot resolve this, so the map identifies
each value by its block and position and checks every value against the
stored estimates.

## What the function does not do

It reads stored values only. It does not evaluate the likelihood, its
gradient or Hessian, and it leaves the model unchanged. A failed or
nonconverged fit can still have finite values, and its stored
diagnostics keep their meaning. Nonfinite or inconsistent values stop
with an error rather than being dropped or repaired.

The list can be saved with
[`saveRDS()`](https://rspatial.github.io/terra/reference/serialize.html)
without saving a fitted model. It cannot restore a fit. Refitting still
needs the original data, transformations, formulas and optimizer
settings.

## See also

[`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md),
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md),
[`expand_gdr_intervals()`](https://wpeterman.github.io/multiScaleR/reference/expand_gdr_intervals.md)

## Examples

``` r
frame <- unmarked::unmarkedFrameGDR(matrix(c(0, 0, 1, 0), 2, byrow = TRUE),
  matrix(c(0, 0, 0, 1, 0, 0), 2, byrow = TRUE),
  dist.breaks = c(0, 25, 75), unitsIn = "m", period.lengths = c(1, 2, 1))
# A fixed point for inspecting coordinates, not a converged fit.
point <- unmarked::gdistremoval(data = frame, engine = "C", K = 30,
  starts = c(log(8), log(40), qlogis(.3)), se = FALSE,
  control = list(maxit = 0))
snapshot <- gdr_start_values(point)
# Log density, log distance scale and logit removal, in starting order
snapshot$starts
#> [1]  2.0794415  3.6888795 -0.8472979
# Component, block position and stored label of each value
snapshot$map[, c("block", "within_block", "component", "kind", "label")]
#>         block within_block component       kind       label
#> 1 beta_lambda            1    lambda regression (Intercept)
#> 2   beta_dist            1      dist regression (Intercept)
#> 3    beta_rem            1       rem regression (Intercept)
```
