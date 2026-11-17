# Predict from a Distance-Removal Model with a Fixed-Range Abundance Field

Predict density, the distance scale or the removal probability from a
model fitted by
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md).
Abundance predictions include the fitted field, at fitted blocks or at
new projected locations, and their standard errors include the field's
uncertainty.

## Usage

``` r
# S4 method for class 'unmarkedFitGDRField'
predict(
  object,
  type = "lambda",
  newdata = NULL,
  backTransform = TRUE,
  na.rm = TRUE,
  appendData = FALSE,
  level = 0.95,
  re.form = NULL,
  coords = NULL,
  crs = NULL,
  ...
)
```

## Arguments

- object:

  An
  [unmarkedFitGDRField](https://wpeterman.github.io/multiScaleR/reference/unmarkedFitGDRField-class.md)
  returned by
  [`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md).

- type:

  Character. The model component to predict: `"lambda"` (density, the
  default), `"dist"` (distance scale) or `"rem"` (removal probability).
  `"state"` is accepted as another name for `"lambda"`.

- newdata:

  A data frame of covariate values, or an `unmarkedFrameGDR` with the
  fitted survey protocol (one visit, the same distance bands and the
  same interval lengths). The default `NULL` predicts for the rows used
  in fitting. Transformations, factor levels and spline bases keep their
  definitions from fitting. With the default, removal predictions have
  one row per survey and removal interval.

- backTransform:

  Logical, default `TRUE`. `TRUE` returns density per square kilometer,
  the distance scale in meters or the removal probability per unit of
  interval time. `FALSE` returns the link scale (log, log or logit).

- na.rm:

  Logical, default `TRUE`. Kept for compatibility with `unmarked`. Rows
  with missing predictor values are not dropped. They stop the
  prediction, so remove them first.

- appendData:

  Logical, default `FALSE`. If `TRUE`, the prediction columns are added
  to the right of the prediction rows.

- level:

  Confidence level between 0 and 1, default 0.95. `NULL` skips the
  limits, which are then returned as `NA`.

- re.form:

  `NULL` (default) includes the fitted field, and any observer
  intercepts in the distance component. `NA` sets all of them to zero. A
  zero-field prediction describes a block with no deviation. It is not
  an average over blocks.

- coords:

  Optional two-column matrix or data frame of projected coordinates for
  abundance predictions at new locations, one row per row of `newdata`.
  The default `NULL` locates each row by its level of the grouping
  factor (for example `block = "b3"`), which must be a fitted level.

- crs:

  The CRS of `coords`, as an EPSG code or other input accepted by
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html).
  Required with `coords`. It must equal the CRS used in fitting, because
  coordinates are never transformed.

- ...:

  Not used. Unrecognized arguments stop the prediction.

## Value

A data frame with one row per prediction row:

- `Predicted`: the prediction on the scale set by `backTransform`.

- `SE`: its standard error on the same scale. Back-transformed SEs use
  the delta method.

- `lower`, `upper`: confidence limits at `level`, computed on the link
  scale and then back-transformed, so they are not symmetric around
  `Predicted`.

The attribute `response_field_uncertainty` explains the SE. Its element
`fitted_joint_predictor_variance` is the link-scale variance from the
fitted coefficients and field. `new_field_residual_predictor_variance`
is the extra link-scale variance at new locations (zero at fitted
blocks). `SE_includes_new_field_residual` is `TRUE` when any row has
that extra variance. `omitted` lists what the SE leaves out: observation
noise and uncertainty in the ecological scales and the range.

## Choose where to predict

For abundance, each prediction needs a field value. A row whose grouping
factor names a fitted block uses that block's fitted deviation. A row
with a location in `coords` uses the field value predicted at that
location from the fitted block deviations. At a fitted block location
this equals the block's own deviation. Far from all blocks, compared
with the range, it shrinks toward zero. For raster maps, use
[`response_field_geometry()`](https://wpeterman.github.io/multiScaleR/reference/response_field_geometry.md)
with
[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md)
instead of calling this method cell by cell.

For `type = "dist"` with an observer random intercept, rows must name
fitted observers unless `re.form = NA`.

## What the standard errors include

The SE combines two link-scale variances. The first comes from the joint
covariance of the fitted coefficients and field values, including their
cross-covariances. The second applies only at new locations. It is the
field variance left after conditioning on the fitted blocks, evaluated
at the estimated field SD. It is zero at a fitted block and grows toward
the full field variance far from all blocks. Neither variance is
observation noise, so these are confidence intervals for the expected
value, not prediction intervals for a new count.

A density prediction is [`exp()`](https://rdrr.io/r/base/Log.html) of
the linear predictor at the fitted field values. It is not the average
density over possible field values. The intervals treat the ecological
scales of effect and the correlation range as known.

## See also

[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md),
[`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md)

## Examples

``` r
# \donttest{
# The field needs unmarked 1.5.2 and TMB 1.9.25.
if (requireNamespace("TMB", quietly = TRUE) &&
    utils::packageVersion("unmarked") == "1.5.2" &&
    utils::packageVersion("TMB") == "1.9.25") {
  # Eight survey blocks 1 km apart, six point counts in each
  block <- factor(rep(paste0("b", 1:8), each = 6))
  hab <- rep(seq(-1, 1, length.out = 6), 8)
  shift <- c(-0.4, -0.3, -0.1, 0, 0.1, 0.2, 0.35, 0.4)[block]
  sim <- sim_dat_gdr(density = exp(log(400) + 0.5 * hab + shift),
    distance_scale = 40, removal = 0.2, dist_breaks = c(0, 25, 50, 75),
    interval_lengths = c(2, 3, 5), interval_units = "minutes",
    site_covs = data.frame(hab = hab, block = block), seed = 1)
  template <- unmarked::gdistremoval(lambdaformula = ~ hab + (1 | block),
    removalformula = ~ 1, distanceformula = ~ 1, data = sim$observed,
    keyfun = "halfnorm", mixture = "P", output = "density",
    unitsOut = "kmsq", K = 60, engine = "TMB",
    starts = c(log(400), 0, log(40), qlogis(0.2), log(0.3)))
  block_xy <- cbind(x = 500000 + 1000 * (as.integer(block) - 1),
                    y = 4400000)
  rownames(block_xy) <- rownames(unmarked::siteCovs(sim$observed))
  cold <- c("beta_lambda[1]" = log(300), "beta_lambda[2]" = 0,
            "beta_dist[1]" = log(30), "beta_rem[1]" = qlogis(0.4),
            "lsigma_lambda[1]" = log(0.5))
  settings <- list(method = "BFGS", se = TRUE, threads = 1L,
                   control = list(maxit = 2500L, reltol = 1e-10))
  field_fit <- fit_gdr_response_field(template, block_xy, crs = 32617,
    range = 2000, cold_starts = list(cold = cold),
    audit_dir = tempfile("gdr-field-"), tmb_settings = settings)

  # Density (birds per square km) at average habitat: in block b3, at a
  # point halfway between blocks b3 and b4, and with the field set to zero
  in_block <- predict(field_fit, type = "lambda",
    newdata = data.frame(hab = 0, block = factor("b3", levels(block))))
  between <- predict(field_fit, type = "lambda",
    newdata = data.frame(hab = 0), coords = cbind(502500, 4400000),
    crs = 32617)
  no_field <- predict(field_fit, type = "lambda",
    newdata = data.frame(hab = 0), re.form = NA)
  dens <- rbind(in_block, between, no_field)
  rownames(dens) <- c("block b3", "between b3 and b4", "no field")
  print(dens)
  # The SE between blocks adds field variance left after conditioning on
  # the fitted blocks (both variances on the log scale)
  u <- attr(between, "response_field_uncertainty")
  print(c(fitted = unname(u$fitted_joint_predictor_variance),
          between_blocks = unname(u$new_field_residual_predictor_variance)))
  # Removal probability per minute, one row per survey and interval
  print(head(predict(field_fit, type = "rem"), 3))
}
#>                   Predicted        SE    lower    upper
#> block b3           378.4863  85.53074 243.0495 589.3939
#> between b3 and b4  423.9838  94.82578 273.5107 657.2402
#> no field           459.1513 100.73197 298.6847 705.8276
#>         fitted between_blocks 
#>     0.03513006     0.01489118 
#>   Predicted         SE     lower     upper
#> 1 0.1637236 0.02684601 0.1176291 0.2233096
#> 2 0.1637236 0.02684601 0.1176291 0.2233096
#> 3 0.1637236 0.02684601 0.1176291 0.2233096
# }
```
