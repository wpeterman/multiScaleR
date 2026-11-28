# Fit a Fixed-Range Spatial Abundance Field to a Distance-Removal Model

Replace the independent random intercepts on abundance in a native TMB
distance-removal (GDR) model with a spatially correlated Gaussian field.
Nearby survey blocks then tend to have similar density deviations. The
correlation between two blocks falls with the distance `d` between them
as `exp(-d / range)`. You supply `range`; it is not estimated.
Ecological scales of effect remain separate parameters, estimated later
with
[`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md)
and
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

## Usage

``` r
fit_gdr_response_field(
  model,
  coords,
  crs,
  range,
  cold_starts,
  audit_dir,
  tmb_settings,
  max_memory_mb = 512
)
```

## Arguments

- model:

  An `unmarkedFitGDR` fitted by
  [`unmarked::gdistremoval()`](https://ecoverseR.github.io/unmarked/reference/gdistremoval.html)
  with `engine = "TMB"`. It serves as the template. Its data, formulas,
  cutoff `K` and estimates are reused, and its independent abundance
  intercepts become the field. See "Check that your model qualifies" for
  the supported model forms.

- coords:

  Numeric matrix or data frame with two projected coordinate columns
  (x, y) and one row per survey, in the row order of the model's frame.
  Row names must equal the frame's site row names, which are
  `rownames(unmarked::siteCovs(model@data))` when the frame has site
  covariates. All surveys in one block (one level of the abundance
  grouping factor) must carry exactly the same coordinates, usually the
  block center, and different blocks must have different coordinates.

- crs:

  The projected coordinate reference system of `coords`, as an EPSG code
  (for example `32617`) or any other input accepted by
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html).
  It must have two horizontal axes and declared linear units.
  Coordinates are never transformed, so `range` and `coords` must use
  these units.

- range:

  Positive number. The fixed correlation range of the exponential field,
  in the linear units of `crs`. Two blocks `range` apart have
  correlation `exp(-1)`, about 0.37, and blocks three ranges apart have
  correlation about 0.05. This is a distance over which density
  deviations are shared. It is not a kernel scale of effect, and it is
  not counted as a parameter.

- cold_starts:

  Named list of one or more extra starting vectors (for example
  `list(cold = cold)`). Each vector holds one value for every fixed
  coefficient and log standard deviation, on the link scales listed by
  [`gdr_start_values()`](https://wpeterman.github.io/multiScaleR/reference/gdr_start_values.md).
  Supply it either unnamed, in the order of
  `gdr_start_values(model)$starts`, or named by parameter block and
  position, such as `"beta_lambda[1]"` (log density intercept),
  `"beta_dist[1]"` (log distance scale), `"beta_rem[1]"` (logit removal
  probability) and `"lsigma_lambda[1]"` (log SD of the abundance
  intercepts). Each vector must differ from the template's estimates and
  from the other vectors. The names `"original"` and `"final"` are
  reserved.

- audit_dir:

  Character string naming a directory for audit records. It is created
  if needed. Every fit writes small `.rds` records there: the request,
  the validated inputs, each attempt (including failures), the optimizer
  result saved before covariance work, and the final outcome. Records
  hold numbers, settings and identities, not fitted models. Keep the
  directory to inspect a failed fit.

- tmb_settings:

  Named list of optimizer settings with elements `method`, `se`,
  `threads` and `control`. This function requires `method = "BFGS"`,
  `se = TRUE`, `threads = 1` and a positive integer `control$maxit`, for
  example
  `list(method = "BFGS", se = TRUE, threads = 1L, control = list(maxit = 2500L, reltol = 1e-10))`.
  The optional elements `random.start`, `inner.control` and
  `record_memory` are described in
  [`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md).
  Pass the identical list to
  [`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md)
  when you later estimate scales of effect with the field.

- max_memory_mb:

  Positive number, default 512. A conservative planning budget, in MiB,
  for the dense field covariance and the TMB tape. The fit stops before
  optimizing when its estimate exceeds the budget. The same budget
  limits later
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  calls on the fitted object. It is an estimate, not an operating-system
  memory cap.

## Value

An object of class
[unmarkedFitGDRField](https://wpeterman.github.io/multiScaleR/reference/unmarkedFitGDRField-class.md).
It extends the `unmarkedFitGDR` class of `unmarked`, so
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`logLik()`](https://rdrr.io/r/stats/logLik.html),
[`unmarked::sigma()`](https://rdrr.io/r/stats/sigma.html) and
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
work as usual, while
[`predict()`](https://rspatial.github.io/terra/reference/predict.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`confint()`](https://rdrr.io/r/stats/confint.html) use the retained
field covariance. The table under "Read the fitted model" says where
each parameter lives and on what scale. The log-likelihood is the native
Laplace approximation of the marginal likelihood, with the field
integrated out. The AIC counts the same parameters as the template
model, because the range is fixed. The attribute `msr_gdr_refit_audit`
summarizes the fit: `status`, the `attempts` from every start with their
diagnostics, the index of the `selected` start,
`maximum_start_value_difference` (the spread of the negative
log-likelihoods among starts) and the audit `path`.

## Details

A random intercept for each survey block lets density differ among
blocks beyond what the covariates explain. In the template model those
block deviations are independent. This function keeps the same distance
and removal observation model, the same cutoff `K` and the same SD of
the deviations, but makes the deviations of nearby blocks correlated.
See the "Model spatial correlation among survey blocks" section of
[`vignette("v12-distance-removal", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v12-distance-removal.md)
for a worked example.

## Check that your model qualifies

The template must have all of the following:

- Poisson abundance (`mixture = "P"`) and a half-normal distance
  function (`keyfun = "halfnorm"`).

- Density output in birds per square kilometer (`output = "density"`,
  `unitsOut = "kmsq"`).

- One visit per point (one primary period) and complete counts for every
  survey, with no sites dropped by `unmarked`.

- Distance bands in meters with boundaries `c(0, 25, 75)` or
  `c(0, 25, 50, 75)`, so the survey circle has a 75 m radius.

- Exactly one random intercept in the abundance formula, such as
  `~ land1 + (1 | block)`. Every level of the grouping factor must have
  at least one survey.

- Optionally, one independent random intercept in the distance formula,
  such as an observer effect. It stays independent. The removal formula
  takes no random effects.

- No [`offset()`](https://rdrr.io/r/stats/offset.html) terms.

- `unmarked` 1.5.2 and `TMB` 1.9.25.

The field correlation matrix must also be numerically positive definite.
Blocks at nearly identical locations relative to `range` can fail this
check. No numerical jitter is added.

## How the fit is checked

The field model is fitted from the template's own estimates (the start
named `"original"`) and from every vector in `cold_starts`. Each fit
starts its block deviations at zero. Every start must return a
converged, warning-free fit with usable covariance and a largest
absolute gradient below 0.001. All starts must also agree, with negative
log-likelihoods within 1e-5 and coefficients and log SDs within 0.001 of
one another. If any start fails or the starts disagree, the function
stops with an error and leaves the records in `audit_dir`. When every
check passes, it returns the fit with the highest likelihood. Agreement
among a few starts does not prove that no better optimum exists. If the
starts disagree, compare the attempts recorded in `audit_dir` before
trying other starting values.

## Read the fitted model

|  |  |  |  |
|----|----|----|----|
| Quantity | Where to find it | Scale | How to read it |
| Density intercept | `coef(fit)["lambda(Int)"]` | log | [`exp()`](https://rdrr.io/r/base/Log.html) gives density per square km when covariates and the field are zero |
| Abundance covariate effect | `coef(fit)["lambda(x)"]` | log | [`exp()`](https://rdrr.io/r/base/Log.html) is the factor by which density changes per unit of `x` |
| Distance scale | `coef(fit)["dist(Int)"]` | log | [`exp()`](https://rdrr.io/r/base/Log.html) gives the half-normal scale in meters |
| Removal probability | `coef(fit)["rem(Int)"]` | logit | [`plogis()`](https://rdrr.io/r/stats/Logistic.html) gives the probability, per unit of interval time, that a bird present is first detected |
| Field SD | `unmarked::sigma(fit)` | log density | Typical size of a block deviation. A block one SD above zero has `exp(SD)` times the density of a block with no deviation |
| Block deviations | `coef(fit, type = "lambda", fixedOnly = FALSE)` | log density | The entries labeled `lambda(b_lambda)`, one fitted field value per block, follow the fixed coefficients |
| Correlation range | `fit@response_field$definition$range` | CRS units | Supplied by you, not estimated |

Internally, TMB works with uncorrelated (whitened) versions of the block
deviations. Every value reported by
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html) and
[`predict()`](https://rspatial.github.io/terra/reference/predict.html)
is converted back to the physical field, the deviation at each block's
location.

## How uncertainty is reported

[vcov()](https://wpeterman.github.io/multiScaleR/reference/vcov-unmarkedFitGDRField-method.md)
returns the covariance retained from the fit. It includes the covariance
between components (abundance, distance, removal) and between the fixed
coefficients and the field values. It is computed once during fitting,
so [`vcov()`](https://rdrr.io/r/stats/vcov.html) does not rerun TMB.

[predict()](https://wpeterman.github.io/multiScaleR/reference/predict-unmarkedFitGDRField-method.md)
uses that covariance for its standard errors. At a block used in
fitting, the SE comes from the fitted coefficients and field values
alone. At a new location, the SE adds the field variance that remains
after conditioning on the fitted blocks. That extra variance describes
how uncertain the field is between blocks. It is not observation noise.
Density predictions are [`exp()`](https://rdrr.io/r/base/Log.html) of
the linear predictor at the fitted field values. They are not averages
over possible field values, which would be larger. With `re.form = NA`,
predictions set the field to zero, which describes a block with no
deviation rather than an average over blocks.

All of these intervals treat the ecological scales of effect and the
correlation range as known. They add no uncertainty for either.

## Estimate scales of effect with the field in place

Build a refit function from the field fit, using the same starts and the
identical `tmb_settings`, and pass it to
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md):

    refit <- make_gdr_refitter(field_fit, list(cold = cold),
      audit_dir = tempfile("gdr-field-refits-"), tmb_settings = settings)
    opt <- multiScale_optim(field_fit, kernel_inputs, refit_fn = refit,
      n_cores = 1)

Each candidate scale refits the field from every start, with the range
held fixed. Expect minutes even for 100 surveys. The search must run on
one core with the default `backend = "refit"` and `PSOCK = FALSE`. The
ordinary
[`update()`](https://rspatial.github.io/terra/reference/update.html)
method stops for this class. Use the refit function instead.

## Map density with the field

A density map needs the field value at every raster cell. Bind the
fitting coordinates to the scale-of-effect result with
[`response_field_geometry()`](https://wpeterman.github.io/multiScaleR/reference/response_field_geometry.md),
then pass the descriptor to
[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md):

    geometry <- response_field_geometry(opt, coords = block_xy, crs = 32617)
    layers <- kernel_scale.raster(land, multiScaleR = opt,
      scale_center = TRUE, match_buffer = TRUE)
    dens <- predict(opt, layers, type = "lambda", se.fit = TRUE,
      uncertainty = "conditional",
      raster_args = list(field_geometry = geometry))

The raster CRS must match `crs`, and prediction runs on one core. Each
cell gets the field value predicted at its center from the fitted block
values. The grouping factor of the field (`block` above) is not needed
as a raster. A factor that also enters a fixed part of the model still
needs an explicit value through `site_values` in
[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md)
or a raster coded with the fitted levels. To map the zero-field
reference instead, omit `field_geometry` and set `re.form = NA`.

## Not supported

Estimating the range, other abundance mixtures or distance functions,
repeat visits with availability, a removal field,
[`simulate()`](https://rdrr.io/r/stats/simulate.html),
[`profile()`](https://rdrr.io/r/stats/profile.html), `ranef()` and
spatial cross-validation are not available for this model.

## See also

[unmarkedFitGDRField](https://wpeterman.github.io/multiScaleR/reference/unmarkedFitGDRField-class.md),
[`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md),
[`gdr_start_values()`](https://wpeterman.github.io/multiScaleR/reference/gdr_start_values.md),
[`response_field_geometry()`](https://wpeterman.github.io/multiScaleR/reference/response_field_geometry.md),
[`vignette("v12-distance-removal", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v12-distance-removal.md)

## Examples

``` r
# \donttest{
# The field needs unmarked 1.5.2 and TMB 1.9.25.
if (requireNamespace("TMB", quietly = TRUE) &&
    utils::packageVersion("unmarked") == "1.5.2" &&
    utils::packageVersion("TMB") == "1.9.25") {
  # Eight survey blocks 1 km apart along a line, six point counts in each
  block <- factor(rep(paste0("b", 1:8), each = 6))
  hab <- rep(seq(-1, 1, length.out = 6), 8)
  shift <- c(-0.4, -0.3, -0.1, 0, 0.1, 0.2, 0.35, 0.4)[block]
  sim <- sim_dat_gdr(density = exp(log(400) + 0.5 * hab + shift),
    distance_scale = 40, removal = 0.2, dist_breaks = c(0, 25, 50, 75),
    interval_lengths = c(2, 3, 5), interval_units = "minutes",
    site_covs = data.frame(hab = hab, block = block), seed = 1)
  umf <- sim$observed

  # Template: an independent random intercept for each block (TMB engine)
  template <- unmarked::gdistremoval(lambdaformula = ~ hab + (1 | block),
    removalformula = ~ 1, distanceformula = ~ 1, data = umf,
    keyfun = "halfnorm", mixture = "P", output = "density",
    unitsOut = "kmsq", K = 60, engine = "TMB",
    starts = c(log(400), 0, log(40), qlogis(0.2), log(0.3)))

  # Every survey carries the projected coordinates (m) of its block center
  block_xy <- cbind(x = 500000 + 1000 * (as.integer(block) - 1),
                    y = 4400000)
  rownames(block_xy) <- rownames(unmarked::siteCovs(umf))

  # A second set of starting values, named by parameter block and position
  cold <- c("beta_lambda[1]" = log(300), "beta_lambda[2]" = 0,
            "beta_dist[1]" = log(30), "beta_rem[1]" = qlogis(0.4),
            "lsigma_lambda[1]" = log(0.5))
  settings <- list(method = "BFGS", se = TRUE, threads = 1L,
                   control = list(maxit = 2500L, reltol = 1e-10))

  field_fit <- fit_gdr_response_field(template, coords = block_xy,
    crs = 32617, range = 2000, cold_starts = list(cold = cold),
    audit_dir = tempfile("gdr-field-"), tmb_settings = settings)

  # Log density, log distance scale (m) and logit removal probability
  print(unmarked::coef(field_fit))
  # SD of the block deviations, on the log density scale
  print(unmarked::sigma(field_fit))
  # The range adds no parameter, so the AIC penalty matches the template
  print(c(template = template@AIC, field = field_fit@AIC))
}
#> lambda(Int) lambda(hab)   dist(Int)    rem(Int) 
#>   6.1293802   0.4423489   3.6600799  -1.6307798 
#>    Model Groups        Name     sigma      lower   upper
#> 1 lambda  block (Intercept) 0.2465782 0.09603057 0.63314
#> template    field 
#> 581.9212 580.3834 
# }
```
