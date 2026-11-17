# Fitted Distance-Removal Model with a Fixed-Range Abundance Field

The class of object returned by
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md).
It extends the `unmarkedFitGDR` class of `unmarked` with one slot,
`response_field`, which holds the field geometry, the covariance
retained from the fit and its diagnostics. Most `unmarked` methods work
as usual. The methods listed below are replaced or stopped because the
inherited versions do not handle the correlated field.

## Usage

``` r
# S4 method for class 'unmarkedFitGDRField'
confint(object, parm, level = 0.95, type, method = c("normal", "profile"))

# S4 method for class 'unmarkedFitGDRField'
update(object, ...)

# S4 method for class 'unmarkedFitGDRField'
simulate(object, nsim = 1, seed = NULL, ...)

# S4 method for class 'unmarkedFitGDRField'
profile(fitted, ...)

# S4 method for class 'unmarkedFitGDRField'
ranef(object, ...)
```

## Arguments

- object, fitted:

  A fitted `unmarkedFitGDRField` object.

- parm, type, level, method:

  Arguments of [`confint()`](https://rdrr.io/r/stats/confint.html).
  `type` names the component (`"lambda"`, `"dist"` or `"rem"`), `parm`
  optionally selects coefficients within it, and `level` is the
  confidence level (default 0.95). `method = "normal"` (the default)
  gives Wald intervals from the retained covariance.
  `method = "profile"` stops for this class.

- ...:

  Not used. The stopped methods ignore extra arguments.

- nsim, seed:

  Not used. [`simulate()`](https://rdrr.io/r/stats/simulate.html) stops
  for this class.

## Slots

- `response_field`:

  A named list. The elements you are most likely to read are:

  - `definition`: the field definition, with `covariance`
    (`"exponential"`), the supplied `range`, `range_estimated` (always
    `FALSE`), the linear `units` of the CRS, `state_count` (the number
    of blocks) and `survey_count`.

  - `coords`: the survey coordinates supplied to
    [`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md),
    one row per survey. `block` gives each survey's block number, and
    `crs` holds the CRS as WKT text.

  - `point`: the complete parameter vector at the optimum, named by TMB
    block. It holds the fixed coefficients (`beta_*`, link scale), the
    log SDs (`lsigma_*`) and the block deviations in their internal
    whitened form (`b_lambda`).

  - `joint_covariance`: the covariance matrix of `point`, the inverse of
    the joint precision returned by
    [`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html).

  - `lower`: the lower Cholesky factor of the block correlation matrix.
    Multiplying the whitened deviations by `lower` gives the physical
    field values reported by
    [`coef()`](https://rdrr.io/r/stats/coef.html).

  - `diagnostics`: the checks of the selected start, such as `pdHess`,
    `maximum_absolute_gradient`,
    `conditional_latent_maximum_absolute_gradient` and
    `conditional_latent_Newton_decrement`.

  - `settings`: the evaluated `tmb_settings`. `audit_file` is the path
    of the final audit record, and `max_memory_mb` is the planning
    budget used by
    [`predict()`](https://rspatial.github.io/terra/reference/predict.html).

  The remaining elements are frozen design recipes and SHA-256
  fingerprints. Every method checks them before use and stops if any has
  changed, so do not edit this slot.

## Methods for this class

|  |  |
|----|----|
| Method | Behavior |
| [`coef()`](https://rdrr.io/r/stats/coef.html), [`logLik()`](https://rdrr.io/r/stats/logLik.html), [`unmarked::sigma()`](https://rdrr.io/r/stats/sigma.html), [`summary()`](https://rspatial.github.io/terra/reference/summary.html) | Inherited from `unmarked`. Field values reported by `coef(fixedOnly = FALSE)` are physical deviations on the log density scale |
| [`predict()`](https://rspatial.github.io/terra/reference/predict.html) | Replaced. See [predict,unmarkedFitGDRField-method](https://wpeterman.github.io/multiScaleR/reference/predict-unmarkedFitGDRField-method.md) |
| [`vcov()`](https://rdrr.io/r/stats/vcov.html) | Replaced. See [vcov,unmarkedFitGDRField-method](https://wpeterman.github.io/multiScaleR/reference/vcov-unmarkedFitGDRField-method.md) |
| [`confint()`](https://rdrr.io/r/stats/confint.html) | Normal (Wald) intervals from the retained covariance. `method = "profile"` stops |
| [`update()`](https://rspatial.github.io/terra/reference/update.html) | Stops. Refit through [`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md), which keeps the field |
| [`simulate()`](https://rdrr.io/r/stats/simulate.html), [`profile()`](https://rdrr.io/r/stats/profile.html), `ranef()` | Stop. They are not available for this class |

See
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md)
for a table that locates each parameter and its scale, and
[`vignette("v12-distance-removal", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v12-distance-removal.md)
for a worked example.

## See also

[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md),
[`make_gdr_refitter()`](https://wpeterman.github.io/multiScaleR/reference/make_gdr_refitter.md),
[`response_field_geometry()`](https://wpeterman.github.io/multiScaleR/reference/response_field_geometry.md)

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

  # TRUE: the methods of unmarked still apply
  print(methods::is(field_fit, "unmarkedFitGDR"))
  # The supplied range, its units, and the number of blocks and surveys
  str(field_fit@response_field$definition[c("range", "units",
                                            "state_count", "survey_count")])
  # Wald limits for the abundance coefficients, on the log scale
  print(confint(field_fit, type = "lambda"))
}
#> [1] TRUE
#> List of 4
#>  $ range       : num 2000
#>  $ units       : chr "metre"
#>  $ state_count : int 8
#>  $ survey_count: int 48
#>                 0.025     0.975
#> lambda(Int) 5.6993885 6.5593710
#> lambda(hab) 0.2067472 0.6779487
# }
```
