# Covariance Matrix for a Distance-Removal Model with a Fixed-Range Field

Return the covariance matrix of the coefficients of a model fitted by
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md).
The matrix is retained from the fit, so this method does not rerun TMB.
It matches [`coef()`](https://rdrr.io/r/stats/coef.html) on the link
scale and keeps the covariances between components and between the fixed
coefficients and the field values.

## Usage

``` r
# S4 method for class 'unmarkedFitGDRField'
vcov(object, type, altNames = TRUE, method = "hessian", fixedOnly = TRUE, ...)
```

## Arguments

- object:

  An
  [unmarkedFitGDRField](https://wpeterman.github.io/multiScaleR/reference/unmarkedFitGDRField-class.md)
  returned by
  [`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md).

- type:

  Optional character naming one component: `"lambda"` (abundance),
  `"dist"` (distance) or `"rem"` (removal). If omitted, all components
  are returned together, in the order of `coef(object)`.

- altNames:

  Logical, default `TRUE`. `TRUE` labels rows and columns with component
  prefixes, such as `lambda(Int)`, as
  [`coef()`](https://rdrr.io/r/stats/coef.html) does.

- method:

  Character. Only `"hessian"`, the default, is available. It returns the
  covariance computed from the curvature of the likelihood during
  fitting. Bootstrap covariance is not available for this model.

- fixedOnly:

  Logical, default `TRUE`. `TRUE` returns the fixed regression
  coefficients only. `FALSE` adds the field value of every block
  (labeled `lambda(b_lambda)`) and any observer intercepts in the
  distance component.

- ...:

  Not used. Unrecognized arguments stop the call.

## Value

A square numeric matrix with row and column names matching
`coef(object, type, altNames = altNames, fixedOnly = fixedOnly)`.
Entries are on the link scale of each coefficient: log for abundance and
distance, logit for removal. Block field values are physical deviations
on the log density scale, not TMB's internal whitened values. The square
roots of the diagonal are the standard errors.

## Details

The log SDs of the random effects are estimated parameters, but they are
not rows of this matrix. Use `unmarked::sigma(object)` to see the field
SD with its confidence limits.

The matrix treats the ecological scales of effect and the correlation
range as known, and it contains no observation noise.
[`confint()`](https://rdrr.io/r/stats/confint.html) uses it for normal
(Wald) intervals, and
[predict()](https://wpeterman.github.io/multiScaleR/reference/predict-unmarkedFitGDRField-method.md)
uses it for standard errors.

## See also

[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md),
[unmarkedFitGDRField](https://wpeterman.github.io/multiScaleR/reference/unmarkedFitGDRField-class.md),
[predict,unmarkedFitGDRField-method](https://wpeterman.github.io/multiScaleR/reference/predict-unmarkedFitGDRField-method.md)

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

  # Abundance coefficients: covariance and standard errors (log scale)
  V <- vcov(field_fit, type = "lambda")
  print(V)
  print(sqrt(diag(V)))
  # All components together: abundance, distance, then removal
  print(round(vcov(field_fit), 4))
  # Adding the eight block field values gives a 10 x 10 matrix
  print(dim(vcov(field_fit, type = "lambda", fixedOnly = FALSE)))
}
#>              lambda(Int)  lambda(hab)
#> lambda(Int)  0.048130797 -0.002926743
#> lambda(hab) -0.002926743  0.014449645
#> lambda(Int) lambda(hab) 
#>   0.2193873   0.1202067 
#>             lambda(Int) lambda(hab) dist(Int) rem(Int)
#> lambda(Int)      0.0481     -0.0029   -0.0104  -0.0126
#> lambda(hab)     -0.0029      0.0144    0.0000   0.0000
#> dist(Int)       -0.0104      0.0000    0.0079   0.0000
#> rem(Int)        -0.0126      0.0000    0.0000   0.0384
#> [1] 10 10
# }
```
