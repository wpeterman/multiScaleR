# Covariance matrices for an analytic multiscale fit

Return the covariance matrix of the estimated parameters from a
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
fit made with an analytic backend. The joint matrix covers the kernel
parameters and the regression coefficients together, and it is what
[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md)
uses for joint uncertainty. The conditional matrix covers the
coefficients with the kernel parameters held at their estimates.
[multiScaleR-prediction](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-prediction.md)
explains when each kind of uncertainty applies.

## Usage

``` r
# S3 method for class 'multiScaleR'
vcov(object, type = c("joint", "conditional"), ...)
```

## Arguments

- object:

  A `multiScaleR` fit from `backend = "profile_glm"` or
  `backend = "joint_nb2"`. A `backend = "joint_mixed"` fit returns only
  native conditional covariance; its scale standard errors are in
  `object$scale_est`. Other backends raise an error. Use
  `vcov(object$opt_mod)` for their coefficient covariance conditional on
  the selected scales.

- type:

  `"joint"` (the default) or `"conditional"`.

  - `"joint"` returns the observed-information covariance of the log
    kernel parameters, the regression coefficients, and the log
    dispersion when the family has one.

  - `"conditional"` returns the coefficient covariance with the kernel
    parameters held fixed and the dispersion estimated. At the exact
    Poisson boundary of a `"joint_nb2"` fit it also holds the dispersion
    at that boundary, and joint covariance is unavailable there.

- ...:

  Passed to native [`vcov()`](https://rdrr.io/r/stats/vcov.html) for
  `"joint_mixed"`. Otherwise unused.

## Value

A numeric covariance matrix with named rows and columns, as described in
Details. For `"joint_mixed"` with `type = "conditional"`, the native
glmmTMB covariance list, one matrix per model component, conditional on
the selected scales.

## Read the joint matrix

Rows and columns are named by parameter:

- `log_scale:<covariate>`: `log(sigma / u)`, the natural log of that
  covariate's sigma after dividing it by the internal distance unit
  `u = object$kernel_inputs$unit_conv`. Dividing by `u` shifts the log
  by a constant, so the variance equals the variance of log(sigma) in
  map units.

- `log_shape:<covariate>`: the log shape parameter, for
  `kernel = "expow"`.

- `beta:<term>`: a regression coefficient on the link scale, named after
  its model-matrix column, for example `beta:(Intercept)`. A
  standardized kernel-weighted covariate's coefficient describes one
  standard deviation of that covariate.

- A final dispersion row when the family has one: `log_theta` for a
  negative binomial model and `log_residual_sd` for a Gaussian model.

The square root of a diagonal element is a standard error on that
parameter's scale. For sigma, a 95% Wald interval in map units is
`sigma x exp(-1.96 x SE)` to `sigma x exp(1.96 x SE)`, where `SE` is the
standard error of `log_scale`. These are the limits that
`summary(object)` reports, and `object$scale_est$SE` is `sigma x SE`.
Off-diagonal elements show how the estimates trade off, and
[`stats::cov2cor()`](https://rdrr.io/r/stats/cor.html) turns them into
correlations.

Joint covariance is unavailable for fits at a parameter bound or when
the identifiability checks fail. In that case `type = "joint"` raises an
error that lists `object$uncertainty$reasons`, rather than reporting a
misleading Wald result. The conditional matrix remains available.

One exception applies to both backends. When kernel parameters at a
search bound are the only failed check, the joint matrix holds those
parameters fixed at their bounds. Their rows and columns are zero, and
the other rows and columns come from the observed information with the
bounded parameters removed. That reduced information must pass the
remaining checks. The matrix then carries a `conditional_on` attribute
that names the bounded parameters, as does
`object$uncertainty$joint_conditional_on`. A zero there means the
parameter was held fixed, not that it is known exactly, and its standard
error is missing in `object$scale_est`. The other standard errors are no
wider than they would be if the bounded scale were free to vary, so
treat them as conditional on that scale. A dispersion bound, every
kernel parameter at a bound, or any other failed check still withholds
joint covariance.

## Read the conditional matrix

The conditional matrix holds only the `beta:` rows. It answers how
uncertain the coefficients are if the scales of effect are exactly their
estimates, which is also what `vcov(object$opt_mod)` and
`summary(object$opt_mod)` describe.

## See also

[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md)
for predictions that use these matrices,
[multiScaleR-prediction](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-prediction.md),
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
and
[`vignette("v09-estimates-and-intervals", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v09-estimates-and-intervals.md).

## Examples

``` r
## Running example: counts at 100 survey points, two landscape layers
data("landscape_counts")
data("surv_pts")
pts <- terra::vect(surv_pts)
land <- terra::rast(system.file("extdata", "landscape.tif",
                                package = "multiScaleR"))[[c("land1", "land2")]]
inputs <- kernel_prep(pts, land, max_D = 1700, kernel = "gaussian",
                      bin = TRUE, verbose = FALSE)
dat <- data.frame(landscape_counts, inputs$kernel_dat)
start <- glm(counts ~ site + land1 + land2, family = poisson(), data = dat)
fit <- multiScale_optim(start, inputs, backend = "profile_glm",
                        n_cores = 1, verbose = FALSE)

## Joint covariance of the two log scales and four coefficients
V <- vcov(fit)
round(V, 4)
#>                  log_scale:land1 log_scale:land2 beta:(Intercept) beta:site
#> log_scale:land1           0.0193          0.0011          -0.0002    0.0002
#> log_scale:land2           0.0011          0.0090          -0.0005   -0.0014
#> beta:(Intercept)         -0.0002         -0.0005           0.0083   -0.0016
#> beta:site                 0.0002         -0.0014          -0.0016    0.0070
#> beta:land1               -0.0009         -0.0020           0.0030   -0.0015
#> beta:land2                0.0014          0.0020          -0.0034    0.0001
#>                  beta:land1 beta:land2
#> log_scale:land1     -0.0009     0.0014
#> log_scale:land2     -0.0020     0.0020
#> beta:(Intercept)     0.0030    -0.0034
#> beta:site           -0.0015     0.0001
#> beta:land1           0.0068    -0.0008
#> beta:land2          -0.0008     0.0058

## 95% Wald limits for sigma, in meters, from the log-scale standard errors.
## They match the limits in summary(fit).
se_log <- sqrt(diag(V))[c("log_scale:land1", "log_scale:land2")]
sigma <- stats::setNames(fit$scale_est$Mean, rownames(fit$scale_est))
round(cbind(sigma, lower = sigma * exp(-qnorm(0.975) * se_log),
            upper = sigma * exp(qnorm(0.975) * se_log)))
#>       sigma lower upper
#> land1   238   182   313
#> land2   500   415   602

## Correlations show how the scale and coefficient estimates trade off
round(cov2cor(V), 2)
#>                  log_scale:land1 log_scale:land2 beta:(Intercept) beta:site
#> log_scale:land1             1.00            0.08            -0.02      0.02
#> log_scale:land2             0.08            1.00            -0.06     -0.17
#> beta:(Intercept)           -0.02           -0.06             1.00     -0.20
#> beta:site                   0.02           -0.17            -0.20      1.00
#> beta:land1                 -0.08           -0.26             0.40     -0.21
#> beta:land2                  0.13            0.28            -0.49      0.02
#>                  beta:land1 beta:land2
#> log_scale:land1       -0.08       0.13
#> log_scale:land2       -0.26       0.28
#> beta:(Intercept)       0.40      -0.49
#> beta:site             -0.21       0.02
#> beta:land1             1.00      -0.12
#> beta:land2            -0.12       1.00

## Coefficient standard errors with and without scale uncertainty
round(cbind(joint = sqrt(diag(V))[-(1:2)],
            conditional = sqrt(diag(vcov(fit, type = "conditional")))), 3)
#>                  joint conditional
#> beta:(Intercept) 0.091       0.091
#> beta:site        0.084       0.082
#> beta:land1       0.082       0.079
#> beta:land2       0.076       0.072
```
