# Standard model methods for multiScaleR fits

The usual model accessors work on a fit from
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
Each reads the final model, the model refit at the estimated scales of
effect and stored in `object$opt_mod`. `coef(fit)`, `fitted(fit)`, and
the other methods therefore describe that model, never the starting
model passed to
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
Step-selection wrappers from
[`amt::fit_clogit()`](https://rdrr.io/pkg/amt/man/fit_clogit.html) and
[`amt::fit_issf()`](https://rdrr.io/pkg/amt/man/fit_clogit.html) are
first unwrapped to their native `survival` model.

A fit from
[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
has no single fitted model, so most of these methods stop with a pointer
to the Bayesian tools (see "Bayesian fits").

## Usage

``` r
# S3 method for class 'multiScaleR'
coef(object, ...)

# S3 method for class 'multiScaleR'
fitted(object, ...)

# S3 method for class 'multiScaleR'
residuals(object, ...)

# S3 method for class 'multiScaleR'
nobs(object, ...)

# S3 method for class 'multiScaleR'
formula(x, ...)

# S3 method for class 'multiScaleR'
family(object, ...)

# S3 method for class 'multiScaleR'
deviance(object, ...)

# S3 method for class 'multiScaleR'
df.residual(object, ...)

# S3 method for class 'multiScaleR'
logLik(object, ...)

# S3 method for class 'multiScaleR_bayes'
coef(object, ...)

# S3 method for class 'multiScaleR_bayes'
fitted(object, ...)

# S3 method for class 'multiScaleR_bayes'
residuals(object, ...)

# S3 method for class 'multiScaleR_bayes'
vcov(object, ...)

# S3 method for class 'multiScaleR_bayes'
logLik(object, ...)

# S3 method for class 'multiScaleR_bayes'
nobs(object, ...)
```

## Arguments

- object, x:

  A `multiScaleR` object from
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
  or for [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`vcov()`](https://rdrr.io/r/stats/vcov.html),
  [`logLik()`](https://rdrr.io/r/stats/logLik.html), and
  [`nobs()`](https://rdrr.io/r/stats/nobs.html) a `multiScaleR_bayes`
  object from
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md).

- ...:

  Passed to the native method of the final model. Ignored by the
  Bayesian methods.

## Value

- [`coef()`](https://rdrr.io/r/stats/coef.html),
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
  [`residuals()`](https://rdrr.io/r/stats/residuals.html),
  [`nobs()`](https://rdrr.io/r/stats/nobs.html),
  [`formula()`](https://rdrr.io/r/stats/formula.html),
  [`family()`](https://rdrr.io/r/stats/family.html),
  [`deviance()`](https://rdrr.io/r/stats/deviance.html), and
  [`df.residual()`](https://rdrr.io/r/stats/df.residual.html): the value
  the same generic returns for the final model.

- [`logLik()`](https://rdrr.io/r/stats/logLik.html): an object of class
  `"logLik"` holding the log-likelihood of the final model, with a `df`
  attribute that also counts the estimated kernel parameters and a
  `nobs` attribute (see "Log-likelihood and information criteria").

- [`nobs()`](https://rdrr.io/r/stats/nobs.html) for a Bayesian fit: the
  number of observations.

## Methods that delegate to the final model

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`nobs()`](https://rdrr.io/r/stats/nobs.html),
[`formula()`](https://rdrr.io/r/stats/formula.html),
[`family()`](https://rdrr.io/r/stats/family.html),
[`deviance()`](https://rdrr.io/r/stats/deviance.html), and
[`df.residual()`](https://rdrr.io/r/stats/df.residual.html) call the
same generic on the final model, so their results, and their errors, are
the native ones. For example, `unmarked` fits support
[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html), and
[`residuals()`](https://rdrr.io/r/stats/residuals.html) through the
`unmarked` methods, while
[`family()`](https://rdrr.io/r/stats/family.html) stops with the native
error. Coefficient standard errors from `summary(fit$opt_mod)` or
`vcov(fit$opt_mod)` condition on the estimated scales; see
[`vcov.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/vcov.multiScaleR.md)
for the joint covariance of the analytic backends.

## Log-likelihood and information criteria

[`logLik()`](https://rdrr.io/r/stats/logLik.html) returns the maximized
log-likelihood of the final model. Its `df` attribute counts the native
model parameters plus each estimated scale (sigma) and each estimated
kernel shape, exactly as `K` in
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md).
`AIC(fit)` therefore equals the AIC reported by
`aic_tab(list(fit), AICc = FALSE)`, and `BIC(fit)` equals the BIC from
[`bic_tab()`](https://wpeterman.github.io/multiScaleR/reference/bic_tab.md).
Both are larger than `AIC(fit$opt_mod)` and `BIC(fit$opt_mod)`, which
leave the kernel parameters out. The `nobs` attribute is the sample size
those tables use, for example the event count for Cox and conditional
logistic models.

## Bayesian fits

[`coef()`](https://rdrr.io/r/stats/coef.html),
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html),
[`vcov()`](https://rdrr.io/r/stats/vcov.html), and
[`logLik()`](https://rdrr.io/r/stats/logLik.html) stop for a
`multiScaleR_bayes` object, because a posterior has no single set of
estimates or maximized likelihood. Use
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
for posterior summaries of the coefficients and scales,
[`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md)
for the scale draws, and
[`predict()`](https://rspatial.github.io/terra/reference/predict.html)
for fitted values with their posterior uncertainty.
[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the number of
observations in the fitted data.

## See also

[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`summary.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/summary.multiScaleR.md),
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md),
[`bic_tab()`](https://wpeterman.github.io/multiScaleR/reference/bic_tab.md),
[`vcov.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/vcov.multiScaleR.md),
[multiScaleR](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-class.md).

## Examples

``` r
data(pts)
data(count_data)
hab <- terra::rast(system.file("extdata", "hab.tif",
                               package = "multiScaleR"))
inputs <- kernel_prep(pts = pts, raster_stack = hab, max_D = 250,
                      kernel = "gaussian", verbose = FALSE)
fit <- multiScale_optim(glm(y ~ hab, family = poisson, data = count_data),
                        inputs, verbose = FALSE)
coef(fit)              # coefficients of the model at the estimated scale
#> (Intercept)         hab 
#>   0.4806770   0.8691288 
head(fitted(fit))      # fitted values of that model
#>         1         2         3         4         5         6 
#> 0.8489975 6.3488214 2.1604611 0.9062028 1.0315259 4.4971443 
logLik(fit)            # df counts the estimated sigma
#> 'log Lik.' -114.1443 (df=3)
AIC(fit)               # equals aic_tab(list(fit), AICc = FALSE)$AIC
#> [1] 234.2887
AIC(fit$opt_mod)       # native AIC, without the sigma parameter
#> [1] 232.2887
```
