# Methods for joint scale fits

Print, summarize, and extract results from a joint fit made with
[`multiScale_joint()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_joint.md).
[`logLik()`](https://rdrr.io/r/stats/logLik.html) makes
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) and
[`stats::BIC()`](https://rdrr.io/r/stats/AIC.html) work, so several
structures fitted to the same data sets can be compared directly.
[`anova()`](https://rdrr.io/r/stats/anova.html) gives likelihood-ratio
tests between nested structures. A test of a `"common"` fit against a
`"footprint_scaled"` fit tests `kappa = 0`, on the boundary of the
parameter space, and uses the matching mixture of chi-square
distributions (half the chi-square p-value for one `kappa`).

## Usage

``` r
# S3 method for class 'multiScaleR_joint'
print(x, digits = 4, ...)

# S3 method for class 'multiScaleR_joint'
summary(object, ...)

# S3 method for class 'summary_multiScaleR_joint'
print(x, digits = 4, ...)

# S3 method for class 'multiScaleR_joint'
logLik(object, ...)

# S3 method for class 'multiScaleR_joint'
nobs(object, ...)

# S3 method for class 'multiScaleR_joint'
coef(object, ...)

# S3 method for class 'multiScaleR_joint'
vcov(object, ...)

# S3 method for class 'multiScaleR_joint'
anova(object, ...)
```

## Arguments

- x, object:

  A `multiScaleR_joint` object.

- digits:

  Number of significant digits to print.

- ...:

  Further `multiScaleR_joint` objects for
  [`anova()`](https://rdrr.io/r/stats/anova.html); ignored by the other
  methods.

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
returns an object of class `summary_multiScaleR_joint` holding the
parameter, scale, comparison, and per-data-set coefficient tables.
[`logLik()`](https://rdrr.io/r/stats/logLik.html) returns a `logLik`
object whose `df` counts every native parameter plus the joint scale
parameters, `kappa`, and shared coefficients.
[`coef()`](https://rdrr.io/r/stats/coef.html) returns the joint
parameters on their natural scales (scales in map units).
[`vcov()`](https://rdrr.io/r/stats/vcov.html) returns their covariance
on the optimizer scale (log scales, `kappa`, coefficients).
[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the total number
of observations. [`anova()`](https://rdrr.io/r/stats/anova.html) returns
a data frame of likelihood-ratio tests.

## See also

[`multiScale_joint()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_joint.md)
