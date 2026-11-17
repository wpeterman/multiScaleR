# Summarize Posterior Landscape Scales

Report the scale of each landscape predictor in projected map units, for
every distinct combination of scale-formula covariates and groups, with
its posterior uncertainty and the posterior mass outside the range of
the prepared inputs.

## Usage

``` r
posterior_scales(
  object,
  level = 0.95,
  draws = NULL,
  scale = c("distance", "log"),
  estimator = c("draws", "weighted")
)
```

## Arguments

- object:

  A fitted
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  object with posterior draws, from sampling or either approximation. A
  prepared object has no draws and gives an error.

- level:

  Probability of the central posterior interval; default 0.95.

- draws:

  Optional positive integer number of stored draws to use, taken at
  evenly spaced positions. The default `NULL` uses every stored draw.
  Fewer draws run faster and add simulation noise.

- scale:

  `"distance"` (default) reports scales in projected map units. `"log"`
  reports natural log distances in those units. Log summaries avoid
  exponentiating very large distances, and their means and SDs stay
  finite under the default group-SD prior (see "Read means of
  long-tailed distances").

- estimator:

  `"draws"` (default) summarizes the stored draws, for every method.
  `"weighted"` summarizes all retained production candidates of an
  importance-corrected fit with their normalized importance weights,
  which avoids the extra noise of resampling. It requires `draws = NULL`
  and a fit made with `control = list(correction = "importance")`.

## Value

A data frame with one row per landscape predictor and distinct
combination of scale-formula covariates and groups.

- `covariate`: the landscape predictor.

- `group`: a label for the combination, such as `"year=2021"`, or
  `"shared"` for a predictor with a single scale.

- `row`: a representative observation row with that combination.

- `mean`, `median`, `sd`, `lower`, `upper`: the posterior mean, median,
  and SD, and the limits of the central `level` credible interval.
  Values are in map units, or in log map units when `scale = "log"`. For
  log summaries, [`exp()`](https://rdrr.io/r/base/Log.html) of the
  median and the limits gives distances, but
  [`exp()`](https://rdrr.io/r/base/Log.html) of the mean gives a
  geometric mean, not the arithmetic mean distance.

- `p_below_min_D`, `p_above_max_D`: the posterior probability that the
  scale lies below `min_D` or above `max_D` of the prepared inputs.
  These flag truncation and resolution limits. They are not evidence of
  a biological boundary.

`mean` and `sd` summarize a finite set of draws. A finite value does not
show that the corresponding posterior moment exists. With
`estimator = "weighted"`, quantiles are empirical inverse-CDF quantiles
of the weighted candidates and `sd` is the SD of that finite weighted
population.

## Read annual scales after partial pooling

For a scale formula such as `~1 + (1 | year)`, each row describes one
fitted year's scale after partial pooling. Its median is not a
separately estimated annual value. The shared log-scale intercept
describes a population median, which differs from the arithmetic mean of
the annual scales. Intervals are conditional on the likelihood, the
priors, the maps, and the binning. With `scale = "log"`, `sd` measures
uncertainty in one year's log distance. It differs from the between-year
SD reported in `summary(object)$group_sd`.

## Read means of long-tailed distances

A random scale effect is \\\tau z\\, with standard-normal z and a
half-normal tau whose prior scale is s (`scale_group_sd`). Its prior
exponential moment of order r exists only when \\r s \< 1\\. The default
`scale_group_sd = 0.5` therefore places the variance of annual distances
at the boundary where it becomes infinite, and `scale_group_sd >= 1`
also allows an infinite arithmetic mean distance. The binned likelihood
can keep these tails, because very large distances approach uniform
weighting within the extracted neighborhood. Report medians and interval
limits for distances, and use `scale = "log"` when you need means and
SDs.

## Summarize importance-corrected fits directly

For a fit made with `correction = "importance"`,
`estimator = "weighted"` summarizes the retained candidate population
with its normalized weights instead of the resampled draws. Both target
the same posterior. The weighted version avoids multinomial resampling
noise, but it neither raises the importance effective sample size nor
recovers a posterior region the candidates missed. Candidates with zero
weight are dropped. The cautions about infinite distance moments still
apply.

[`plot.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/plot.multiScaleR_bayes.md)
uses the weighted estimator by default for corrected fits, while this
function and
[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
default to the draws. Request `estimator = "weighted"` in both to keep a
table and a plot consistent.

Weighted summaries are only as precise as the candidate population.
`importance_max` in the `control` list of
[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
is a cap, and a population that passes its checks at a smaller size
stops there. To require more candidates for precise quantiles, raise
`importance_min` as well, as in
`control = list(correction = "importance", importance_min = 65536, importance_max = 131072)`.
Larger populations cost more and do not guarantee a chosen quantile
error. Check the retained count in `fit$diagnostics$candidates`, and
compare independently seeded fits.

## See also

[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
[multiScaleR_bayes](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md),
[`plot.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/plot.multiScaleR_bayes.md),
[`predict_scales()`](https://wpeterman.github.io/multiScaleR/reference/predict_scales.md),
[`diagnose_bayes_scales()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_scales.md),
[`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md);
[`vignette("v13-bayesian-scales", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v13-bayesian-scales.md)
and
[`vignette("v15-bayesian-computation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v15-bayesian-computation.md).

## Examples

``` r
example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
  package = "multiScaleR"))
fit <- example$fits$pooled  # compact draws from an annual-scale NUTS fit
annual <- posterior_scales(fit)
annual[c("group", "median", "lower", "upper", "p_below_min_D")]
#>       group   median    lower    upper p_below_min_D
#> 1 year=2021 173.0828 39.08530 431.8495        0.1935
#> 2 year=2022 173.7422 41.71146 423.5986        0.1840
#> 3 year=2023 143.6104 34.01918 376.3533        0.2665
#> 4 year=2024 175.0161 41.81391 441.1051        0.1665
#> 5 year=2025 155.9366 49.68400 373.9298        0.1870
# Sigma for each year in meters, with central 95% credible intervals.
# Part of each posterior lies below min_D (100 m, one raster cell).

posterior_scales(fit, scale = "log")[c("group", "mean", "sd")]
#>       group     mean        sd
#> 1 year=2021 5.091107 0.6186922
#> 2 year=2022 5.104301 0.6013087
#> 3 year=2023 4.928907 0.5995917
#> 4 year=2024 5.122512 0.5886036
#> 5 year=2025 5.016202 0.5348244
# The log-distance SD measures uncertainty in one year's scale, not the
# variation among years.

# Convert the median sigma to the radius holding 90% of the kernel weight
kernel_dist(sigma = annual$median[1], kernel = "gaussian")
#> [1] 371.43
```
