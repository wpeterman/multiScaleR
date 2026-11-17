# Predict Landscape Scales for New Covariates and Groups

Carry joint posterior uncertainty into the scale of each landscape
predictor at new values of the scale-formula covariates, including years
that were not surveyed when the scale formula has a random year effect.
This predicts the scale itself, so no prepared landscape is needed.

## Usage

``` r
predict_scales(
  object,
  newdata,
  new_levels = c("error", "sample", "zero"),
  re.form = NULL,
  level = 0.95,
  draws = NULL,
  summary = TRUE,
  scale = c("distance", "log"),
  seed = NULL
)
```

## Arguments

- object:

  A fitted
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  object with stored joint posterior draws, from sampling or either
  approximation.

- newdata:

  A nonempty data frame that holds the fixed scale-formula covariates
  and the grouping variables of any included scale random effects.
  Responses, response-only covariates, and prepared landscape bins are
  not needed.

- new_levels:

  What to do with a scale group level that was not in the fitted data.
  `"error"` (default) stops and asks for a decision. `"sample"` draws a
  new standardized group effect for each joint posterior vector.
  `"zero"` sets the new deviation to zero. Known levels keep their
  fitted effects. Neither choice creates a coefficient for an unseen
  level of a fixed factor.

- re.form:

  `NULL` (default) includes the scale random effects. `NA` sets every
  scale deviation to zero, a reference value rather than an average over
  group variation. Response effects are never used here.

- level:

  Probability of the central interval; default 0.95.

- draws:

  Optional positive integer number of evenly spaced stored posterior
  vectors to use. The default `NULL` uses every vector.

- summary:

  Logical. `TRUE` (default) returns interval summaries and computes each
  repeated scale combination once, which avoids a large matrix when many
  rows share one year's scale. `FALSE` returns one matrix of draws by
  new observations for each landscape predictor.

- scale:

  `"distance"` (default) returns projected map units. `"log"` returns
  natural log distances in those units.

- seed:

  Optional positive integer seed for sampled new effects. Your own
  random-number state is restored afterward, even when prediction fails.

## Value

With `summary = TRUE`, a data frame with one row per landscape predictor
and distinct active combination of scale covariates and groups.

- `covariate`: the landscape predictor.

- `group`: a label for the covariates and groups of the row, such as
  `"year=2026"`.

- `row`: a representative row of `newdata`.

- `new_level`: whether the row contains an unseen scale group.

- `distribution`: which question the row answers. See "Read the
  distribution column".

- `mean`, `median`, `sd`, `lower`, `upper`: summaries of the draws on
  the requested scale, with the central `level` interval. Prefer the
  median and the limits for distances. See
  [`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md)
  for why the arithmetic mean and SD of a distance can be poorly
  defined.

- `p_below_min_D`, `p_above_max_D`: the probability that the scale lies
  below `min_D` or above `max_D` of the fitted inputs. See "Read the
  support columns".

With `summary = FALSE`, a named list with one numeric matrix per
landscape predictor. Each matrix has one row per posterior vector and
one column per row of `newdata`, in the supplied order. Rows in the same
new group share its effect.

Both forms carry a `new_groups` attribute with `component`, `covariate`,
`group`, `level`, and `action`, listing the unseen levels that were
sampled or set to zero. Summaries also carry a `prediction_contract`
attribute (see
[prediction_contract](https://wpeterman.github.io/multiScaleR/reference/prediction_contract.md))
that describes the target, the included scale effects, and the
uncertainty. It does not establish diagnostic eligibility or
calibration.

## Read the distribution column

- `"conditional_posterior"`: a surveyed group, using its estimated
  deviation.

- `"new_group_predictive"`: the range of scales a new group could have.
  Its interval carries uncertainty in the typical scale and in the
  between-group SD, plus the variation among groups. It is a marginal
  prediction interval, not a simultaneous band and not an estimate of a
  new year's scale after its responses are observed.

- `"zero_deviation_reference"`: the new group placed exactly at the
  typical scale (with `new_levels = "zero"`), or every group with
  `re.form = NA`. Its interval keeps uncertainty in the typical scale
  but leaves out variation among groups, so it is narrower without being
  a more precise forecast.

The zero-deviation reference is not the average scale over years. With
\\\log D_y = a + \tau z_y\\, even for fixed a and tau the mean of
\\D_y\\ over years is \\\exp(a + \tau^2/2)\\, while the zero-deviation
scale is \\\exp(a)\\. Averaging over the uncertainty in a and tau also
depends on their posterior dependence.

## Read the support columns

`p_below_min_D` and `p_above_max_D` refer to the range of the fitted
inputs. They say nothing about the coverage of a future map you have not
supplied. The distance priors have no upper or lower bound, so these
probabilities matter when much of the posterior lies outside the
distances your maps resolve. A scale below `min_D` is finer than one
raster cell.

To check a future map, compare the scale draws with that map's `max_D`.
The scale parameter is not the reach of the kernel. Under an ideal,
untruncated planar kernel, the radius that holds 95% of the weight is
`kernel_dist(prob = 0.95, sigma = 1, kernel = fit$specification$kernel)`
times the scale. The examples compute both probabilities. Neither
measures actual cell coverage near raster edges or gaps, which you
should inspect separately. A smaller future buffer renormalizes a more
truncated kernel. If larger neighborhoods are scientifically needed,
rebuild comparable fitting and prediction inputs and refit, rather than
changing the predictor's support only at prediction time.

## Predict new groups jointly

A new year's sampled effect is shared by all of that year's rows within
one posterior vector. Different new years are independent given the
population parameters, whose posterior dependence is kept. Effects for
separate model components follow their specified independence. Separate
calls draw separate new effects, so put rows that belong together in one
batch. Exchangeable annual effects impose no temporal autocorrelation. A
fitted numeric trend extrapolates under that trend assumption, and an
unseen level of a fixed year factor is an error.

## Know the limits

The method uses the stored joint draws, which for an
importance-corrected fit are stored resamples. It adds no posterior
information and does not validate a failed fit. A nonfinite distance
stops with an error, so use `scale = "log"` for extreme scales. Fitted
contrasts and prediction-aware bases such as
[`scale()`](https://rspatial.github.io/terra/reference/scale.html),
[`poly()`](https://rdrr.io/r/stats/poly.html), and the splines keep
their training constants, including inside nested expressions such as
`I(scale(temp))`. Fits saved by older versions need verified recovery
from the original training data. An expression such as `I(x - mean(x))`
takes its mean from the new batch by ordinary R formula rules, so use a
prediction-aware basis or precomputed training constants when the
original centering matters. Formula transforms must come from base R or
stats, or carry an explicit package namespace.

## See also

[`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md),
[prediction_contract](https://wpeterman.github.io/multiScaleR/reference/prediction_contract.md),
[`predict.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR_bayes.md),
[`simulate.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/simulate.multiScaleR_bayes.md),
[`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md),
[multiScaleR_bayes](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md);
[`vignette("v14-bayesian-prediction", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v14-bayesian-prediction.md).

## Examples

``` r
example <- readRDS(system.file("extdata", "bayesian-prediction-example.rds",
  package = "multiScaleR"))
fit <- example$fit   # habitat scale varies among years 2021 to 2025
years <- data.frame(year = c("2025", "2026"))   # 2026 was not surveyed

# A surveyed year, then a new year drawn from the spread among years
predict_scales(fit, years, new_levels = "sample", draws = 500,
  seed = 801)[c("group", "distribution", "median", "lower", "upper")]
#>       group          distribution   median    lower    upper
#> 1 year=2025 conditional_posterior 152.7981 43.31996 375.6791
#> 2 year=2026  new_group_predictive 161.1544 37.54491 513.8587
# The reference: the new year placed at the typical scale
predict_scales(fit, years[2, , drop = FALSE], new_levels = "zero",
  draws = 500)[c("group", "distribution", "median", "lower", "upper")]
#>       group             distribution   median    lower    upper
#> 1 year=2026 zero_deviation_reference 157.4092 49.96857 408.1584

# Can a future map of the same extent hold the new year's kernel?
log_D <- predict_scales(fit, years[2, , drop = FALSE], new_levels = "sample",
  scale = "log", summary = FALSE, draws = 500, seed = 806)$habitat
max_D <- example$kernel_inputs$max_D   # 1200 m for the future inputs
radius_per_D <- kernel_dist(prob = 0.95, sigma = 1,
  kernel = fit$specification$kernel)
mean(log_D > log(max_D))                       # scale itself beyond max_D
#> [1] 0.002
mean(log_D + log(radius_per_D) > log(max_D))   # 95% radius beyond max_D
#> [1] 0.03
# The second probability is larger, because the kernel reaches past its
# scale parameter.
```
