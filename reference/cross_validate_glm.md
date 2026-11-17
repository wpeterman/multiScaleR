# Cross-Validate an Ordinary GLM Reference Model

Fit a GLM to the analysis rows of each fold and score its predictions
for the held-out rows. No scales of effect are estimated. Use this
function to build null and nonspatial reference models on the same folds
you give
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
so that a scale-of-effect model has a baseline to beat.

## Usage

``` r
cross_validate_glm(
  formula,
  data,
  folds,
  family = c("gaussian", "poisson", "binomial", "nbinom2", "Gamma"),
  weights = NULL,
  glm_control = stats::glm.control(),
  metrics = NULL,
  keep_fits = FALSE,
  keep_selection = FALSE,
  theta = NULL,
  dispersion = NULL
)
```

## Arguments

- formula:

  Two-sided model formula, written as for
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html). The response must
  be a column of `data`, optionally wrapped in
  [`I()`](https://rdrr.io/r/base/AsIs.html), or
  `cbind(successes, failures)` for grouped binomial data. Predictors can
  include factors, interactions, offsets, and the data-dependent terms
  [`scale()`](https://rspatial.github.io/terra/reference/scale.html),
  [`poly()`](https://rdrr.io/r/stats/poly.html),
  [`splines::ns()`](https://rdrr.io/r/splines/ns.html), and
  [`splines::bs()`](https://rdrr.io/r/splines/bs.html). Use `y ~ 1` for
  a null model. Every variable must be a column of `data`.

- data:

  Data frame with one row per observation, holding the raw responses and
  every predictor. Its row names are the row IDs that the folds refer
  to, and they must be unique. Do not supply predictors standardized or
  transformed with all rows (see "What each fold fits").

- folds:

  Named list of explicit splits, for example from
  [`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md)
  or
  [`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md).
  Each split is a list of two character vectors of row IDs, `analysis`
  and `assessment`, which must not overlap. A row can be left out of a
  split entirely (for example, by a spatial buffer) and can be assessed
  in more than one split. Keep rows that share a site, animal, or other
  sampling unit on the same side of every split.

- family:

  Character response family. Choose `"gaussian"` (identity link, the
  default), `"poisson"` (log), `"binomial"` (logit), `"nbinom2"`
  (negative binomial, log), or `"Gamma"` (log). Lowercase `"gamma"` is
  accepted. For `"nbinom2"` the variance is mu + mu^2 / theta, where mu
  is the expected count and theta the size parameter, so a small theta
  means strong extra variation. For `"Gamma"` the variance is dispersion
  x mu^2. Other links and quasi-likelihood families are not supported.

- weights:

  Optional case frequencies, one nonnegative whole number per row of
  `data` in the same order. If the vector has names, they must equal the
  row names. A weight of 3 counts the row three times. A weight of 0
  keeps the row in the records but out of the fit and the scores. The
  default `NULL` gives every row a weight of one. Gaussian models accept
  only weights of one. For binomial data, put the trial counts in a
  `cbind(successes, failures)` response, not in `weights`. Weighted
  proportions are rejected. Some scale-fitting backends add their own
  weight requirements.

- glm_control:

  Named list of
  [`stats::glm.control()`](https://rdrr.io/r/stats/glm.control.html)
  settings for the GLM fitted in each fold. The default keeps the
  settings of [`glm()`](https://rdrr.io/r/stats/glm.html).

- metrics:

  Optional character vector of scores to record. The defaults depend on
  `family`. They are `c("rmse", "mae")` for Gaussian,
  `c("rmse", "mae", "neg_log_score", "deviance")` for Poisson,
  `c("neg_log_score", "brier")` for binomial, and
  `c("rmse", "mae", "neg_log_score")` for NB2 and Gamma. You can request
  a subset of the defaults, but not other scores. NB2 and Gamma leave
  out the deviance because its value depends on theta or the dispersion,
  which makes it unsuitable for comparing models. Lower is better for
  every score, and the
  [multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
  page defines each one.

- keep_fits:

  Logical. The default `FALSE` keeps the predictions, scores, fold
  estimates, and compact diagnostics. `TRUE` also keeps the full fitted
  model of every fold in `fits`, which can use much more memory. You
  need full fits only to inspect the fold models or to have
  [`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md)
  return the selected fits.

- keep_selection:

  Logical. `TRUE` keeps a compact record of each fold's training
  likelihood and parameter count in `selection_records`, which is what
  [`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md)
  needs when `keep_fits = FALSE`. The records come from the fits already
  run, so they add no fitting time. The default is `FALSE`.

- theta:

  Optional NB2 size parameter fixed in advance, one positive number. The
  default `NULL` estimates theta in each fold with
  [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html). A fixed
  theta is fitted with `glm(family = MASS::negative.binomial(theta))`
  and adds no parameter to K. Fix it only from information independent
  of the assessment outcomes. Other families reject it.

- dispersion:

  Optional Gamma dispersion fixed in advance, one positive number. The
  default `NULL` estimates it in each fold as the training deviance
  divided by the sum of the frequency weights, the value native Gamma
  AIC uses. A fixed dispersion adds no parameter to K and is built into
  the Gamma density before fitting. Fix it only from information
  independent of the assessment outcomes. Other families reject it.

## Value

A list of class `multiScaleR_cv`, with class
`multiScaleR_distribution_cv` added for `"nbinom2"` and `"Gamma"`. It
has the elements of a
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md)
result, with these differences. `settings$backend` is `"glm"`.
`training` holds only the coefficients of each fold, plus the frozen
nuisance parameter for NB2 and Gamma. `diagnostics` holds the
convergence record of each fold's GLM, and `fits` holds ordinary `glm`
objects when `keep_fits = TRUE`. See
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
for every element and how to read each score.

## What each fold fits

Each fold fits [`glm()`](https://rdrr.io/r/stats/glm.html) to its
analysis rows with the declared family and frequency weights, then
predicts the expected response of the assessment rows. The families,
response forms, weights, and scores are those of
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md).
Factors, interactions, and offsets use the native GLM machinery. The
constants of
[`scale()`](https://rspatial.github.io/terra/reference/scale.html),
[`poly()`](https://rdrr.io/r/stats/poly.html),
[`splines::ns()`](https://rdrr.io/r/splines/ns.html), and
[`splines::bs()`](https://rdrr.io/r/splines/bs.html) come from the
analysis rows, even inside expressions and offsets. Any other
preprocessing that learns from the data must also happen within the
analysis rows of each fold, so do not supply predictors standardized
with all rows. The section "Formula terms that learn from the data" in
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md)
lists the supported terms and their limits. Missing values are never
dropped silently, and native warnings and failed folds are kept.

## Binomial responses

A binomial response can be numeric or logical 0/1 values, a two-level
factor (the second level counts as a success), or
`cbind(successes, failures)` for several trials per row. Keep trial
counts in the response and case frequencies in `weights`, because the
two are scored differently.

- The negative log score of a row is the negative log binomial
  probability of its observed number of successes, given its trials.
  Rows are averaged with their frequency weights.

- The Brier score of a row is the squared error of the predicted success
  probability, averaged over the row's trials. Rows are averaged with
  weights equal to frequency weight x trials, so the pooled Brier score
  is a mean per trial. A row with many trials counts more in the Brier
  score than in the log score.

- A row with zero trials or zero weight stays in `predictions` but
  contributes to neither score.

- An observed outcome given zero probability produces an infinite log
  score. It is kept as is, not clipped.

## Negative binomial and Gamma responses

`"nbinom2"` fits
[`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) when theta
is estimated and `glm(family = MASS::negative.binomial(theta))` when it
is fixed, so it needs the MASS package. `"Gamma"` uses a log link. In
each fold the nuisance parameter (theta, or the Gamma dispersion) is
estimated from the analysis rows, frozen, and used to score the
assessment rows. Coefficients and their native standard errors are those
of the native fit. A fixed Gamma `dispersion` is built into a normalized
Gamma family before fitting. The section "Negative binomial and Gamma
responses" in
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md)
describes how these parameters enter K and the rules they follow.

### How to read the Gamma dispersion

`distribution[[fold]]$dispersion` is the training deviance divided by
the training sum of case frequencies, the value native Gamma AIC uses.
It is held fixed when the held-out Gamma densities are scored. Native
coefficient standard errors use the Pearson dispersion instead, so they
can rest on a different number. The density dispersion does not replace
that coefficient uncertainty. `convention` records the rule that
produced the value, and `estimated` is `FALSE` when you fixed it with
`dispersion`. Supply `theta` or `dispersion` only when the value is
known or was fixed independently of the assessment outcomes.

## What the scores do and do not include

Predictions use each fold's fitted parameters as if they were known.
They do not integrate parameter uncertainty, provide prediction
intervals, or establish calibration. Use
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md)
to estimate scales of effect within the same folds.

## See also

[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md),
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md),
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md),
and
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md).

## Examples

``` r
data("landscape_counts")   # counts and a site-level covariate
data("surv_pts")           # survey locations
# Three spatial folds. The simulated coordinates are in meters.
folds <- spatial_cv_folds(crds(vect(surv_pts)), block_size = 1000, k = 3,
                          crs = 32617)

# Null and site-only reference models on the same folds
null_cv <- cross_validate_glm(counts ~ 1, landscape_counts, folds,
                              family = "poisson")
site_cv <- cross_validate_glm(counts ~ site, landscape_counts, folds,
                              family = "poisson")
site_cv   # pooled held-out scores, lower is better
#> Cross-validation: 3 of 3 folds scored; 100 of 100 input rows assessed.
#>         metric estimate        unit n_assessment n_scored weight_sum complete
#>           rmse 2.232857 observation          100      100        100     TRUE
#>            mae 1.661678 observation          100      100        100     TRUE
#>  neg_log_score 2.169778 observation          100      100        100     TRUE
#>       deviance 2.372818 observation          100      100        100     TRUE
cv_tab(list(null = null_cv, site = site_cv), metric = "neg_log_score")
#>   Modnames Rank    Score       Delta        Metric        Unit N_assessment
#> 1     null    1 2.163750 0.000000000 neg_log_score observation          100
#> 2     site    2 2.169778 0.006027989 neg_log_score observation          100
#>   N_scored Weight_sum Warnings
#> 1      100        100        0
#> 2      100        100        0

# Negative binomial: theta is estimated in each fold and frozen for scoring
if (requireNamespace("MASS", quietly = TRUE)) {
  nb_cv <- cross_validate_glm(counts ~ site, landscape_counts, folds,
                              family = "nbinom2")
  sapply(nb_cv$distribution, function(d) d$theta)   # small theta: overdispersed
}
#>    fold1    fold2    fold3 
#> 1.561918 1.438781 1.939947 

# Gamma model for a strictly positive response
i <- 1:96
positive <- data.frame(
  y = exp(0.2 + 0.35 * sin(i * 0.73)) *
    rep(c(0.3, 0.6, 1.1, 2.4, 0.8, 1.8, 0.4, 1.6), 12),
  x = sin(i * 0.73))
gamma_folds <- grouped_cv_folds(rep(letters[1:12], each = 8), k = 3)
gamma_cv <- cross_validate_glm(y ~ x, positive, gamma_folds, family = "Gamma")
gamma_cv$summary
#>          metric  estimate        unit n_assessment n_scored weight_sum complete
#> 1          rmse 0.9077540 observation           96       96         96     TRUE
#> 2           mae 0.7704337 observation           96       96         96     TRUE
#> 3 neg_log_score 1.1450485 observation           96       96         96     TRUE
# Training deviance / sum of case frequencies, frozen for scoring each fold
sapply(gamma_cv$distribution, function(d) d$dispersion)
#>     fold1     fold2     fold3 
#> 0.4401327 0.4399048 0.4391975 
```
