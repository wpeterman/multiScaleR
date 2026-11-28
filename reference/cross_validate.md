# Cross-Validate Estimated Scales with Explicit Folds

Estimate the scales of effect and the GLM coefficients again inside
every fold, using only that fold's training rows, then score the
predictions for the rows it holds out. A *fold* is one split of the data
into *analysis* rows, which fit the model, and *assessment* rows, which
are held out and predicted. Use this function to ask how well a
scale-of-effect model predicts observations it never saw.

## Usage

``` r
cross_validate(
  formula,
  data,
  kernel_inputs,
  folds,
  family = c("gaussian", "poisson", "binomial", "nbinom2", "Gamma"),
  backend = c("profile_glm", "refit", "joint_nb2"),
  weights = NULL,
  optim_args = list(),
  glm_control = stats::glm.control(),
  metrics = NULL,
  keep_fits = FALSE,
  keep_selection = FALSE,
  theta = NULL,
  dispersion = NULL
)

# S3 method for class 'multiScaleR_cv'
print(x, ...)
```

## Arguments

- formula:

  Two-sided model formula, written as for
  [`stats::glm()`](https://rdrr.io/r/stats/glm.html). The response must
  be a column of `data`, optionally wrapped in
  [`I()`](https://rdrr.io/r/base/AsIs.html), or
  `cbind(successes, failures)` for grouped binomial data. Create any
  other transformed response as a column before the call. Predictors can
  include factors, interactions, offsets, and the data-dependent terms
  described in "Formula terms that learn from the data". Every variable
  must be a column of `data` or a covariate prepared in `kernel_inputs`,
  and at least one predictor must be a covariate whose scale is
  optimized.

- data:

  Data frame with one row per observation, holding the raw responses and
  every nonspatial predictor. Its row names are the row IDs. They must
  be unique and must be the same set of IDs as the row names of
  `kernel_inputs$kernel_dat`. The folds refer to rows by these IDs.
  Spatial covariate columns are rebuilt inside every fold, so any values
  `data` holds for them are ignored. Rows with missing values are never
  dropped silently, so a fold that contains one fails.

- kernel_inputs:

  A `multiScaleR_data` object from
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  or
  [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md)
  for the same observations. It must keep the cell-level distances and
  values (`store_cell_data = TRUE`, the default), because every fold
  rebuilds its covariates from the raw cells. Lean objects are rejected.
  Choose `max_D`, the starting scales, and any metric definitions
  without looking at the assessment outcomes.

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

- backend:

  Fitting route for the scales, passed to
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
  The default `"profile_glm"` is the fast analytic route. `"refit"`
  refits the GLM at every trial scale and accepts formula terms the
  default rejects, such as a spline of an optimized covariate without
  fixed knots. `"joint_nb2"` estimates the scales, coefficients, and
  theta together, and applies only to `family = "nbinom2"` with `theta`
  estimated. If you leave `backend` unset, `"Gamma"` and an NB2 model
  with a fixed `theta` use `"refit"`, the only route they support.
  Naming another backend for them is an error. `"profile_glm"` and
  `"joint_nb2"` accept only smooth kernels, optimized kernel means, and
  optimized weighted surface SQ. A formula with an estimated hard radius
  (a landscape metric or an unweighted surface metric) or another
  weighted surface metric stops before any fold is fitted, so choose
  `backend = "refit"` for those predictors. The usual starting-value,
  convergence, and kernel or metric checks of
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  apply in every fold. See
  [multiScaleR-fitting-methods](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-fitting-methods.md).

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

- optim_args:

  Named list of extra
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  arguments used in every fold, such as `par` or `profile_control`.
  Supply `radius_grid` with `backend = "refit"` to profile hard radii,
  together with any smooth scales, in each training split. An estimated
  hard radius usually needs it, because without it folds often fail (see
  "How each fold is fitted"). Supply `covariate_transform = "raw"` with
  `backend = "refit"` to keep optimized covariates in their own units
  (see "Keep optimized covariates in their own units"). Arguments that
  supply data, the starting model, a custom refit function, parallel
  workers, `backend`, or `verbose` are reserved. Every fold runs on one
  worker. The default [`list()`](https://rdrr.io/r/base/list.html) keeps
  the backend's own defaults. Choose starting values without using the
  assessment outcomes.

- glm_control:

  Named list of
  [`stats::glm.control()`](https://rdrr.io/r/stats/glm.control.html)
  settings for the initial GLM fitted to the analysis rows of each fold.
  The default keeps the settings of
  [`glm()`](https://rdrr.io/r/stats/glm.html).

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
  default `NULL` estimates theta separately from the analysis rows of
  each fold. A fixed theta adds no parameter to K (the parameter count
  used by AIC) and requires the `"refit"` backend. Fix it only from
  information independent of the assessment outcomes. Other families
  reject it.

- dispersion:

  Optional Gamma dispersion fixed in advance, one positive number. The
  default `NULL` estimates it in each fold as the training deviance
  divided by the sum of the frequency weights, the value native Gamma
  AIC uses. A fixed dispersion adds no parameter to K and is built into
  the Gamma density before fitting. Fix it only from information
  independent of the assessment outcomes. Other families reject it.

- x:

  A `multiScaleR_cv` result, for
  [`print()`](https://rdrr.io/r/base/print.html).

- ...:

  Not used. Reserved for future print controls.

## Value

`cross_validate()` returns a list of class `multiScaleR_cv`, with class
`multiScaleR_distribution_cv` added for `"nbinom2"` and `"Gamma"`. Read
`complete` first, then `summary` (pooled scores), `metrics` (scores per
fold), `predictions` (one row per held-out observation), `fold_results`,
and `training` (the scales and coefficients estimated in each fold). See
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
for every element and how to read each score.
[`print()`](https://rdrr.io/r/base/print.html) shows how many folds were
scored and the pooled scores, and returns `x` invisibly.

## How each fold is fitted

Suppose you estimated sigma, the kernel scale, once from all of the data
and then cross-validated only the regression. Every held-out row would
already have helped choose the scale that defines its own predictors,
and the scores would flatter the model. `cross_validate()` avoids this
by repeating the whole fit inside each fold, using the analysis rows
alone.

1.  Rebuild each spatial covariate from the raw cells at the starting
    scales in `kernel_inputs`, then center and scale it with the
    analysis-row mean and standard deviation. Each covariate must have a
    positive, finite standard deviation in the analysis rows.

2.  Fit the initial GLM to the analysis rows.

3.  Estimate the scales and coefficients with
    [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
    No fit to the full data is used as a starting point.

4.  Rebuild the assessment covariates with the fitted scales, the
    training standardization, and the training distance bins, then
    predict the expected response and score it.

The fold fails if the initial GLM, the final GLM, or the scale optimizer
does not converge. Folds run one after another, and a failed fold is
recorded and never retried.

An estimated hard radius, such as a landscape metric whose radius is
optimized, changes the likelihood only in steps as the radius crosses
cell distances. A gradient search on such a likelihood often stops with
a nonzero convergence code, and then the fold fails even when the same
model fitted to all of the data looks acceptable. With
`backend = "refit"`, supply candidate radii through
`optim_args = list(radius_grid = ...)` so each fold profiles the hard
radius over a grid instead (see
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)).
The message of a fold that fails this way names the hard-radius
predictors and this remedy.

## Build the folds

Row names tie the pieces together. Before any fitting, the function
checks that every split names known, unique IDs and that its analysis
and assessment rows do not overlap. It cannot check ecological
independence from the IDs you supply. Keep observations that share a
site, animal, year, or other sampling unit on the same side of each
split, and choose the design before you look at any scores.

Rows left out of a split, such as those removed by a buffer, are neither
fitted nor scored in that split. A row assessed in several splits is
scored each time, every occurrence enters the pooled scores, and
`assessment_counts` records how often each row was assessed. A factor
level present in the assessment rows but absent from the analysis rows
makes that fold fail at prediction. No zero effect is substituted.
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md)
compares random, grouped, spatial, and buffered fold designs on the same
model.

## Choose a family and scores

|  |  |  |  |
|----|----|----|----|
| `family` | Link | Default `metrics` | Backend when unset |
| `"gaussian"` | identity | `rmse`, `mae` | `"profile_glm"` |
| `"poisson"` | log | `rmse`, `mae`, `neg_log_score`, `deviance` | `"profile_glm"` |
| `"binomial"` | logit | `neg_log_score`, `brier` | `"profile_glm"` |
| `"nbinom2"` | log | `rmse`, `mae`, `neg_log_score` | `"profile_glm"`, or `"refit"` with a fixed `theta` |
| `"Gamma"` | log | `rmse`, `mae`, `neg_log_score` | `"refit"` |

Lower is better for every score. RMSE and MAE judge only the predicted
mean, in the units of the response. The negative log score judges the
whole predicted distribution, so it also penalizes a model that is too
confident. Binomial responses can be 0/1 numbers, logical values, a
two-level factor, or `cbind(successes, failures)`.
[`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md)
explains how binomial log and Brier scores are weighted, and the
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
page defines every score.

### Negative binomial and Gamma responses

Both families carry a nuisance parameter, a parameter of the response
distribution that is not part of the mean model. It is theta for NB2 and
the dispersion for Gamma. Unless you fix it with `theta` or
`dispersion`, each fold estimates it from the analysis rows and freezes
it to score the assessment rows. The negative log score therefore uses a
complete, normalized density that the held-out outcomes did not help
shape. The frozen value of each fold is stored in `distribution` and
repeated in the `theta` and `dispersion` columns of `predictions`.

- An estimated nuisance parameter adds one to the native K, and a fixed
  one adds zero. The native `nobs` counts rows with positive weight, not
  the sum of the weights.

- Estimated Gamma dispersion is the training deviance divided by the sum
  of the frequency weights, the value native Gamma AIC uses. It differs
  from the Pearson dispersion that
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html)
  uses for coefficient standard errors, and it does not replace that
  standard-error convention. A fixed Gamma dispersion is built into the
  Gamma density before fitting, so the likelihood, AIC, and K agree.
  Setting only the `dispersion` field of a native Gamma family would
  leave its AIC estimating the dispersion from the deviance.

- NB2 responses must be nonnegative whole numbers and Gamma responses
  strictly positive, including rows with zero weight.

- A theta warning from the native NB2 fit, such as a reached iteration
  limit, makes the fold fail. The warning is kept in `warnings`.

- `backend = "joint_nb2"` starts the coefficients from a GLM with theta
  = 1 and then estimates theta jointly with the scales. That starting
  value is neither a fixed theta nor a prior. A fold reports an infinite
  theta, the exact Poisson limit, only when the joint fit verifies and
  flags that limit. A large finite theta never stands in for it.

A Gamma model for a strictly positive response uses the same call, and
the backend switches to `"refit"` on its own, for example
`cross_validate(mass ~ habitat, dat, inputs, folds, family = "Gamma")`.
Compare NB2 and Poisson candidates with the negative log score on the
same counts.
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
lists the comparisons each family allows.

## Formula terms that learn from the data

Some formula terms learn constants from the data.
[`scale()`](https://rspatial.github.io/terra/reference/scale.html)
learns a mean and standard deviation,
[`poly()`](https://rdrr.io/r/stats/poly.html) learns polynomial
coefficients, and [`splines::ns()`](https://rdrr.io/r/splines/ns.html)
and [`splines::bs()`](https://rdrr.io/r/splines/bs.html) learn knots. In
each fold these constants come from the analysis rows and are reused
unchanged for the assessment rows. This holds with or without a
`base::`, `stats::`, or `splines::` prefix, and inside arithmetic,
[`I()`](https://rdrr.io/r/base/AsIs.html), pointwise functions such as
[`log()`](https://rdrr.io/r/base/Log.html), and offsets. For spatial
covariates the constants come from the final fitted training values, not
the starting-scale values. Arguments you fix in the formula, such as
`knots`, stay fixed. Coefficients, standard errors, likelihoods, and K
are those of the native fit, and fits kept with `keep_fits = TRUE` carry
the same training constants.

Only these four functions are recognized, and only with numeric inputs
and supported arguments. An unsupported argument makes the fold fail.
Any other operation that learns from the data, such as `mean(x)` inside
a formula or your own transformation function, is not detected. Compute
it within each analysis split and apply it unchanged to the assessment
rows, or the scores will leak assessment information. For the same
reason, do not supply predictors standardized with all rows. Every
variable must still be a column of `data` or a prepared covariate.
Variables are never taken from the formula's environment.

With the default `"profile_glm"` backend, a spline of an optimized
covariate needs fixed `knots` and `Boundary.knots` and no `df`.
Otherwise the fold fails with a message that names the fix.
`backend = "refit"` accepts adaptive splines.

## Keep optimized covariates in their own units

By default each fold centers and scales every optimized covariate with
its analysis-row mean and standard deviation, so a coefficient is the
effect of a one standard deviation change. With `backend = "refit"`, the
setting `optim_args = list(covariate_transform = "raw")` skips that
final step for the optimized covariates. They keep the units of their
kernel or metric calculation, such as a kernel-weighted proportion of
habitat. Use this when a fixed spline knot or another formula constant
is meant in those units.

The initial fit and every scale refit within a fold use the same units,
and prediction applies a center of zero and a scale of one. Covariates
whose scales are not optimized are still standardized with analysis-row
values, and normalization applied upstream, such as the reference
normalization of surface metrics, is unchanged. Assessment rows never
choose the transformation, the knots, or the scale. Decide on raw units,
and on any knots in those units, before you look at the results. Because
a change of units can alter a constrained or nonlinear model, compare
the two versions by their held-out scores rather than assume they are
equivalent.

## Use separate maps for groups of surveys

Inputs from
[`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md),
which take each row's covariates from its own map (for example, the map
of its survey year), are supported. Before fitting, the function checks
that the named `raster_group` labels align with the prepared row IDs. It
records the map labels of each fold's analysis and assessment rows in
`raster_groups`. Standardization and distance bins still come from the
pooled analysis rows. The map labels are fixed inputs. They add no
random effects, no extra parameters, and no scales that vary among maps.
Group effects you write in the formula, such as a year factor, keep
their native fit, parameter count, and new-level check. See
[`vignette("v06-year-specific-rasters", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v06-year-specific-rasters.md).

## Distance bins and the extraction envelope

When `kernel_inputs` holds distance-binned summaries
([`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
with `bin = TRUE`, the default), each fold rebuilds them. The bin edges
span `max_D` plus half a cell diagonal, so cells that intersect the
buffer but whose centers lie just beyond `max_D` keep their values.
These edges depend on the declared geometry only, not on any data. The
representative distance of each bin is the mean cell distance in that
bin over the analysis rows, and the assessment rows reuse it. Inputs
saved by older versions of
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
lack the cell geometry and keep the stricter rule that every cell center
lies within `max_D`. If such inputs are rejected, prepare them again
with the current
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).

## Check every fold before comparing models

Read `complete` first. `TRUE` means every planned fold was fitted and
scored. When a fold fails, `complete` is `FALSE` and the pooled scores
in `summary` are `NA`. The scores of the folds that finished, the failed
assessment rows, the diagnostics, and the warnings all stay available,
and `errors` names the stage and message of each failure. The function
also gives one warning when it returns an incomplete result, with the
number of failed folds and their first distinct reasons.
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
and
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md)
refuse an incomplete result rather than rank a model on its easier
folds. A complete result says only that the calculations finished, so
inspect `diagnostics` and `warnings` before you rank models.

## What the scores do and do not include

The predictions and scores treat each fold's fitted parameters as known.
They do not average over uncertainty in the coefficients or scales, give
prediction intervals, or establish calibration. Compare candidates only
on identical folds, assessment rows, responses, and weights, which
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
checks for you. The scores complement
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md)
but are not AIC values or Akaike weights. If you use the scores to
choose a model, kernel, or metric, the winning score is optimistic.
Assess the chosen procedure with an outer validation layer or an
independent test set.
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md)
scores a predeclared training-AIC rule without this bias, because the
assessment outcomes never choose the winner within a fold.

`cross_validate()` fits GLMs. For other model structures use
[`cross_validate_mixed()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_mixed.md)
(a random intercept for repeated surveys),
[`cross_validate_clogit()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_clogit.md)
(step selection),
[`cross_validate_field()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_field.md)
(a spatial random field), or
[`cross_validate_bayes()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_bayes.md)
(Bayesian models). GAMs and `unmarked` models have no cross-validation
function.
[`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md)
scores reference models without scales on the same folds.

## See also

[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
for the returned object,
[`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md)
for reference models,
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
and
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md)
to compare candidates,
[`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md)
and
[`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md)
to build folds, and
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md)
for the full workflow.

## Examples

``` r
data("landscape_counts")   # counts and one site-level covariate
data("surv_pts")           # 100 survey locations
pts <- vect(surv_pts)

# Folds are lists of row IDs, so build and check them before any fitting.
# The simulated coordinates are in meters, and crs = 32617 declares that.
folds <- spatial_cv_folds(crds(pts), block_size = 1000, k = 3, crs = 32617)
folds   # rows fitted (analysis) and held out (assessment) in each fold
#> Cross-validation folds: spatial 
#>   fold analysis assessment excluded
#>  fold1       66         34        0
#>  fold2       67         33        0
#>  fold3       67         33        0

# \donttest{
land <- rast(system.file("extdata", "landscape.tif", package = "multiScaleR"))
inputs <- kernel_prep(pts = pts, raster_stack = land[[c("land1", "land2")]],
                      max_D = 1000, kernel = "gaussian", verbose = FALSE)
# Row names of dat match those of inputs$kernel_dat ("1" to "100")
dat <- data.frame(landscape_counts, inputs$kernel_dat)

# Scales and coefficients are estimated from the analysis rows of each fold
cv <- cross_validate(counts ~ site + land1 + land2, data = dat,
                     kernel_inputs = inputs, folds = folds,
                     family = "poisson", keep_selection = TRUE)
cv            # pooled held-out scores, lower is better
#> Cross-validation: 3 of 3 folds scored; 100 of 100 input rows assessed.
#>         metric estimate        unit n_assessment n_scored weight_sum complete
#>           rmse 1.648342 observation          100      100        100     TRUE
#>            mae 1.237010 observation          100      100        100     TRUE
#>  neg_log_score 1.628450 observation          100      100        100     TRUE
#>       deviance 1.290161 observation          100      100        100     TRUE
cv$complete   # TRUE: every fold was fitted and scored
#> [1] TRUE

# Sigma (map units) estimated in each fold. Similar values across folds
# suggest that no single block of surveys drives the estimated scales.
round(sapply(cv$training, function(fold) fold$scale_est$Mean))
#>      fold1 fold2 fold3
#> [1,]   236   246   239
#> [2,]   558   541   566

# Compare with a reference model that has no landscape covariates
site_cv <- cross_validate_glm(counts ~ site, dat, folds, family = "poisson",
                              keep_selection = TRUE)
cv_tab(list(landscape = cv, site_only = site_cv), metric = "rmse")
#>    Modnames Rank    Score     Delta Metric        Unit N_assessment N_scored
#> 1 landscape    1 1.648342 0.0000000   rmse observation          100      100
#> 2 site_only    2 2.232857 0.5845151   rmse observation          100      100
#>   Weight_sum Warnings
#> 1        100        0
#> 2        100        0
# }
```
