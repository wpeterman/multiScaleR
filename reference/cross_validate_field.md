# Cross-Validate a Spatial Response Field

Cross-validate a model that adds a spatial random field to its fixed
predictors. A *spatial random field* is a smooth surface of spatially
correlated variation that the predictors do not explain. In each fold
the function estimates the field, its correlation range, the regression
coefficients, and optionally the scales of effect from the analysis
rows. It then scores the joint distribution of each held-out block, so
that correlation among nearby held-out observations counts. Gaussian
responses use exact Gaussian conditioning at the fitted parameters.
Count, binary, and positive responses use Laplace integration.

## Usage

``` r
cross_validate_field(
  formula,
  data,
  coords,
  folds,
  audit_dir,
  crs = NULL,
  covariance = c("exponential", "gaussian"),
  kernel_inputs = NULL,
  optim_args = list(),
  control = list(iter.max = 1000L, eval.max = 1500L, rel.tol = 1e-10, x.tol = 1e-08),
  metrics = c("rmse", "mae", "neg_log_score"),
  max_refits = 200L,
  max_memory_mb = 512,
  keep_fits = FALSE,
  family = c("gaussian", "poisson", "binomial", "nbinom2", "Gamma")
)
```

## Arguments

- formula:

  Two-sided fixed-effects formula whose response is a numeric column of
  `data`. Products, interactions, factors, offsets, and the
  data-dependent terms
  [`scale()`](https://rspatial.github.io/terra/reference/scale.html),
  [`poly()`](https://rdrr.io/r/stats/poly.html),
  [`splines::ns()`](https://rdrr.io/r/splines/ns.html), and
  [`splines::bs()`](https://rdrr.io/r/splines/bs.html) keep their usual
  meanings, with their constants learned from the analysis rows. Do not
  add random effects or penalized smooths, because the function adds the
  one shared field itself.

- data:

  Data frame of raw responses and predictors, with unique row names that
  match the fold IDs, in the same row order as `coords`. Supply raw
  predictors, not transformations learned from all rows. Every row has a
  case weight of one.

- coords:

  Locations of the observations, one per row of `data`, as projected
  `sf` points, a
  [`terra::SpatVector`](https://rspatial.github.io/terra/reference/SpatVector-class.html)
  of points, or a numeric two-column table. Rows with identical
  coordinates share one field value, so repeated surveys of a location
  are linked.

- folds:

  Named list of explicit splits, usually from
  [`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md).
  A location cannot appear in both the analysis and the assessment rows
  of one split, so keep repeated surveys of a location together. Rows
  excluded from a split, for example by a buffer, enter neither its fit
  nor its score.

- audit_dir:

  Directory where compact, write-once records of each native fit and its
  starting values are saved, one subfolder per fold. Records of failed
  fits are kept too. No models, tapes, or dense covariance matrices are
  saved. Give each candidate model its own directory to keep the records
  apart.

- crs:

  Projected coordinate reference system, in any form accepted by
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html).
  It is required when `coords` is a coordinate table. The default `NULL`
  uses the CRS of spatial `coords`, which must be projected, and a
  supplied `crs` must match it. Coordinates are never transformed, so
  the correlation range is reported in the CRS units.

- covariance:

  Correlation function of the field, either `"exponential"` (default),
  with correlation exp(-d / range) at distance d, or `"gaussian"`, with
  correlation exp(-(d / range)^2). The range and the field SD are
  estimated in every fold, as is the observation dispersion for
  Gaussian, NB2, and Gamma models.

- kernel_inputs:

  Optional `multiScaleR_data` object from
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  with the cell-level data kept and row names matching `data`. The
  default `NULL` fits a field model without estimating scales. When it
  is supplied, each fold also estimates the scales of the optimized
  covariates that appear in `formula` (list them with
  [`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md))
  by refitting the field model.

- optim_args:

  Named list of extra
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  arguments for the scale fit, used only with `kernel_inputs`. The
  default is an empty list. This function fixes the backend (`"refit"`),
  the data, the refit function, and the single worker, so those cannot
  be replaced.

- control:

  Named list of native [`nlminb()`](https://rdrr.io/r/stats/nlminb.html)
  optimizer settings for each field fit. Accepted names are `iter.max`,
  `eval.max`, `rel.tol`, `x.tol`, and `abs.tol`, each a positive number.
  The defaults are `iter.max = 1000`, `eval.max = 1500`,
  `rel.tol = 1e-10`, and `x.tol = 1e-8`. The settings are stored with
  each fit, and the convergence and curvature checks apply whatever they
  are.

- metrics:

  Scores to record, any of `"rmse"`, `"mae"`, and `"neg_log_score"` (all
  three by default). The negative log score uses one joint density, or
  one joint probability for counts and 0/1 data, per assessment block.
  The pooled score sums the block losses and divides by the number of
  assessment observations.

- max_refits:

  Positive whole number, the most native refits allowed in one fold
  during the scale search (default 200). Each refit checks two starting
  values. A fold that exceeds the budget fails.

- max_memory_mb:

  Positive number, a conservative budget in MiB for the covariance
  calculations of one fold (default 512), checked before any fitting.
  Repeated surveys at one location share one field value, so they do not
  enlarge the training covariance matrix. This is a planning check, not
  an operating-system memory limit.

- keep_fits:

  Logical. `TRUE` keeps the complete training model of each fold in
  `fits`. The default `FALSE` saves memory. The compact training records
  that
  [`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md)
  reads are kept either way.

- family:

  Response family and link. Use `"gaussian"` (identity, the default) for
  continuous observations, `"poisson"` or `"nbinom2"` (log) for counts,
  `"binomial"` (logit) for numeric 0/1 observations, or `"Gamma"` (log)
  for strictly positive continuous observations. All families other than
  `"gaussian"` are experimental.

## Value

A list of class `multiScaleR_field_cv`, which inherits from
`multiScaleR_cv`. `summary` holds the pooled scores, `predictions` the
predicted mean, errors, and variance components of each held-out
observation, `joint_scores` the joint loss of each held-out block, and
`training` the field SD, correlation range (in coordinate units),
dispersion, and coefficients of each fold. Compact `selection_records`
are always kept. The spatial field section of
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
describes every element and its scale.

## Details

The count, binary, and positive response families are experimental.
Their arguments and results may change in a later release without a
deprecation period.

## What each fold estimates

The native fits use glmmTMB and TMB with maximum likelihood. The model
has one shared latent Gaussian field with exponential or Gaussian
correlation, and a constant observation dispersion where the family has
one. From the analysis rows each fold estimates the regression
coefficients, the field SD, the correlation range, the observation
dispersion, and, with `kernel_inputs`, the scales of effect. Distance
bins, standardization, and the constants of learned formula terms also
come from the analysis rows only. Custom preprocessing that depends on
the data must be learned inside each training fold too, because formula
evaluation cannot detect outside information passed in through your own
functions. No missing rows are dropped silently.

The correlation range is the distance, in coordinate units, over which
the unexplained variation stays correlated. It is separate from the
kernel scale of effect, which describes how far the landscape around a
point affects it.

Every fit tries two deterministic starting values chosen from the
training data and records the package versions, convergence, marginal
score, and curvature checks of each. For Gaussian responses both starts
must pass and agree, or the fold fails. For the other families the
passing start with the highest likelihood is used, and the outcome of
each start is kept. With few locations in a training fold, the field's
standard deviation and range are often weakly identified, and the two
Gaussian starts can stop at different points of a nearly flat
likelihood. Those folds fail by design ("starts disagree") rather than
report one arbitrary solution. Use more locations per fold, or fewer
folds, and inspect `audit_dir` to see how far apart the starts ended.

## Build folds for a field model

Hold out complete locations, so that no location has observations in
both parts of a split. Spatial blocks with a buffer test prediction at
new places. The field locations in `coords` can differ from the focal
points used for the landscape covariates, but each location must be
given by its actual, consistent coordinates, and `coords` is required
even when `kernel_inputs` is supplied. A template that compares two
correlation functions on the same folds:

    folds <- spatial_cv_folds(field_xy, block_size = 10000, k = 5,
      buffer = 1000, crs = 32617, ids = rownames(raw_data))
    covariances <- c(exponential = "exponential", gaussian = "gaussian")
    candidates <- lapply(covariances, function(covariance) {
      cross_validate_field(y ~ habitat * cover + offset(log(exposure)),
        raw_data, field_xy, folds, crs = 32617, covariance = covariance,
        kernel_inputs = inputs,
        audit_dir = file.path("field-cv-audits", covariance))
    })
    cv_tab(candidates, metric = "neg_log_score")
    selected <- select_cv(candidates, criterion = "AICc")
    selected$selection$chosen   # chosen by the training criterion in each fold

## How the joint score works

The negative log score treats the held-out observations of each fold as
one block. For Gaussian responses it is the negative log of one
multivariate normal density. Its covariance holds the field variation
that remains after conditioning on the analysis responses, plus the
observation noise. For the other families the score is the Laplace log
marginal likelihood of the analysis responses alone minus that of the
analysis and assessment responses together. Both calculations use
parameters estimated from the analysis rows only, including the field
covariance, coefficients, dispersion, and scales, so the score keeps the
dependence among held-out locations and repeated observations. For
counts and 0/1 data, a block probability above one beyond numerical
roundoff is recorded as a failed score.

The pooled score is the total block loss divided by the number of
assessment observations, and smaller is better. RMSE and MAE summarize
the errors of the predicted response means. Multiplying separate
marginal densities would ignore the field dependence and define a
different score, so
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
keeps these scores apart from the per-observation scores of ordinary
GLMs.

## Predicted means and spreads

For 0/1 responses, predicted means and observation SDs integrate over a
Gaussian approximation to the field given the analysis responses. For
log-link families, the first two response moments come from Laplace
ratios of exponentially tilted analysis likelihoods, with two additional
latent solves per distinct assessment location and the training
covariance factor held fixed. This improves the accuracy of the moments
in low-dimensional numerical comparisons, particularly for sparse
observations with strong fields. The joint block score uses its own
integration.

`observation_predictive_SE` combines the conditional observation noise
with the variation of the response mean. For a nonlinear link, square it
to get the total response variance. The two field variance columns
describe the field on the linear-predictor scale. All of these condition
on the estimated parameters. Use them as predictive spreads for new
observations, not as standard errors of the predicted mean.

## Check calibration separately

A scored fold has an eligible native fit and a valid score. It does not
show that the predictive spreads are reliable. Assess that on the
held-out locations, and count failed folds in the denominator when you
report how often the method succeeds. In simulations with buffered
held-out locations, squared prediction errors averaged 1.5 (Gaussian),
2.7 (Poisson), and 1.9 (NB2) times the predicted observation variance,
against 1 for a calibrated prediction. Treat the returned spreads as
optimistic when parameters are estimated.

## Compare and select field models

Use
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
to rank field results fitted to the same responses, folds, and
coordinates. Poisson and NB2 results share a count target. Use
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md)
to score a training AIC, AICc, or BIC rule from the compact records. The
coefficients, field SD, estimated range, and any observation dispersion
count toward K, while the integrated field values add nothing. Field
scores do not integrate parameter or scale uncertainty and do not
establish calibration.
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md)
explains fold design and how to read the scores, and
[`vignette("v11-spatial-fields", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v11-spatial-fields.md)
covers spatial field models.

## See also

[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
[`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md),
[`response_field_geometry()`](https://wpeterman.github.io/multiScaleR/reference/response_field_geometry.md),
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md),
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md),
and
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md).

## Examples

``` r
# \donttest{
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  data("landscape_counts")
  data("surv_pts")
  xy <- crds(vect(surv_pts))   # simulated coordinates in meters
  folds <- spatial_cv_folds(xy, block_size = 1000, k = 3, crs = 32617)

  # Poisson counts with a site covariate and an exponential field. The
  # Poisson family is experimental.
  field_cv <- cross_validate_field(counts ~ site, landscape_counts, xy,
    folds, audit_dir = file.path(tempdir(), "field-cv"), crs = 32617,
    family = "poisson")
  print(field_cv)   # the log score is a joint score per held-out survey

  # Field SD (log-link scale) and correlation range (m) in each fold
  print(sapply(field_cv$training, function(fold)
    unlist(fold$parameters[c("field_SD", "range")])))
}
#> Cross-validation: 3 of 3 folds scored; 100 of 100 input rows assessed.
#>         metric estimate n_assessment                             unit complete
#>           rmse 1.760629          100                         response     TRUE
#>            mae 1.328896          100                         response     TRUE
#>  neg_log_score 1.703037          100 joint_block_mass_per_observation     TRUE
#>                 fold1      fold2       fold3
#> field_SD    0.9617028   0.812557   0.8259306
#> range    1346.5835180 624.511124 859.2548280
# }
```
