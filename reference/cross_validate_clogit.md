# Cross-Validate Conditional Choices in Complete Matched Sets

Cross-validate a step-selection model or another conditional logistic
model. In each fold, refit the coefficients, and optionally the scales
of effect, using the analysis rows, then score the probability the model
gives the used alternative among all alternatives of each held-out
choice set. A *choice set* is one used alternative together with the
available alternatives it was chosen from, for example one observed step
and its random steps.

## Usage

``` r
cross_validate_clogit(
  formula,
  data,
  folds,
  kernel_inputs = NULL,
  backend = c("conditional_logit", "refit"),
  method = c("exact", "efron", "breslow"),
  cox_control = survival::coxph.control(),
  optim_args = list(),
  keep_fits = FALSE,
  keep_selection = FALSE
)

# S3 method for class 'multiScaleR_choice_cv'
print(x, ...)
```

## Arguments

- formula:

  Two-sided formula. The response is a numeric or logical column coded 1
  for the used alternative and 0 for the available ones. The right side
  holds the predictors and one bare `strata(choice_id)` term naming the
  column that identifies each choice set. Factors, interactions,
  offsets, and the data-dependent terms
  [`scale()`](https://rspatial.github.io/terra/reference/scale.html),
  [`poly()`](https://rdrr.io/r/stats/poly.html),
  [`splines::ns()`](https://rdrr.io/r/splines/ns.html), and
  [`splines::bs()`](https://rdrr.io/r/splines/bs.html) are allowed. A
  dot (`.`) expands to the columns of `data`, so remove the set ID from
  it, for example `case_ ~ . - step + strata(step)`.

- data:

  Data frame with one row per alternative, holding the case column, the
  set ID, and every predictor. Its row names are the row IDs used by
  `folds` and, when supplied, by `kernel_inputs`. They must be unique. A
  tibble, such as the output of
  [`amt::random_steps()`](https://rdrr.io/pkg/amt/man/random_steps.html),
  is used as the equivalent data frame, whose row IDs are `"1"` to the
  number of rows.

- folds:

  Named list of explicit splits of row IDs, as in
  [`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md).
  Every choice set must lie wholly in the analysis rows, wholly in the
  assessment rows, or wholly outside the split.
  [`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md)
  applied to the set IDs, or to animal IDs, does this. See "Build folds
  that keep choice sets whole".

- kernel_inputs:

  Optional `multiScaleR_data` object from
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  or
  [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md)
  with the cell-level data kept, one prepared row per alternative. The
  default `NULL` fits an ordinary conditional logistic reference model
  without scales. When it is supplied, each fold also estimates the
  scales, and learns the standardization and distance bins, from its
  analysis rows only.

- backend:

  Scale-fitting route used with `kernel_inputs`, either
  `"conditional_logit"` (default) or `"refit"`. `"conditional_logit"`
  accepts Gaussian or exponential kernels with kernel means or weighted
  surface SQ, and other kernels or hard-radius predictors stop before
  any fold is fitted. Supply `optim_args = list(radius_grid = ...)` with
  `"refit"` to profile hard radii, together with any smooth scales, in
  each training split. The usual kernel, metric, and numerical checks of
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  apply in every fold. Ignored when `kernel_inputs` is `NULL`.

- method:

  Partial-likelihood method passed to
  [`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html),
  one of `"exact"` (default), `"efron"`, or `"breslow"`. These are rules
  for tied event times. With one used alternative per set they give the
  same likelihood, but the method is recorded, and
  [`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
  and
  [`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md)
  compare only candidates that declare the same method.

- cox_control:

  Named list of
  [`survival::coxph.control()`](https://rdrr.io/pkg/survival/man/coxph.control.html)
  settings for the native training fit. The default keeps the settings
  of `coxph.control()`.

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

- x:

  A `multiScaleR_choice_cv` result, for
  [`print()`](https://rdrr.io/r/base/print.html).

- ...:

  Not used. Reserved for future print controls.

## Value

A list of class `multiScaleR_choice_cv`, which inherits from
`multiScaleR_cv`. Its `summary` holds the negative log score of the used
alternative, averaged over held-out choice sets, and `predictions` holds
one row per alternative with its probability within its set. The other
elements follow
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md).
The conditional-choice section of
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
lists every difference. [`print()`](https://rdrr.io/r/base/print.html)
shows how many folds were scored and the pooled score, and returns `x`
invisibly.

## Build folds that keep choice sets whole

Before any fitting, the function rejects a split that divides a choice
set. Holding out whole choice sets measures transfer to new sets from
the same animals. To measure transfer to new animals, build the folds
from animal IDs instead. If step labels repeat across animals, first
create a globally unique set ID, such as
`interaction(animal, step, drop = TRUE)`, and use that column in
`strata()`. For spatial blocks, pass the set IDs or the animal IDs as
`group` to
[`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md),
so that no set crosses a split boundary.

## What the model can contain

Each set needs exactly one used alternative and at least one available
alternative. Weights, missing values, robust clustering (`cluster()`),
frailty and other special or penalized Cox terms, and general censored
survival responses are outside this function, and no rows are dropped
silently. Predictors that are constant within every set have no
identifiable coefficient. Native alias and convergence problems are kept
as failures rather than resolved by dropping terms or alternatives, and
a warning from the native training fit makes the fold fail. Every scored
fold checks that its native fit is unweighted and unpenalized, uses the
analysis rows and complete sets, and used the declared `method`, whether
or not fits or records are kept. Learned transformations follow the
rules described in
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md).
Their training constants add no stratum effects and do not change the
event count. A failed fold is recorded in `errors`, and when any fold
fails the function gives one warning with the number of failed folds and
their first distinct reasons.

## Estimate scales within the folds

With `kernel_inputs`, each fold rebuilds the spatial covariates from the
raw cells of its analysis rows, learns the standardization and the
representative bin distances there, and estimates the scales. The
assessment rows reuse those frozen references. Inputs from
[`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md)
are fixed covariate maps and add no random effects. Raw covariate units
(`covariate_transform = "raw"` in `optim_args`) require
`backend = "refit"`.
[`vignette("v04-landscape-metrics", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v04-landscape-metrics.md)
shows a search over a grid of radii on a full data set.

## Read the scores

The score of a held-out set is the negative natural log of the
probability the fitted model gives its used alternative, where the
probabilities of a set's alternatives sum to 1. The pooled score
averages over sets. A model without covariates or offsets gives each of
k alternatives a probability of 1/k, so its score is log(k), or 1.099
for three alternatives. A useful covariate lowers the score. Scores are
computed on the log scale, so a probability that prints as zero still
has a finite score. Predictions use the saved training terms and assign
no effect to the labels of new sets. These within-set probabilities
differ from Bernoulli probabilities of use and from unnormalized
relative selection strengths.

For a model with no estimated coefficients, survival's
[`logLik()`](https://rdrr.io/r/stats/logLik.html) can return `NA`. The
function then checks survival's stored null log-likelihood against the
exact complete-set likelihood before using it, with K = 0, for
training-criterion selection. Missing or inconsistent values remain
failures. Other native fits are not altered.

## Compare and select models

Use
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
with `metric = "neg_log_score"` to rank choice models fitted to the same
sets, folds, and `method`. Use
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md)
to score a training AIC, AICc, or BIC rule. Its sample size is the
number of choice sets (one event each), and K counts the coefficients
and any estimated scales and shapes. Choice scores cannot be mixed with
scores from Bernoulli GLMs. They treat fitted parameters as known, so
they do not integrate parameter or selection uncertainty, give
intervals, or establish calibration.
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md)
includes a step-selection example.

## See also

[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
[`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md),
[`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md),
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md),
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md),
and
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md).

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  # 24 choice sets of three alternatives. The first of each set was used.
  i <- 1:72
  choices <- data.frame(step = (i - 1) %/% 3 + 1,
                        used = as.numeric(i %% 3 == 1),
                        habitat = sin(i * 1.7))
  folds <- grouped_cv_folds(choices$step, k = 3)   # whole sets held out

  models <- list(
    null = cross_validate_clogit(used ~ strata(step), choices, folds,
                                 keep_selection = TRUE),
    habitat = cross_validate_clogit(used ~ habitat + strata(step), choices,
                                    folds, keep_selection = TRUE))
  print(models$habitat)   # natural-log loss per held-out choice set

  # The null model scores log(3) = 1.099. Here habitat is unrelated to use,
  # so it does not lower the score.
  print(cv_tab(models, metric = "neg_log_score"))

  # Model picked by training AICc in each fold
  select_cv(models, criterion = "AICc")$selection$chosen
}
#> Conditional-choice CV: 3 of 3 folds scored.
#>         metric estimate       unit n_assessment n_scored weight_sum complete
#>  neg_log_score 1.126067 choice_set           24       24         24     TRUE
#>  n_alternatives
#>              72
#>   Modnames Rank    Score      Delta        Metric       Unit N_assessment
#> 1     null    1 1.098612 0.00000000 neg_log_score choice_set           24
#> 2  habitat    2 1.126067 0.02745496 neg_log_score choice_set           24
#>   N_scored Weight_sum Warnings
#> 1       24         24        0
#> 2       24         24        0
#>  fold1  fold2  fold3 
#> "null" "null" "null" 
```
