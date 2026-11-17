# Validate Training-Only Model Selection Within CV Folds

Score a model-selection rule rather than a single model. In each fold
the rule picks the candidate with the lowest training AIC, AICc, or BIC,
and the result keeps the held-out predictions that candidate already
made. No models are refitted.

## Usage

``` r
select_cv(
  mod_list,
  criterion = c("AICc", "AIC", "BIC"),
  ties = c("error", "first"),
  keep_fits = FALSE
)
```

## Arguments

- mod_list:

  Named list of complete cross-validation results, all scored on the
  same folds and observations. Each must keep its training records
  (`keep_selection = TRUE`) or full fits (`keep_fits = TRUE`). Results
  from
  [`cross_validate_field()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_field.md)
  and
  [`cross_validate_mixed()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_mixed.md)
  always keep their records.

- criterion:

  Training criterion, declared in advance: `"AICc"` (default), `"AIC"`,
  or `"BIC"`. It is computed from each candidate's training
  log-likelihood and parameter count K with the same rules as
  [`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md)
  and
  [`bic_tab()`](https://wpeterman.github.io/multiScaleR/reference/bic_tab.md).
  AICc needs more training observations than K + 1.

- ties:

  What to do when two candidates have exactly the same training
  criterion. The default `"error"` marks the fold as failed. `"first"`
  picks the tied candidate listed first in `mod_list`. No tolerance is
  applied, so only exact ties count.

- keep_fits:

  Logical. The default `FALSE` returns compact output. `TRUE` also
  returns the selected training fit of each fold in `fits`. This
  requires every candidate to have been run with `keep_fits = TRUE`,
  because records cannot rebuild a fit.

## Value

A cross-validation result of the same class as the candidates that
scores the selection procedure. Each fold's predictions, training
values, and diagnostics come from the candidate chosen in that fold, and
`settings$backend` is `"training_selection"`. The result adds
`selection`, which holds `chosen` (the candidate picked in each fold,
`NA` where selection failed), `tables` (each fold's training table, with
every candidate's K, log-likelihood, and criterion), and the
`criterion`, `ties`, and `candidates`. Results from
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
[`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md),
and
[`cross_validate_clogit()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_clogit.md)
also record the `metrics` scored and `metrics_dropped`, and each
candidate's `candidate_settings`. See
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
for the other elements.

## How the procedure is scored

Suppose your plan is to fit several candidate models and keep the one
with the lowest AICc. Scoring only the winner on held-out data would
flatter it, because it was chosen after seeing all of the data.
`select_cv()` scores the whole plan instead. In each fold it computes
the training criterion of every candidate from its fit to that fold's
analysis rows, picks the lowest, and keeps the held-out predictions that
candidate made. Assessment outcomes and CV scores never choose the
winner. The pooled score therefore evaluates the declared candidate set
and rule. If you then tune the candidate set or the rule using these
results, evaluate the final choice with an independent test set or an
outer validation layer.

## What the candidates must keep

Run every candidate with `keep_selection = TRUE` (compact records) or
`keep_fits = TRUE` (full fits). A record holds the native training
log-likelihood, the native and ecological parameter counts, the rows
that entered the likelihood, their responses, trials, and frequency
weights, the formula as text, the controls, and the original
diagnostics. NB2 and Gamma records also hold the estimated or fixed
nuisance parameter. A record holds no fitted model, design matrix,
raster, or formula environment, so it is much smaller than a fit,
although its row metadata still grows with the number of observations
and folds. Keeping records adds no fitting and does not change which
folds pass their diagnostics.

Before selecting, the function checks each record or fit against the
fold's analysis rows and the recorded responses and weights. If a
candidate keeps both records and fits, they must agree exactly. A record
that could not be built makes selection fail in that fold, without
changing the candidate's own result. Results that kept fits work without
records. A result that kept neither must be run again with
`keep_selection = TRUE`.

## Which candidates can be combined

The candidates must pass the same checks as in
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md),
with identical folds and observations and comparable likelihoods. NB2
and Poisson candidates can share a count target. Gamma and Gaussian
candidates can share a continuous target, scored with RMSE or MAE. An
information criterion cannot compare a count likelihood with a
continuous density. Conditional-choice candidates must also share their
partial-likelihood `method`.

For results from
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
[`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md),
and
[`cross_validate_clogit()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_clogit.md),
the procedure is scored on the metrics that every candidate recorded. A
default Poisson candidate records `"deviance"` and a default NB2
candidate does not, so their selection keeps `"rmse"`, `"mae"`, and
`"neg_log_score"`. `selection$metrics` lists the metrics kept and
`selection$metrics_dropped` those left out. Field and mixed-model
candidates must record the same metrics. Bayesian results have no
training likelihood and are rejected.

Penalties follow the native parameter counts, and scales and shapes
estimated by
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
count toward K.

- Conditional choice: the sample size in AICc and BIC is the number of
  choice sets (one event each), not the number of alternatives.

- Spatial fields: the coefficients, field SD, estimated range, and any
  observation dispersion count toward K. The integrated field values add
  nothing. The held-out joint block score is that of the model chosen in
  the fold.

- Mixed models: the sample size is the number of analysis observations.
  The covariance parameters count toward K, and the integrated group
  effects add nothing.

## Read the result

Read `complete` and `selection$chosen` first. A fold fails, and the
pooled scores become `NA`, when a criterion cannot be computed, is not
finite, or ties under `ties = "error"`. `selection$tables` shows the
training table behind each choice. For GLM and conditional-choice
results, the shared `settings` are copied from the first candidate,
named in `settings$controls_from`, and `selection$candidate_settings`
holds every candidate's controls. NB2 and Gamma selections also keep the
settings (`settings$fold_settings`) and `distribution` of the model
chosen in each fold, including when a Poisson or Gaussian reference
wins. GLM, conditional-choice, and mixed-model selections also keep
every candidate's warnings and diagnostics in `selection`, and
grouped-map labels (`raster_groups`) are kept for the chosen model. The
score does not integrate parameter or selection uncertainty, add
intervals, or establish calibration.

## See also

[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md),
[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
[`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md),
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md),
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md),
[`bic_tab()`](https://wpeterman.github.io/multiScaleR/reference/bic_tab.md),
and
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md).

## Examples

``` r
data("landscape_counts")
data("surv_pts")
# Three spatial folds. The simulated coordinates are in meters.
folds <- spatial_cv_folds(crds(vect(surv_pts)), block_size = 1000, k = 3,
                          crs = 32617)

# keep_selection = TRUE stores the training likelihoods that select_cv() reads
candidates <- list(
  null = cross_validate_glm(counts ~ 1, landscape_counts, folds,
                            family = "poisson", keep_selection = TRUE),
  site = cross_validate_glm(counts ~ site, landscape_counts, folds,
                            family = "poisson", keep_selection = TRUE))
selected <- select_cv(candidates, criterion = "AICc")
selected$selection$chosen         # model picked by training AICc in each fold
#>  fold1  fold2  fold3 
#> "site" "site" "null" 
selected$selection$tables$fold1   # the training table behind the first pick
#> 
#> Model selection based on AICc:
#> 
#>      K   AICc Delta_AICc AICcWt Cum.Wt      LL
#> site 2 287.62       0.00   0.53   0.53 -141.72
#> null 1 287.86       0.24   0.47   1.00 -142.90
#> 
selected$summary                  # held-out scores of the whole procedure
#>          metric estimate        unit n_assessment n_scored weight_sum complete
#> 1          rmse 2.236579 observation          100      100        100     TRUE
#> 2           mae 1.668719 observation          100      100        100     TRUE
#> 3 neg_log_score 2.174640 observation          100      100        100     TRUE
#> 4      deviance 2.382543 observation          100      100        100     TRUE
```
