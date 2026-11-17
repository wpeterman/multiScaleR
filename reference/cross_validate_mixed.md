# Cross-Validate Mixed-Effects Models

Refit a mixed model in every analysis split and score the observations
held out from it. A *random intercept* gives each group, such as a site,
its own baseline drawn from a common distribution, which accounts for
the similarity among repeated surveys of that group. Choose whether
prediction is for groups represented in training or for wholly new
groups. Response means and joint group log scores integrate over the
shared group effect.

## Usage

``` r
cross_validate_mixed(
  formula,
  data,
  folds,
  family = c("gaussian", "poisson", "binomial", "nbinom2", "Gamma"),
  prediction = c("new_group", "existing_group"),
  kernel_inputs = NULL,
  optim_args = list(),
  metrics = NULL,
  control = NULL,
  integration_tol = 1e-08,
  keep_fits = FALSE
)
```

## Arguments

- formula:

  Two-sided native mixed-model formula with one grouping column, for
  example `count ~ habitat + (1 | site)`. Gaussian models also support
  one block of numeric random slopes such as `(1 + habitat | site)`.
  Other families allow one random intercept. Use explicit predictors (no
  `.`) and formula offsets. The response is a data column, or
  `cbind(successes, failures)` for grouped binomial data with positive
  trial counts.

- data:

  Raw data frame with unique row names. The names identify rows in the
  split plan and, if supplied, in the prepared raw-cell buffers.

- folds:

  Named list of explicit analysis and assessment row-ID splits. Use
  [`grouped_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/grouped_cv_folds.md)
  or
  [`spatial_cv_folds()`](https://wpeterman.github.io/multiScaleR/reference/spatial_cv_folds.md)
  to hold out entire groups. Use explicit splits within groups, such as
  later visits, for existing-group prediction.

- family:

  Character response family. Choose `"gaussian"` (identity, default),
  `"poisson"` (log), `"binomial"` (logit), `"nbinom2"` (log), or
  `"Gamma"` (log).

- prediction:

  Character prediction target. `"new_group"` (default) predicts groups
  absent from training, so each assessment group must have no analysis
  rows in that split, and its prediction averages over the fitted
  distribution of group effects. `"existing_group"` predicts new
  observations of groups seen in training, so every assessment group
  must have analysis rows in that split, and its prediction conditions
  on them.

- kernel_inputs:

  Optional raw-cell `multiScaleR_data` object for estimating the scales
  of effect from the analysis rows by native refits. The default `NULL`
  fits the supplied formula without scale optimization.

- optim_args:

  Named list of extra serial
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  controls. Use only with `kernel_inputs`. `profile_control` and the
  arguments this function sets itself are reserved.

- metrics:

  Character vector of `"neg_log_score"`, `"rmse"`, `"mae"`, or `"brier"`
  (binomial only). By default, report the joint group log score and
  RMSE, or the joint group log score and Brier score for binomial
  responses. The log score is summed across assessment groups and
  divided by the number of assessment rows. Binomial RMSE and MAE
  compare the expected success probability with the observed fraction of
  successes. The Brier score averages the squared probability error over
  trials.

- control:

  A
  [`glmmTMB::glmmTMBControl()`](https://rdrr.io/pkg/glmmTMB/man/glmmTMBControl.html)
  list for the native analysis fits. The default `NULL` uses
  `glmmTMBControl()`.

- integration_tol:

  Relative tolerance for the one-dimensional adaptive integration over
  the group effect in non-Gaussian models. The default is `1e-8`.

- keep_fits:

  Logical. Retain the fitted model of each successfully scored fold.
  Defaults to `FALSE`.

## Value

A list of class `multiScaleR_mixed_cv`, which inherits from
`multiScaleR_cv`. `summary` and `metrics` hold the pooled and per-fold
scores, `predictions` the integrated response mean of each held-out row,
`joint_scores` the joint log score of each held-out group, and
`training` the coefficients, random-effect covariance, dispersion, and
any fitted scales of each fold. Compact `selection_records` are always
kept. A failed fold makes the pooled scores `NA`. The mixed model
section of
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md)
describes every element.

## Details

This function is experimental. Its arguments and results may change in a
later release without a deprecation period.

## How predictions and scores are made

Mixed CV uses maximum-likelihood fits with constant dispersion and a
Gaussian distribution for the group effects. Given the fitted
parameters, new groups use the estimated population distribution of
effects. Existing groups condition that distribution on their analysis
responses. Gaussian response models use exact multivariate normal
conditioning. For other families the function integrates over one group
intercept with adaptive quadrature. Each assessment group contributes
one joint predictive density, because all of its rows share one unknown
group effect. Parameter, ecological-scale, and model-selection
uncertainty are outside these plug-in predictions. Native fitting,
preprocessing, and selection use the analysis observations only.

## Compare and select models

Use
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
to compare complete results with the same observation support, group
assignments, prediction target, and splits. A joint group score is not
on the same footing as the per-observation scores of
[`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md),
so the two cannot be mixed. Use
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md)
to evaluate a training information-criterion selection rule. Its sample
size is the native number of analysis observations. Estimated covariance
parameters count toward the penalty, and integrated group effects do not
count individually.
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md)
shows new-site and later-visit designs.

## See also

[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
[`cross_validate_field()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_field.md),
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md),
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md),
and
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md).

## Examples

``` r
if (requireNamespace("glmmTMB", quietly = TRUE)) {
  surveys <- glmmTMB::Salamanders
  # Hold out whole sites, so every prediction is for a new site
  folds <- grouped_cv_folds(surveys$site, k = 3, ids = rownames(surveys))
  # One fold keeps the example fast
  cv <- cross_validate_mixed(count ~ spp + (1 | site), surveys,
    folds[1], family = "poisson")
  cv$summary   # joint log score per held-out row, and RMSE in counts
}
#>          metric estimate                        unit n_assessment n_scored
#> 1 neg_log_score 1.930772 joint_group_density_per_row          224      224
#> 2          rmse 3.293890                    response          224      224
#>   weight_sum complete
#> 1        224     TRUE
#> 2        224     TRUE
```
