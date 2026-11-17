# Compare Cross-Validated Prediction Scores

Rank candidate models by one held-out score, after checking that every
candidate was scored on the same folds and the same observations. Lower
scores are better.

## Usage

``` r
cv_tab(mod_list, metric, by_fold = FALSE)
```

## Arguments

- mod_list:

  Named list of complete cross-validation results to compare, for
  example `list(landscape = cv, null = null_cv)`. The names label the
  rows of the table and must be unique. Every candidate must use
  identical folds, responses, trials, and frequency weights, and the
  candidates must be of compatible kinds (see "What can be compared").

- metric:

  One score name, chosen explicitly: `"rmse"`, `"mae"`,
  `"neg_log_score"`, `"deviance"`, or `"brier"`. Every candidate must
  have recorded it. There is no default, because the scores answer
  different questions.

- by_fold:

  Logical. The default `FALSE` ranks the pooled scores. `TRUE` ranks the
  candidates separately within each fold, which shows whether a ranking
  is consistent across folds.

## Value

A data frame with one row per candidate (per candidate and fold when
`by_fold = TRUE`), sorted from best to worst. Its columns are:

- `Modnames`: the candidate's name in `mod_list`.

- `Rank`: 1 for the best score. Tied scores share the lowest rank.

- `Score`: the pooled score, or the fold score when `by_fold = TRUE`.

- `Delta`: `Score` minus the best score, in the units of the score. It
  is `NA` when every score is infinite.

- `Metric`, `Unit`: the score and what one score refers to, for example
  `"observation"`, `"trial"`, or `"choice_set"`. See
  [multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md).

- `N_assessment`, `N_scored`, `Weight_sum`: the planned assessment
  occurrences, the number that entered the score, and their total
  weight.

- `Warnings`: the number of warnings recorded while fitting the
  candidate (in that fold when `by_fold = TRUE`). Read them in the
  candidate's `warnings` element.

- `Fold`: the fold, when `by_fold = TRUE`.

The candidate results themselves are not modified.

## How to read the table

The best candidate is listed first. `Delta` is each candidate's distance
from it, in the units of the score. For RMSE and MAE that is a
difference in response units, such as animals per survey. For the
negative log score, a `Delta` of D means the best model gave the
observed outcomes exp(D) times the probability, as a geometric mean per
observation (per choice set for conditional choice). The per-fold table
from `by_fold = TRUE` is descriptive. Use it to see whether one fold
drives a ranking, never to pick a model for that fold from its own
assessment score.

## What can be compared

Every candidate must be complete, with no failed fold. An incomplete
candidate stops the comparison, because a model scored only on its
easier folds would look better than it is. No failed model or fold is
removed silently. Results of different kinds cannot be mixed.

|  |  |
|----|----|
| Results | Comparisons allowed |
| Gaussian, Poisson, and binomial from [`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md) and [`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md) | RMSE and MAE across families on identical observations. The log score, deviance, and Brier score need the same family and link. |
| NB2 or Gamma, alone or with Gaussian or Poisson | RMSE and MAE across these families on identical observations. Log scores between Poisson and NB2 on the same counts. A Gamma log score compares only with another Gamma. |
| [`cross_validate_clogit()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_clogit.md) | Other conditional-choice results with the same complete sets and partial-likelihood `method`, using `neg_log_score` only. |
| [`cross_validate_field()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_field.md) | Other field results with the same responses, folds, and projected coordinates. Poisson and NB2 share a count target. |
| [`cross_validate_mixed()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_mixed.md) | Other mixed-model results with the same groups, prediction target, and folds. |
| [`cross_validate_bayes()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_bayes.md) | Other Bayesian results with the same targets, folds, and blocks, using `neg_log_score` only. |

A binomial log score refers to the probability of the observed count of
successes, while a Brier score is a mean over individual trials. Field,
mixed-model, and Bayesian log scores are joint scores for blocks of
held-out observations, so they are kept apart from the per-observation
scores of ordinary GLMs.

## What the table does not do

The scores carry no penalty for model complexity, and `Delta` is not an
AIC difference, an Akaike weight, or a model probability. The spread of
fold scores is not a standard error, because folds share training data
and are few. If you use the table to choose a model or a tuning rule,
the winning score is optimistic. Evaluate that choice with an
independent test set or an outer validation layer. To score a
predeclared rule that picks a model by training AIC, AICc, or BIC, use
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md).

## See also

[`cross_validate()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate.md),
[`cross_validate_glm()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_glm.md),
[`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md),
[multiScaleR_cv](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_cv-class.md),
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md),
and
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md).

## Examples

``` r
data("landscape_counts")
data("surv_pts")
# Three spatial folds. The simulated coordinates are in meters.
folds <- spatial_cv_folds(crds(vect(surv_pts)), block_size = 1000, k = 3,
                          crs = 32617)
candidates <- list(
  null = cross_validate_glm(counts ~ 1, landscape_counts, folds,
                            family = "poisson"),
  site = cross_validate_glm(counts ~ site, landscape_counts, folds,
                            family = "poisson"))

# Pooled ranking. Delta is the gap from the best model in animals per survey.
cv_tab(candidates, metric = "rmse")
#>   Modnames Rank    Score      Delta Metric        Unit N_assessment N_scored
#> 1     null    1 2.210604 0.00000000   rmse observation          100      100
#> 2     site    2 2.232857 0.02225306   rmse observation          100      100
#>   Weight_sum Warnings
#> 1        100        0
#> 2        100        0

# exp(Delta) is how many times more probability the best model gave the
# observed counts, as a geometric mean per survey
cv_tab(candidates, metric = "neg_log_score")
#>   Modnames Rank    Score       Delta        Metric        Unit N_assessment
#> 1     null    1 2.163750 0.000000000 neg_log_score observation          100
#> 2     site    2 2.169778 0.006027989 neg_log_score observation          100
#>   N_scored Weight_sum Warnings
#> 1      100        100        0
#> 2      100        100        0

# Is the ranking the same in every fold?
cv_tab(candidates, metric = "rmse", by_fold = TRUE)
#>   Modnames Rank    Score       Delta Metric        Unit N_assessment N_scored
#> 1     null    1 2.017854 0.000000000   rmse observation           34       34
#> 2     site    2 2.091114 0.073260576   rmse observation           34       34
#> 3     null    1 1.608284 0.000000000   rmse observation           33       33
#> 4     site    2 1.620645 0.012360460   rmse observation           33       33
#> 5     site    1 2.824229 0.000000000   rmse observation           33       33
#> 6     null    2 2.833142 0.008912299   rmse observation           33       33
#>   Weight_sum Warnings  Fold
#> 1         34        0 fold1
#> 2         34        0 fold1
#> 3         33        0 fold2
#> 4         33        0 fold2
#> 5         33        0 fold3
#> 6         33        0 fold3
```
