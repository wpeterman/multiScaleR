# Assess Fitted Scale Probabilities from an Uncorrected Laplace Fit

Check the probability that each fitted scale lies at or below a distance
you name, as reported by a joint Laplace or uncorrected marginal Laplace
fit, against a fresh importance-weighted estimate of the same posterior.
Read this check alongside the fit's own optimizer and integration
diagnostics.

## Usage

``` r
diagnose_bayes_laplace(
  object,
  thresholds,
  scale = c("distance", "log"),
  tolerance = 0.05,
  candidates = 16384L,
  seed = NULL
)
```

## Arguments

- object:

  A fitted
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  object from `method = "laplace"`, or from `method = "marginal"`
  without importance correction, with its approximation draws and
  geometry intact. The check needs the RTMB and posterior (version 1.6.0
  or later) packages.

- thresholds:

  Named list of scale thresholds, as in
  [`diagnose_bayes_probabilities()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_probabilities.md).
  Each name selects a fitted landscape predictor, and each vector holds
  distinct finite thresholds, such as `list(habitat = c(150, 350))`.
  Every threshold is checked for every distinct fitted combination of
  scale covariates and groups, such as each year.

- scale:

  Units of `thresholds`. `"distance"` (default) uses map units and needs
  positive values. `"log"` uses natural log map units. The event is a
  fitted kernel scale less than or equal to the threshold.

- tolerance:

  Absolute difference in probability that matters for your inference,
  strictly between 0 and 1; default 0.05. Choose it, and the thresholds,
  before looking at the results.

- candidates:

  Number of fresh importance candidates, at least 2000; default 16384.
  The validation study used 65536. More candidates can reduce simulation
  error when the proposal overlaps the posterior well. Every candidate
  is kept for the weight and event diagnostics, so memory use grows with
  this number and with the number of parameters.

- seed:

  Optional positive integer seed for the assessment. Your own
  random-number state is restored afterward.

## Value

A list with these elements.

- `available`: `TRUE` when importance ratios and reference summaries
  could be calculated. A reference that fails its overlap checks is
  still available, but its statuses stay `"unresolved"`.

- `reason`: why the reference is unavailable, why it failed its overlap
  checks, or that the original fit failed its own diagnostics. An empty
  string otherwise.

- `raw_diagnostics`: the fit's original `diagnostics`, unchanged.

- `overlap`: the importance checks of the reference (`passed`,
  `pareto_k`, `ess`, `min_tail_ess`, `max_weight`, `split_cdf_error`),
  or `NULL` when the reference is unavailable.

- `target_evaluation`: counts of posterior density evaluations and of
  rechecks of invalid values.

- `probabilities`: a data frame with one row per fitted scale
  combination and threshold. `covariate`, `group`, `row`, `threshold`,
  and `threshold_scale` identify the event. `approximation` is the
  probability from the fit's draws and `reference` the
  importance-weighted probability, with its Monte Carlo standard error
  `reference_mcse`. `difference` is the approximation minus the
  reference, and `difference_lower` and `difference_upper` bound it
  numerically. `indicator_ess` is the effective sample size of the
  reference event indicator, `split_reference_difference` the absolute
  difference between the two halves of the reference, and `status` the
  verdict described in "Read the status".

- `candidates`, `seed`, `tolerance`: the settings used.

- `coefficient_proposal`: the number of fallback coefficient proposals
  (`fallbacks`) and the largest remaining Newton step (`max_step`).

- `defensive_draws`: the number of candidates drawn from the defensive
  prior component.

- `metadata_establishes_eligibility`: always `FALSE`.

## Details

This function is experimental. Its arguments and results may change in a
later release without a deprecation period.

## Read the status

|  |  |  |
|----|----|----|
| `status` | Meaning | What to do |
| `"discrepancy"` | The numerical bounds on the difference lie wholly beyond `tolerance` on one side | Do not report the approximate probability. Sample the model |
| `"within_tolerance"` | Both bounds lie inside `tolerance` | You can report the approximate probability, for this threshold only |
| `"unresolved"` | Every other case, including a reference too imprecise or overlapping the posterior too poorly to decide | Treat the probability as unchecked. Raise `candidates`, or sample |

An event can be resolved only when the reference passes the package's
importance overlap checks, the original fit passed its own diagnostics,
the event has an indicator ESS of at least 400, and the two independent
halves of the reference differ by no more than the smaller of 0.05 and
`tolerance`. The bounds combine three standard errors of the weighted
reference with Wilson bounds for the original approximate draws. A
constant reference indicator also leaves the event unresolved. Expect
`"unresolved"` often, and read it as "not checked", not as "close
enough".
[`vignette("v15-bayesian-computation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v15-bayesian-computation.md)
reports how often each status occurred in simulations.

## How the reference is built

Candidates come 80% from heavy-tailed local proposals around the
approximation and 20% from a defensive component based on the
latent-effect and scale priors. Response coefficients are drawn from
their exact Gaussian conditional, or from a normalized heavy-tailed
conditional proposal for other families. The full posterior density and
the actual mixture density determine every weight. Neither the
approximation nor its draws are replaced.

This exploratory check targets only the requested scale probabilities. A
heavy-tailed proposal can still miss posterior regions. The report does
not certify approximation accuracy, general calibration, held-out
predictions, or the fit's original diagnostics. Compare with a full
posterior fit whenever your conclusions depend on these probabilities.

## See also

[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
[`diagnose_bayes_probabilities()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_probabilities.md),
[`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md),
[multiScaleR_bayes](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md);
[`vignette("v15-bayesian-computation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v15-bayesian-computation.md).

## Examples

``` r
# \donttest{
if (requireNamespace("RTMB", quietly = TRUE) &&
    requireNamespace("posterior", quietly = TRUE)) {
  example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
    package = "multiScaleR"))
  # A joint Laplace fit of the annual model, which takes about a second
  fast_fit <- multiScale_bayes(count ~ habitat + factor(year) + (1 | site),
    example$inputs, example$data, scales = list(habitat = ~1 + (1 | year)),
    row_id = "site_year", family = "nbinom2",
    priors = list(scale_median = 250, scale_group_sd = 0.5),
    method = "laplace", seed = 527001)

  # Probability that each year's sigma is 150 m or less. The smallest
  # candidate budget keeps this example short. The validation study used
  # 65536 candidates.
  check <- diagnose_bayes_laplace(fast_fit,
    thresholds = list(habitat = 150), candidates = 2000, seed = 127)
  print(check$reason)   # empty when the reference and the fit passed
  print(check$probabilities[c("group", "approximation", "reference",
                              "status")])
  # "unresolved" means the probability was not checked, not that it is
  # close enough.
}
#> [1] "The reference did not pass the posterior importance overlap checks."
#>       group approximation   reference     status
#> 1 year=2021        0.0140 0.001640765 unresolved
#> 2 year=2022        0.0175 0.002140476 unresolved
#> 3 year=2023        0.2050 0.108813215 unresolved
#> 4 year=2024        0.0175 0.003687778 unresolved
#> 5 year=2025        0.1465 0.002778998 unresolved
# }
```
