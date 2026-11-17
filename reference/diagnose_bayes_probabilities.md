# Check Monte Carlo Precision of Scale Probabilities

Estimate the posterior probability that each fitted kernel scale lies at
or below thresholds you name, and report how precisely the original
sampling chains pin down each probability.

## Usage

``` r
diagnose_bayes_probabilities(
  object,
  thresholds,
  scale = c("distance", "log"),
  mcse_target = NULL
)
```

## Arguments

- object:

  A fitted
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  object that keeps its original, unpermuted post-warmup chains from
  NUTS (`sampler`) or slice sampling (`chains`). Compact, thinned, and
  approximate fits return an unavailable report, because permuted draws
  cannot recover chain order.

- thresholds:

  A nonempty named list of numeric vectors. Each name selects a fitted
  landscape predictor, and each vector holds distinct finite scale
  thresholds. Every threshold is applied to every distinct fitted
  combination of scale covariates and groups for that predictor. For
  example, `list(habitat = c(100, 250))` checks two scales in projected
  map units.

- scale:

  `"distance"` (default) means `thresholds` are positive kernel scales
  in map units. `"log"` means they are natural logs of those scales.
  Thresholds refer to the scale parameters summarized by
  [`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md),
  not to transformed kernel-weight radii such as those from
  [`kernel_dist()`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md).

- mcse_target:

  Optional target for the Monte Carlo standard error (MCSE) on the
  probability axis, strictly between 0 and 1. For example, `0.005` asks
  for an MCSE no larger than half a percentage point. The default `NULL`
  reports precision without judging it against a target.

## Value

A list with these elements.

- `available`, `reason`: whether the original chains could be read, and
  an explanation when they could not. A single event can lack precision
  estimates even when the report is available.

- `raw_diagnostics`: the fit's original diagnostics, unchanged.

- `probabilities`: a data frame with one row per fitted scale
  combination and threshold. `covariate`, `group`, `row`, `threshold`,
  and `threshold_scale` identify the event. `estimate` is the posterior
  probability that the scale is at or below the threshold. `mcse`,
  `ess_mean`, and `rhat` are the Monte Carlo standard error, mean
  effective sample size, and rank-normalized Rhat of the event
  indicator. `constant` flags an indicator that never changed.
  `precision_available` is `TRUE` when the MCSE and mean ESS are finite
  and positive. `precision_target_met` compares `mcse` with
  `mcse_target`, and is `NA` when no target was given. `reason` explains
  missing precision. `rhat_available` and `rhat_reason` report
  separately whether Rhat could be computed. Undefined values are `NA`.

- `chain_probabilities`: the same event identifiers with `chain`, that
  chain's `estimate`, and whether its indicator is `constant`. Large
  differences among chains signal poor mixing.

- `chains`, `iterations`: the number of chains and of post-warmup
  iterations per chain, or `NA` when unavailable.

- `mcse_target`, `event`, `notes`: the requested target, the event
  definition (`"kernel_scale <= threshold"`), and the limits of
  interpretation.

- `metadata_establishes_eligibility`: always `FALSE`. The report does
  not change the fit's sampling or derived-scale diagnostic gates.

## Read the precision

The MCSE is computational uncertainty in the probability estimate, not
posterior uncertainty about the scale and not observation variation.
`precision_target_met` only compares that MCSE with your target. It does
not override Rhat, ESS, the raw sampler checks, or model checks, and a
small MCSE cannot rule out an unvisited tail or mode. Read
[`diagnose_bayes_scales()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_scales.md)
alongside this report, and compare independent runs when a consequential
probability remains uncertain. Longer chains reduce the MCSE.

## How the events are computed

Each event is evaluated on fitted log scales in map units, which avoids
overflow from exponentiating distances. The posterior package computes
the mean ESS, MCSE, and rank-normalized Rhat of each binary indicator.
At least two chains with four retained iterations each are needed for
precision diagnostics. A constant indicator has no estimable precision,
even when its probability is exactly zero or one. Indicator Rhat can be
undefined for a balanced binary event even when the MCSE is finite. It
is then reported as `NA` with its own flag and reason, which alone does
not make the precision unavailable. A nonfinite derived log scale is an
error, and iterations are never dropped.

The function does not fit a model, change priors, simulate new groups,
or treat importance candidates as MCMC chains. Reading a full NUTS
sampler can take substantial memory, but the report keeps only
summaries.
[`vignette("v15-bayesian-computation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v15-bayesian-computation.md)
shows how to report a probability with its MCSE.

## See also

[`diagnose_bayes_scales()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_scales.md),
[`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md),
[`predict_scales()`](https://wpeterman.github.io/multiScaleR/reference/predict_scales.md),
[`diagnose_bayes_laplace()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_laplace.md),
[multiScaleR_bayes](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md)

## Examples

``` r
# A stored slice fit that keeps its four original chains
shared <- readRDS(system.file("extdata", "cache-bayesian-scales.rds",
  package = "multiScaleR"))$fit
if (requireNamespace("posterior", quietly = TRUE)) {
  report <- diagnose_bayes_probabilities(shared,
    thresholds = list(land1 = c(200, 300)), mcse_target = 0.01)
  # Probability that sigma for land1 is at most 200 m and at most 300 m
  print(report$probabilities[c("threshold", "estimate", "mcse",
                               "precision_target_met")])
  print(report$chain_probabilities[c("threshold", "chain", "estimate")])
}
#>   threshold estimate        mcse precision_target_met
#> 1       200  0.08175 0.007034019                 TRUE
#> 2       300  0.90425 0.008701902                 TRUE
#>   threshold chain estimate
#> 1       200     1    0.059
#> 2       200     2    0.078
#> 3       200     3    0.089
#> 4       200     4    0.101
#> 5       300     1    0.898
#> 6       300     2    0.903
#> 7       300     3    0.893
#> 8       300     4    0.923

# A compact fit keeps thinned, permuted draws and cannot be checked
example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
  package = "multiScaleR"))
compact <- diagnose_bayes_probabilities(example$fits$pooled,
  thresholds = list(habitat = c(100, 250)))
compact$available   # FALSE
#> [1] FALSE
compact$reason
#> [1] "Full unpermuted MCMC chains are not retained; stored permuted or thinned draws cannot recover them."
```
