# Check Sampling Precision of Fitted Bayesian Scales

Check chain mixing and simulation precision for the quantities you
report from a sampled fit, the log scale of each landscape predictor and
the between-group SDs. These derived quantities can mix differently from
the individual parameters checked by the sampler's own diagnostics.

## Usage

``` r
diagnose_bayes_scales(object, probs = c(0.025, 0.5, 0.975))
```

## Arguments

- object:

  A fitted
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  object. New diagnostics need the original unpermuted post-warmup
  chains from NUTS (`sampler`) or slice sampling (`chains`). Prepared,
  approximate, and compact fits return an unavailable report.

- probs:

  Numeric vector of distinct probabilities strictly between 0 and 1. The
  default `c(0.025, 0.5, 0.975)` gives the median and the limits of a
  central 95% interval. Changing it changes which quantiles have their
  simulation error reported, not any diagnostic threshold.

## Value

A list with these elements.

- `available`, `reason`: whether full chain diagnostics could be
  computed, and an explanation when they could not.

- `raw_diagnostics`: the sampler or approximation diagnostics stored in
  the fit, unchanged.

- `parameters`: a data frame with one row per derived quantity.
  `variable` is an internal label. `quantity` is `"log_distance"`,
  `"scale_group_sd"`, or `"response_group_sd"`. `covariate`, `group`,
  and a representative observation `row` identify it, and `units` states
  its units. `rhat`, `ess_bulk`, and `ess_tail` give the rank-normalized
  Rhat and the bulk and tail effective sample sizes.

- `quantiles`: a data frame with the same identifying columns plus
  `probability`, the quantile `estimate`, and its Monte Carlo standard
  error `mcse`, both in the stated units. Undefined values are `NA`.

- `chains`, `iterations`: the number of retained chains and of
  post-warmup iterations per chain, or `NA` when unavailable.

- `derived_passed`: `TRUE` when there are at least two chains and every
  derived quantity has Rhat at most 1.01 and bulk and tail ESS of at
  least 100 per chain. `NA` when unavailable.

- `passed`: `TRUE` when both the original sampler checks and the derived
  checks pass. `NA` when unavailable. For a slice fit sampled without
  the posterior package, the original parameter checks are computed here
  from the stored chains, while the fit's stored diagnostics stay as
  recorded.

- `criteria`, `notes`: the fixed derived checks and their limits.

## Read the report

|  |  |  |
|----|----|----|
| Column | Acceptable | If it fails |
| `rhat` | 1.01 or less | Refit with another `seed` and compare. A weak landscape effect leaves its scale poorly identified, and longer chains do not cure that |
| `ess_bulk`, `ess_tail` | At least 100 per chain (400 for four chains) | Increase `iter`, or switch from slice sampling to NUTS |
| `mcse` | No fixed cutoff. Small next to the precision your question needs | Increase `iter`, and compare independent runs when an interval limit matters |

Landscape scales are diagnosed on the natural log axis in map units,
which avoids overflow and the infinite arithmetic distance moments that
the hierarchical priors allow (see
[`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md)).
Exponentiate a log-distance quantile to express it in map units. Its
MCSE is on the log axis, so an MCSE of 0.01 is about 1% of the distance.
The MCSE is computational uncertainty, not a posterior interval or a
biological SD. Scale-group SDs measure population variation in log
distance. Response-group SDs use the response linear-predictor axis.

The posterior package computes rank-normalized Rhat, effective sample
sizes, and quantile MCSE. A passed screen guarantees neither a precise
tail endpoint nor discovery of every mode. Compare the MCSE with the
precision your scientific decision needs, inspect independent runs when
tails matter, and assess the likelihood and priors separately.

## Why compact fits cannot be checked

Permuted or thinned stored draws cannot rebuild the original chain order
or full-chain precision. Compact examples therefore keep their original
diagnostic evidence but cannot provide new chain diagnostics. Laplace
draws and weighted importance candidates are not MCMC chains, so use the
stored numerical checks for those fits.

## Draw chain plots

For chain plots in bayesplot, pass `fit$chains` for a slice fit, or
`rstan::extract(fit$sampler, permuted = FALSE, inc_warmup = FALSE)` for
a NUTS fit.
[`vignette("v13-bayesian-scales", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v13-bayesian-scales.md)
shows trace plots alongside these diagnostics, and
[`vignette("v15-bayesian-computation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v15-bayesian-computation.md)
works through a fit that fails them.

## See also

[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
[`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md),
[`plot.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/plot.multiScaleR_bayes.md),
[`diagnose_bayes_probabilities()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_probabilities.md),
[multiScaleR_bayes](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md)

## Examples

``` r
# A stored slice fit that keeps its four original chains
shared <- readRDS(system.file("extdata", "cache-bayesian-scales.rds",
  package = "multiScaleR"))$fit
if (requireNamespace("posterior", quietly = TRUE)) {
  report <- diagnose_bayes_scales(shared)
  print(report$passed)   # sampler and derived checks together
  # Agreement among chains and effective sample size for each log sigma
  print(report$parameters[c("covariate", "rhat", "ess_bulk", "ess_tail")])
  # Interval limits of log sigma, with their Monte Carlo standard errors
  q <- report$quantiles
  print(data.frame(q[c("covariate", "probability")],
                   sigma_m = exp(q$estimate), mcse_log = q$mcse))
}
#> [1] TRUE
#>   covariate     rhat  ess_bulk ess_tail
#> 1     land1 1.002665  905.3584  937.121
#> 2     land2 1.003975 1033.0974 1611.410
#>   covariate probability  sigma_m    mcse_log
#> 1     land1       0.025 180.7061 0.012290191
#> 2     land1       0.500 246.6455 0.005139396
#> 3     land1       0.975 344.4668 0.027488507
#> 4     land2       0.025 412.5682 0.006516360
#> 5     land2       0.500 501.1973 0.003116837
#> 6     land2       0.975 599.8042 0.005060544

# A compact fit keeps thinned, permuted draws and cannot be checked
example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
                              package = "multiScaleR"))
compact <- diagnose_bayes_scales(example$fits$pooled)
compact$available  # FALSE: the compact example contains no full sampler
#> [1] FALSE
compact$reason
#> [1] "Full unpermuted MCMC chains are not retained. Stored permuted or thinned draws cannot recover original chain diagnostics."
```
