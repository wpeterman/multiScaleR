# Simulate Prior Predictive Responses and Landscape Scales

Check what the complete set of priors implies before you fit. Each
simulated data set draws every parameter from its prior and then
generates responses for the prepared observations, including the scatter
of individual observations around their means.

## Usage

``` r
simulate_bayes_prior(object, nsim = 100L, seed = NULL)
```

## Arguments

- object:

  A
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  object. An object from `method = "prepare"` is enough. For a fitted
  object, the posterior draws are ignored and its specification supplies
  the same prior model.

- nsim:

  Positive integer number of prior predictive data sets; default 100.
  More data sets explore the priors more thoroughly and use more memory.

- seed:

  Optional positive integer simulation seed. Your own random-number
  state is restored afterward, even when simulation fails.

## Value

A list with these elements.

- `responses`: a data frame with one row per prepared observation and
  one column per simulated data set, named `sim_1`, `sim_2`, and so on.
  Values are counts for Poisson and NB2 models, success counts out of
  the prepared trial totals for binomial models, and response units for
  Gaussian models.

- `means`: a data frame of the matching expected responses. Binomial
  values are success probabilities, so multiply by the trial totals for
  expected counts. Other families use the response scale. Each mean is
  conditional on its own prior draw and excludes observation noise. It
  is not an expectation over the whole prior.

- `log_scales`: a named list with one matrix per landscape predictor.
  Rows are its distinct combinations of scale covariates and groups,
  labeled as in
  [`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md),
  and columns match the simulated data sets. Entries are natural log
  distances in projected map units, so exponentiate medians and interval
  limits to read them as distances.

- `observation_ids`: the prepared landscape ID of each observation, in
  the row order of `responses` and `means`.

- `priors`: the complete prior settings used, with defaults filled in.

## Read the simulations

Compare the simulated magnitudes, fractions of zeros, and scale
quantiles with ranges that are scientifically plausible for your design.
Extreme simulations show that a prior is loose for your design and may
call for a tighter setting, such as a smaller `beta_sd`. They are not
evidence of a poor posterior fit, and they do not choose a prior for
you.
[`vignette("v13-bayesian-scales", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v13-bayesian-scales.md)
works through this check.

Broad distance priors can have infinite arithmetic moments, so summarize
scales with log-scale summaries and distance quantiles. Broad response
group-SD priors with a log link can also give infinite prior predictive
count moments. A finite simulation does not show that a mean over prior
replicates exists.

## How the simulation works

Each data set uses one complete, independent prior parameter vector. All
observations in a group share that vector's group effect, which keeps
the dependence that the scale and response formulas specify. Scale
effects and response effects are drawn separately, as the model
specifies.

Simulation conditions on the prepared maps, fixed covariates, offsets,
group identities, binomial trial totals, and the frozen predictor
standardization. It uses no observed responses and no posterior draws,
and it creates no new sites or years. No optional Bayesian package is
needed. A nonfinite value stops with an error that names the replicate.
Values are never clipped and draws are never replaced.

## See also

[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
[`simulate.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/simulate.multiScaleR_bayes.md),
[`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md);
[`vignette("v13-bayesian-scales", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v13-bayesian-scales.md).

## Examples

``` r
example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
  package = "multiScaleR"))
prepared <- multiScale_bayes(count ~ habitat + factor(year) + (1 | site),
  example$inputs, example$data, scales = list(habitat = ~1 + (1 | year)),
  row_id = "site_year", priors = list(scale_median = 250), method = "prepare")
prior <- simulate_bayes_prior(prepared, nsim = 20, seed = 63)

# Fraction of zero counts in each prior data set, against the observed 0.04
range(vapply(prior$responses, function(y) mean(y == 0), numeric(1)))
#> [1] 0.38 0.96
# The default priors are loose and allow far more zeros than were observed.

# Prior range of sigma for each year, in meters
exp(apply(prior$log_scales$habitat, 1, stats::quantile, c(.025, .5, .975)))
#>        year=2021 year=2022 year=2023  year=2024  year=2025
#> 2.5%    25.29771  29.68201  41.06089   20.14086   16.99735
#> 50%    153.85049 124.09453 236.44030  142.68248  208.76500
#> 97.5% 1240.57462 636.67590 678.33440 1929.74700 1334.93399
# These are prior scale intervals, not posterior estimates or evidence of fit.
```
