# Rank models by BIC

Build a model-selection table ranked by BIC (the Bayesian information
criterion) for a list of fitted models. As in
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md),
each estimated scale of effect (sigma) and kernel shape counts as a
parameter, and lists can mix `multiScaleR` objects with plain fitted
models.

## Usage

``` r
bic_tab(mod_list,
        mod_names = NULL,
        verbose = FALSE,
        ...)
```

## Arguments

- mod_list:

  A list of fitted models: `multiScaleR` objects from
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
  plain fitted models, or both. Every model must be fitted to the same
  observations with a comparable likelihood, including matching
  responses, weights, and likelihood measure. Bayesian fits,
  cross-validation results, and private spatial prototypes are rejected.

- mod_names:

  Optional character vector of model names, one per model in `mod_list`.
  By default, `multiScaleR` objects are labeled
  `[kernel=...; family=...]` followed by the right-hand side of the
  model formula.

- verbose:

  Logical. If `TRUE`, the table is also printed to the console. Default
  `FALSE`.

- ...:

  Not used.

## Value

A data frame of class `"bictab"` (from the AICcmodavg package), sorted
from the best (smallest) BIC, with one row per model:

- `Modnames`: model name.

- `K`: number of estimated parameters, counted as in
  [`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md):
  the fitted likelihood parameters (including ancillary ones such as
  negative binomial dispersion, or native effective degrees of freedom
  for GAMs) plus each estimated scale and shape of a `multiScaleR`
  model.

- `BIC`: the Bayesian information criterion. Smaller is better.

- `Delta_BIC`: difference from the best model.

- `BICWt`: BIC weight, the relative support for each model within this
  set. Weights sum to 1.

- `Cum.Wt`: cumulative BIC weight.

- `LL`: maximized log-likelihood.

A `likelihood_adjustments` attribute is added when a likelihood was
adjusted, as described in
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md).

## Details

BIC penalizes each parameter by the log of the sample size rather than
by 2, so it penalizes complexity more heavily than AIC once there are
more than about 7 observations. It is often preferred when the goal is
the single best-supported model rather than model averaging.

Every rule in
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md)
applies here: the parameter count (including all estimated components of
`unmarked` models, spaMM degrees of freedom, and the checked negative
binomial BAM adjustment), the comparability checks (responses, rows,
weights, binomial outcomes, ML and REML, partial likelihoods, used
`pcount()` visits, and distance-removal protocols), the sample-size
conventions (retained site rows for `unmarked` models and native event
counts for Cox and conditional logistic models), and the rejected
models. Inapplicable hard-radius shape positions are excluded from `K`.
Calling the criterion BIC does not make a posterior target a sampling
likelihood, so Bayesian fits are rejected here too.

The table itself is built with
[`AICcmodavg::bictabCustom()`](https://rdrr.io/pkg/AICcmodavg/man/bictabCustom.html).

## See also

[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md),
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
and
[`vignette("v02-user-guide", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v02-user-guide.md).

## Author

Bill Peterman

## Examples

``` r
## Bundled counts simulated with a habitat scale of 75 map units
data("pts")
data("count_data")
hab <- terra::rast(system.file("extdata", "hab.tif", package = "multiScaleR"))
prepared <- kernel_prep(pts, hab, max_D = 300, sigma = 75, verbose = FALSE)
dat <- data.frame(y = count_data$y, prepared$kernel_dat)
fixed_model <- glm(y ~ hab, family = poisson(), data = dat)
null_model <- glm(y ~ 1, family = poisson(), data = dat)
spatial_model <- multiScale_optim(fixed_model, prepared, par = 75 / 300,
                                  n_cores = 1, verbose = FALSE)

## BIC charges log(75), about 4.3, for the estimated sigma, versus 2 in AIC
bic_tab(list(spatial_model, fixed_model, null_model),
        mod_names = c("spatial", "fixed", "null"))
#> 
#> Model selection based on BIC:
#> 
#>         K    BIC Delta_BIC BICWt Cum.Wt      LL
#> fixed   2 236.98      0.00   0.9    0.9 -114.17
#> spatial 3 241.29      4.31   0.1    1.0 -114.17
#> null    1 349.01    112.03   0.0    1.0 -172.34
#> 
```
