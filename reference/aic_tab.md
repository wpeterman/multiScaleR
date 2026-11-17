# Rank models by AIC or AICc

Build a model-selection table ranked by AICc (the default) or AIC for a
list of fitted models. Each estimated scale of effect (sigma), and each
estimated kernel shape, counts as a parameter, so a scale-optimized
model competes fairly with fixed-scale or non-spatial alternatives in
the same list. Lists can mix `multiScaleR` objects with plain fitted
models such as `glm` objects.

[`vignette("v02-user-guide", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v02-user-guide.md)
compares kernels and model structures with these tables, including an
N-mixture model.

## Usage

``` r
aic_tab(mod_list,
        AICc = TRUE,
        mod_names = NULL,
        verbose = FALSE,
        ...)
```

## Arguments

- mod_list:

  A list of fitted models: `multiScaleR` objects from
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
  plain fitted models (for example `glm` or `lm`), or both. Every model
  must be fitted to the same observations with a comparable likelihood
  (see "Check that models are comparable"). Bayesian fits,
  cross-validation results, and private spatial prototypes are rejected.

- AICc:

  Logical. If `TRUE` (the default), rank by the second-order corrected
  AIC (AICc). Set `FALSE` for standard AIC. AICc is recommended when the
  number of observations is small relative to the number of parameters
  (n / K \< 40 is a common rule of thumb). AICc requires n \> K + 1 for
  every model and raises an error otherwise. Use `AICc = FALSE` in that
  case.

- mod_names:

  Optional character vector of model names, one per model in `mod_list`.
  By default, `multiScaleR` objects are labeled
  `[kernel=...; family=...]` followed by the right-hand side of the
  model formula. `unmarked` labels name the model class and, where
  present, the mixture, so Poisson, negative binomial, and ZIP
  alternatives get distinct labels.

- verbose:

  Logical. If `TRUE`, the table is also printed to the console. Default
  `FALSE`.

- ...:

  Not used.

## Value

A data frame of class `"aictab"` (from the AICcmodavg package), sorted
from the best (smallest) criterion value, with one row per model:

- `Modnames`: model name, from `mod_names` or the default label.

- `K`: number of estimated parameters. It counts the fitted likelihood
  parameters, including ancillary ones such as negative binomial
  dispersion (or native effective degrees of freedom for GAMs), plus
  each estimated scale and shape of a `multiScaleR` model.

- `AICc` or `AIC`: the information criterion. Smaller is better.

- `Delta_AICc` or `Delta_AIC`: difference from the best model. Models
  within about 2 units of the best are generally considered to have
  comparable support.

- `ModelLik`: relative likelihood of the model, `exp(-0.5 * Delta)`.

- `AICcWt` or `AICWt`: Akaike weight, the relative support for each
  model within this set. Weights sum to 1.

- `LL`: maximized log-likelihood.

- `Cum.Wt`: cumulative Akaike weight.

When a likelihood was adjusted (see "Count parameters"), the table
carries a `likelihood_adjustments` attribute.

## Count parameters

For a `multiScaleR` object, the log-likelihood and parameter count come
from the refit model in `opt_mod`, and each estimated sigma (and shape,
for `kernel = "expow"`) is added to `K`. The AIC printed by the native
model in `opt_mod` leaves those parameters out. Hard-radius metrics in
an `"expow"` model add their estimated radius but no shape, and
fixed-radius metrics add neither. Plain models in the list get no sigma.
Otherwise the native log-likelihood degrees of freedom are used when
available, so negative binomial dispersion counts as one parameter.
Missing or inconsistent parameter metadata and non-finite
log-likelihoods raise an error rather than guessing from predictor
names.

- **`unmarked` models.** The native AIC and log-likelihood determine the
  count. It includes every abundance or occupancy coefficient, every
  detection coefficient, factor contrasts, mixture parameters, and
  estimated random-effect standard deviations. Optimizer and TMB
  free-parameter metadata are checked for agreement when present. Latent
  abundances and realized random effects add nothing. For example, four
  abundance coefficients, three detection coefficients, and one
  random-effect standard deviation give eight native parameters, or nine
  with one estimated Gaussian scale. A table of fixed coefficients can
  therefore show fewer parameters than the fitted count.

- **The `K` argument of
  [`unmarked::pcount()`](https://ecoverseR.github.io/unmarked/reference/pcount.html)**
  is an upper bound for the abundance sum in the likelihood, not a
  parameter count. Check that the log-likelihood does not change when
  you refit with a larger bound.

- **spaMM `HLfit` models.**
  [`spaMM::DoF()`](https://rdrr.io/pkg/spaMM/man/DoF.html) supplies the
  full count, including estimated residual, covariance, and correlation
  parameters. Realized random effects are not counted. The native
  default log-likelihood target is kept (`p_v` for ordinary maximum
  likelihood, `p_bv` for REML), so absolute REML criteria can differ
  from spaMM's marginal AIC, which uses `p_v`.

- **mgcv negative binomial BAMs.** The table checks the native fit
  against the normalized count likelihood. A retained outer scale flag
  can otherwise add one to both the reported log-likelihood and the
  parameter count, even though the family scale is one. Verified fits
  use the normalized likelihood and native effective degrees of freedom
  plus an estimated theta, with no extra scale parameter. A fixed theta
  adds no parameter. The `likelihood_adjustments` attribute records the
  native and used likelihoods and regression-model counts for the
  affected models, before the scale and shape parameters are added. With
  large case weights, verification checks the native AIC arithmetic and
  a rescaled response-density comparison at the original tolerance, and
  records that method. Reported likelihoods keep the original weights.
  Unknown inconsistencies raise an error, and fitted coefficients,
  covariance, and diagnostics are unchanged.

## Check that models are comparable

Information criteria compare models only when they describe the same
observations with comparable likelihoods. The table checks the number of
observations and, when native metadata can be recovered, the used row
IDs, response values, binomial trials, likelihood weights, and the
likelihood's response measure. Zero-weight GLM rows and removed
`unmarked` sites are excluded from the used observations. Binomial
counts and proportions with their trial weights are compared on the same
count scale, and logical `FALSE`/`TRUE` outcomes, including success and
failure columns, compare as 0 and 1 when trial totals and weights agree.
Fitted responses, covariance, and parameter counts are not changed.

- **N-mixture models.** Mixed
  [`unmarked::pcount()`](https://ecoverseR.github.io/unmarked/reference/pcount.html)
  models keep their native marginal likelihood. Two models can share
  sites and the original count matrix yet use different visits when a
  detection covariate has missing values. The table recovers the visits
  each model used from its component formulas, including
  formula-specific missing observations and omitted sites, and stops
  when they differ. Equal count matrices and site totals alone do not
  establish comparability, and missing or inconsistent used-response
  metadata raises an error without falling back on raw data. Null
  abundance models and detection-only models without site covariates
  keep the native site count and their original response-row identities.

- **Distance-removal models.** `gdistremoval()` comparisons check both
  distance and removal histories and their stored survey protocol,
  including the primary-period layout, distance bins and units, and
  removal durations. Equal paired totals are not enough. Distances in
  meters and kilometers are converted to meters in the comparison copy,
  and the fitted objects keep their native units. Harmonize time units
  and event definitions before comparing, because a native frame does
  not store a removal-time unit. Mixture and distance-key alternatives
  may share the same observed target. Check sensitivity to the
  integration cutoff separately.

- **Maximum likelihood and REML.** ML and REML fits cannot be mixed, and
  REML comparisons require the same recoverable fixed-effect design.

- **Cox and conditional logistic models.** Partial-likelihood
  comparisons require matching methods and equivalent stratum
  partitions.

These checks cannot certify metadata that a native provider does not
store. You still need to establish that the likelihoods are
scientifically comparable.

## Sample size for AICc and BIC

- **`unmarked` models.** AICc and BIC use the retained site rows, not
  visits or random-effect realizations. When groups repeat across rows,
  this convention does not establish a valid small-sample correction or
  an effective independent sample size.

- **Cox and conditional logistic models.** AICc and BIC use the event
  count from the native log-likelihood metadata, so each supported
  choice set with one used step contributes one event, while every
  alternative row is still checked for comparability. This follows the
  native survival BIC convention. It does not adjust for dependence
  among steps or animals or establish that AICc is valid for a
  particular sampling design. Missing event-count metadata raises an
  error.

A finite criterion and an optimizer code of zero do not establish
convergence, identifiability, or interval calibration. Read the native
fit diagnostics, including the gradient and curvature when available.

## Models the table rejects

- **Bayesian fits** from
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md).
  A posterior mode or penalized target is not a maximized sampling
  likelihood for AIC or BIC. Posterior predictive comparison needs a
  separate workflow with its own target and diagnostics.

- **Cross-validation results**, which hold held-out scores rather than a
  maximized likelihood. Use
  [`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
  or
  [`select_cv()`](https://wpeterman.github.io/multiScaleR/reference/select_cv.md).

- **Private spatial prototypes**, which have no qualified maximized
  marginal likelihood. Conditional response densities and posterior
  objectives are not substitutes.

- **`gdistremoval()` models with active formula offsets or TMB zero
  inflation.** These are rejected before fitting, model selection, or
  prediction. In the audited `unmarked` 1.5.2 provider, offsets did not
  enter the fitted likelihood, and the TMB ZIP likelihood ignored its
  zero-inflation parameter. These combinations need separate provider
  qualification. multiScaleR does not remove offsets, switch engines,
  alter parameter counts, or substitute uncertainty estimates.
  Offset-free C-engine ZIP models keep their native likelihood and
  parameter count.

The table itself is built with
[`AICcmodavg::aictabCustom()`](https://rdrr.io/pkg/AICcmodavg/man/aictabCustom.html).

## See also

[`bic_tab()`](https://wpeterman.github.io/multiScaleR/reference/bic_tab.md),
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`cv_tab()`](https://wpeterman.github.io/multiScaleR/reference/cv_tab.md)
for comparisons by held-out prediction, and
[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md).

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
models <- list(spatial = spatial_model, fixed = fixed_model, null = null_model)

## The fixed model already uses the true 75-unit scale, so the spatial model
## reaches the same log-likelihood and pays for one extra parameter (sigma).
aic_tab(models, AICc = FALSE, mod_names = names(models))
#> 
#> Model selection based on AIC:
#> 
#>         K    AIC Delta_AIC AICWt Cum.Wt      LL
#> fixed   2 232.34      0.00  0.73   0.73 -114.17
#> spatial 3 234.33      1.99  0.27   1.00 -114.17
#> null    1 346.69    114.35  0.00   1.00 -172.34
#> 
aic_tab(models, AICc = TRUE, mod_names = names(models))
#> 
#> Model selection based on AICc:
#> 
#>         K   AICc Delta_AICc AICcWt Cum.Wt      LL
#> fixed   2 232.51       0.00   0.75   0.75 -114.17
#> spatial 3 234.67       2.16   0.25   1.00 -114.17
#> null    1 346.74     114.23   0.00   1.00 -172.34
#> 

## K for a stored N-mixture fit: three abundance coefficients, one
## detection coefficient, and two estimated scales
nmix <- readRDS(system.file("extdata", "cache-user-guide-nmix.rds",
                            package = "multiScaleR"))
aic_tab(list(nmix), mod_names = "bin1 + cont2")
#> 
#> Model selection based on AICc:
#> 
#>              K    AICc Delta_AICc AICcWt Cum.Wt      LL
#> bin1 + cont2 6 1001.89          0      1      1 -494.33
#> 
```
