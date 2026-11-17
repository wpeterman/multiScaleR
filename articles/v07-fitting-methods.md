# 07. Choose a Fitting Method for Your Model

## What this vignette covers

`multiScaleR` can estimate a scale of effect in several ways. A *fitting
route* is the combination of function and settings that does the
estimation, either a `backend` of
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
or a `method` of
[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md).
The routes differ in which models they accept, which kinds of landscape
covariate and formula term they allow, and which intervals they can
report. This vignette answers one question. Which route should you use
for your model, and what does it support?

Two tables give the answer, and four short examples on the Quick Start
data show what changes when you fit the same model by a different route.
Read
[`vignette("v01-quickstart", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v01-quickstart.md)
first. What each interval includes and leaves out is the subject of
[`vignette("v09-estimates-and-intervals", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v09-estimates-and-intervals.md),
which this vignette does not repeat.

## Choose a route for your model

Find your model in the first column of the table below. Every route
listed beside it accepts that model, although a route can still reject a
particular term or covariate, which the next section covers.

| Your starting model | Routes that fit it |
|:---|:---|
| Gaussian (identity link), Poisson (log), or binomial (logit) model fitted with [`glm()`](https://rdrr.io/r/stats/glm.html) | `refit`, `profile_glm`, [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md) |
| Negative binomial model fitted with [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) | `refit`, `profile_glm`, `joint_nb2`, [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md) |
| [`lm()`](https://rdrr.io/r/stats/lm.html) model, or a GLM with another family or link such as Gamma with a log link | `refit` |
| Generalized additive model fitted with [`mgcv::gam()`](https://rdrr.io/pkg/mgcv/man/gam.html) | `refit` |
| Mixed model fitted with `lme4` or `glmmTMB` | `refit`. A Poisson or negative binomial `glmmTMB` model with one random intercept, optionally zero-inflated, can also use `joint_mixed`, and [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md) accepts random intercepts |
| Step-selection model fitted with [`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html), one used step per choice set | `refit`, `conditional_logit` |
| Occupancy or abundance model fitted with `unmarked` | `refit` |
| Any count, binary, or Gaussian response whose scale of effect varies among years or depends on a covariate | [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md) |

The names in the second column are values of the `backend` argument of
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
except for
[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
which is a separate function. Three rules of thumb cover most analyses.

- **Start with `refit`.** It is the default. It refits your model at
  each candidate scale, so it works with every model class in the table
  and every covariate type.
- **Move to `profile_glm` or `joint_nb2` when a prediction interval must
  carry scale uncertainty.** These backends estimate the scales and the
  coefficients in one likelihood, which makes a joint interval possible.
- **Use
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  when the scale itself is the object of study.** It is the only route
  that lets a scale differ among years or follow a covariate.

[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
takes a formula in place of a starting model, and its `method` argument
selects how the posterior is computed. Apart from the requirement noted
for `"marginal"`, every method accepts the same model structures.

| `method` | What it does | Optional packages it needs |
|:---|:---|:---|
| `"sample"` (default) | Samples the full posterior with the NUTS algorithm | `RTMB` and `tmbstan` |
| `"slice"` | Samples the same posterior with a slice sampler written in R | None. Install `posterior` to compute the sampling diagnostics |
| `"laplace"` | Fits a Gaussian approximation to the posterior | `RTMB` |
| `"marginal"` | Integrates one between-group standard deviation of scale numerically. The model needs exactly one random intercept in its scale formulas | `RTMB` |
| `"prepare"` | Checks the model and returns without fitting | None |

## Check what each route supports

The next table lists the covariate types, formula terms, and outputs
that each route allows. Each entry was checked by running the route on
the Quick Start data or by reading the error it returns. `joint_mixed`
is left out because it follows the same rules as `conditional_logit` for
kernels and covariate types, and its own requirements appear under
[Current restrictions](#current-restrictions).

| Feature | `refit` | `profile_glm`, `joint_nb2` | `conditional_logit` | [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md) |
|:---|:---|:---|:---|:---|
| Kernel in [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md) | Gaussian, exponential, exponential power, fixed | Gaussian, exponential, exponential power | Gaussian, exponential | Gaussian, exponential |
| Binned inputs (`bin = TRUE`, the default) | Optional, and much faster | Required | Required | Required |
| Kernel-weighted mean, [`kernel_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md) | Yes | Yes | Yes | Yes |
| Kernel-weighted root-mean-square roughness, `surface_var(..., "sq", weighted = TRUE)` | Yes | Yes | Yes | Yes |
| Other kernel-weighted surface metrics, and any surface or landscape metric with an estimated radius | Yes | No | No | No |
| Search over a list of candidate radii, `radius_grid` | [`lm()`](https://rdrr.io/r/stats/lm.html), [`glm()`](https://rdrr.io/r/stats/glm.html), and step-selection models | No | No | No |
| Metric computed at a fixed radius | Yes | Yes | Yes | Yes |
| Interactions with numeric or factor predictors, and between landscape covariates | Yes | Yes | Yes | Yes |
| Arithmetic inside [`I()`](https://rdrr.io/r/base/AsIs.html), such as `I(land2^2)` | Yes | Yes | Yes | Yes |
| [`poly()`](https://rdrr.io/r/stats/poly.html) of a landscape covariate | Yes | No | No | No |
| [`splines::ns()`](https://rdrr.io/r/splines/ns.html) or [`splines::bs()`](https://rdrr.io/r/splines/bs.html) with fixed `knots` and `Boundary.knots` | Yes | Yes | Yes | Yes |
| [`splines::ns()`](https://rdrr.io/r/splines/ns.html) or [`splines::bs()`](https://rdrr.io/r/splines/bs.html) with `df` | Yes | No | Use `refit` | No |
| Smooth terms such as `s()` in `mgcv` | Yes | No | No | No |
| Random effects | Whatever the model class allows | No | No | Random intercepts, `(1 | group)` |
| Offsets | Yes | [`offset()`](https://rdrr.io/r/stats/offset.html) in the formula | [`offset()`](https://rdrr.io/r/stats/offset.html) in the formula | [`offset()`](https://rdrr.io/r/stats/offset.html) in the formula |
| Case weights | Whatever the model class allows | Whole-number frequency weights for count and binomial models | No | No |
| Covariates kept in their own units, `covariate_transform = "raw"` | Yes | No | No | No |
| Interval for sigma | Wald, or profile likelihood with `summary(fit, profile = TRUE)` | Wald, or profile likelihood with `summary(fit, profile = TRUE)` | None | Posterior credible interval |
| Uncertainty from [`predict()`](https://rspatial.github.io/terra/reference/predict.html) | Conditional standard errors | Conditional or joint standard errors and confidence intervals | Conditional standard errors | Posterior intervals |
| Model comparison with [`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md) and [`bic_tab()`](https://wpeterman.github.io/multiScaleR/reference/bic_tab.md) | Yes | Yes | Yes | No |

**How to read the table:** A *conditional* standard error treats the
estimated scales as known. A *joint* one lets the scales and the
coefficients vary together. A landscape covariate here means a column
that
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
builds from a raster layer. Four rows need a word of explanation.

- **Covariate types.**
  [`kernel_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  is the distance-weighted mean used in the Quick Start.
  [`surface_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  summarizes how much a continuous layer varies around a survey point,
  and
  [`landscape_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
  summarizes the composition or arrangement of land-cover classes. See
  [`vignette("v05-surface-metrics", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v05-surface-metrics.md)
  and
  [`vignette("v04-landscape-metrics", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v04-landscape-metrics.md).
- **Fixed-radius metrics.** A metric with `radius` set has no scale to
  estimate.
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  computes it once and every route treats it as an ordinary numeric
  predictor.
- **A grid of radii.** With `radius_grid`, an experimental option,
  `refit` fits the model at every combination of the radii you list and
  re-estimates any kernel sigma at each.
  `diagnostics(fit)$radius_profile` keeps the likelihood of every
  combination. An interval for a kernel sigma then treats the chosen
  radii as known.
  [`vignette("v04-landscape-metrics", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v04-landscape-metrics.md)
  works an example.
- **Raw units.** By default each landscape covariate is centered and
  scaled again at every candidate scale, so a coefficient describes a
  change of one standard deviation. With `covariate_transform = "raw"`
  the covariate keeps its own units, which for the binary `land1` layer
  is the kernel-weighted proportion of habitat.
  [`?multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  describes how to set up the starting model.

## Fit one model by two routes

### Prepare the data once

Every route reads the same
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
object, so you prepare the inputs once. The data are the counts and
landscape layers from the Quick Start, where the true scales of effect
are 250 m for `land1` and 500 m for `land2`.

``` r

data("landscape_counts")
data("surv_pts")
pts <- vect(surv_pts)
land_rast <- rast(system.file("extdata", "landscape.tif",
                              package = "multiScaleR"))

# bin = TRUE (the default) stores distance-binned summaries for every route
inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                      kernel = "gaussian", bin = TRUE, verbose = FALSE)

# Starting model: its formula and family carry through to the fit
df    <- data.frame(landscape_counts, inputs$kernel_dat)
start <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)
```

### Compare the scale estimates

``` r

# Default backend: refit the GLM at every candidate scale
time_refit <- system.time(
  fit_refit <- multiScale_optim(start, inputs, verbose = FALSE)
)

# Analytic backend: scales and coefficients in one likelihood
time_profile <- system.time(
  fit_profile <- multiScale_optim(start, inputs, backend = "profile_glm",
                                  n_cores = 1, verbose = FALSE)
)

# Sigma, its standard error, and the default 95% interval, in meters
est_refit   <- summary(fit_refit)$opt_scale
est_profile <- summary(fit_profile)$opt_scale
scale_table <- data.frame(
  backend   = rep(c("refit", "profile_glm"), each = 2),
  covariate = c(rownames(est_refit), rownames(est_profile)),
  rbind(est_refit, est_profile),
  check.names = FALSE
)
knitr::kable(scale_table, digits = 1, row.names = FALSE,
             col.names = c("Backend", "Covariate", "Sigma (m)", "SE",
                           "Lower", "Upper"))
```

| Backend     | Covariate | Sigma (m) |   SE | Lower | Upper |
|:------------|:----------|----------:|-----:|------:|------:|
| refit       | land1     |     238.5 | 33.1 | 172.7 | 304.2 |
| refit       | land2     |     499.6 | 47.5 | 405.3 | 593.8 |
| profile_glm | land1     |     238.5 | 33.1 | 181.7 | 313.1 |
| profile_glm | land2     |     499.6 | 47.5 | 414.7 | 601.8 |

``` r


# Seconds each fit took
round(c(refit = time_refit[["elapsed"]],
        profile_glm = time_profile[["elapsed"]]), 2)
#>       refit profile_glm 
#>        0.21        0.18
```

**How to read the output:** Both backends maximize the same likelihood
for this model, so the sigma estimates and standard errors agree to the
first decimal place. For `land1` both report 238 m, and for `land2` both
report 500 m. The intervals differ slightly. `refit` builds a symmetric
interval around sigma. `profile_glm` builds its interval on the log of
sigma, so the upper limit sits farther from the estimate than the lower
limit does. Both intervals for each covariate contain the true value.

Both fits finish quickly, because both read binned inputs. Binning is
what matters for speed. The same `refit` call on inputs prepared with
`bin = FALSE` took more than 20 seconds when this vignette was written.

### Compare the prediction uncertainty

``` r

# refit: standard errors that hold the estimated scales fixed
se_refit <- predict(fit_refit, se.fit = TRUE, uncertainty = "conditional")
round(head(data.frame(expected = se_refit$fit, se = se_refit$se.fit), 3), 2)
#>   expected   se
#> 1     0.91 0.14
#> 2     2.58 0.45
#> 3     4.69 0.47
```

``` r

# refit cannot return a joint interval
predict(fit_refit, interval = "confidence", uncertainty = "joint")
#> Error:
#> ! Refit predictions support native conditional SEs only; use uncertainty = 'conditional' and interval = 'none'.
```

``` r

# profile_glm: scales and coefficients vary together
joint <- predict(fit_profile, interval = "confidence", uncertainty = "joint")
round(head(data.frame(expected = joint$fit, se = joint$se.fit,
                      lower = joint$lower, upper = joint$upper), 3), 2)
#>   expected   se lower upper
#> 1     0.91 0.14  0.67  1.23
#> 2     2.58 0.60  1.64  4.07
#> 3     4.69 0.64  3.59  6.12
```

**How to read the output:** Each row is a survey point, and `expected`
is the expected count per survey. The two fits predict the same counts
and differ only in the uncertainty attached to them. At the second
survey point the conditional standard error is 0.45 and the joint
standard error is 0.6, because the joint value also carries the
uncertainty in the two scales. The default backend stops with an error
when you ask it for an interval, and its message names the one setting
it supports.

> **What this means for your study.** Choose the backend before you need
> the interval. If a report or a map must show uncertainty that includes
> the scale of effect, fit with `profile_glm` or `joint_nb2` from the
> start, or use
> [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md).

## Add a curved response

A spline lets the response bend along a landscape covariate. Two recipes
are common. One names the knots, or the covariate values where the curve
is allowed to change shape. The other gives only `df`, the number of
basis columns, and lets R place the knots at quantiles of the data.

``` r

# Recipe 1: knots and boundaries written out
start_knots <- glm(counts ~ site + land1 +
                     splines::ns(land2, knots = c(-0.5, 0.5),
                                 Boundary.knots = c(-4, 4)),
                   family = poisson(), data = df)

# Recipe 2: R chooses the knots from the data
start_df <- glm(counts ~ site + land1 + splines::ns(land2, df = 3),
                family = poisson(), data = df)

# Written-out knots work with the analytic backend
fit_knots <- multiScale_optim(start_knots, inputs, backend = "profile_glm",
                              n_cores = 1, verbose = FALSE)
```

``` r

# The df recipe does not
multiScale_optim(start_df, inputs, backend = "profile_glm",
                 n_cores = 1, verbose = FALSE)
#> Error:
#> ! Dynamic splines require explicit fixed knots and Boundary.knots, without adaptive df; use the refit backend otherwise.
```

``` r

# The default backend accepts either recipe
fit_df <- multiScale_optim(start_df, inputs, verbose = FALSE)

# Compare the straight-line model with the two curved ones
aic_tab(list(fit_profile, fit_knots, fit_df),
        mod_names = c("linear", "spline, fixed knots", "spline, df = 3"))
#> 
#> Model selection based on AICc:
#> 
#>                     K   AICc Delta_AICc AICcWt Cum.Wt      LL
#> linear              6 319.15       0.00   0.62   0.62 -153.12
#> spline, df = 3      8 321.36       2.21   0.20   0.82 -151.89
#> spline, fixed knots 8 321.64       2.49   0.18   1.00 -152.03
```

**How to read the output:** The quantile knots of the `df` recipe move
each time the optimizer tries a new scale, because the covariate values
change with the scale. `profile_glm` needs a basis that stays put, so it
stops and points you to `refit`. In the table, `K` counts the regression
coefficients plus the two scales, and a lower `AICc` means better
support. The straight-line model ranks first. Each spline spends two
more parameters and gains about one unit of log-likelihood (`LL`), which
matches how the data were simulated, with a straight-line effect of
`land2` on the log of the expected count. The two spline fits differ
little from each other.

Knots are in the units of the modeled covariate. Landscape covariates
are standardized by default, so the knot at 0.5 sits half a standard
deviation above the mean kernel-weighted value of `land2`. It is not a
fixed amount of `land2`. All spline columns share the one scale of
effect estimated for `land2`.

## Mix covariate types in one model

[`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md)
names each covariate and says how to compute it. The model below adds
the roughness of `land3` within a fixed 300 m radius to the two
kernel-weighted means. The `"sa"` metric is the mean absolute deviation
of the layer from its local mean.

``` r

# Name each covariate and how to summarize its source layer
features <- msr_vars(
  land1  = kernel_var("land1"),                      # mean, scale estimated
  land2  = kernel_var("land2"),                      # mean, scale estimated
  rough3 = surface_var("land3", "sa", radius = 300)  # roughness within 300 m
)

inputs_mix <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                          kernel = "gaussian", scale_vars = features,
                          bin = TRUE, verbose = FALSE)

df_mix    <- data.frame(landscape_counts, inputs_mix$kernel_dat)
start_mix <- glm(counts ~ site + land1 + land2 + rough3,
                 family = poisson(), data = df_mix)

# Either backend accepts the fixed-radius metric
fit_mix_refit   <- multiScale_optim(start_mix, inputs_mix, verbose = FALSE)
fit_mix_profile <- multiScale_optim(start_mix, inputs_mix,
                                    backend = "profile_glm",
                                    n_cores = 1, verbose = FALSE)

round(fit_mix_profile$scale_est, 1)                  # estimated scales (m)
#>        Mean   SE
#> land1 234.3 34.2
#> land2 493.9 48.5
round(coef(summary(fit_mix_profile$opt_mod)), 3)     # regression table
#>             Estimate Std. Error z value Pr(>|z|)
#> (Intercept)    0.391      0.091   4.278    0.000
#> site           0.172      0.083   2.063    0.039
#> land1         -0.477      0.080  -5.954    0.000
#> land2          0.630      0.074   8.480    0.000
#> rough3         0.035      0.074   0.477    0.633
```

**How to read the output:** The scale table has two rows, one for each
kernel-weighted mean. `rough3` has a regression coefficient and no
scale, because its radius was fixed in advance. Its coefficient is 0.035
with a standard error of 0.074, which is consistent with no effect. That
is the right answer, since `land3` played no part in simulating the
counts. `fit_mix_refit` returns the same estimates.

To predict at new points, prepare them with the same `features` and call
[`predict()`](https://rspatial.github.io/terra/reference/predict.html)
with `include_fixed = TRUE`. The fixed-radius column is then rebuilt and
standardized with the mean and standard deviation of the fitted data.

Had `rough3` been declared without `radius`, its radius would be
estimated, and only `refit` could fit the model. The same holds for a
kernel-weighted metric other than `"sq"`. Fits that estimate a
kernel-weighted roughness also take longer than fits of kernel-weighted
means, on the order of 5 to 40 seconds for these data against less than
1 second.

## Check a Bayesian model before you sample

Sampling a posterior can take minutes, so it pays to learn first whether
[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
accepts your model. `method = "prepare"` runs every check and returns
without fitting. The function rebuilds the landscape covariates from
`inputs`, so `data` holds only the response and the covariates measured
at the survey point.

``` r

# Validate the formula, rows, priors, and bins. Nothing is sampled.
bayes_check <- multiScale_bayes(counts ~ site + land1 + land2,
                                inputs = inputs, data = landscape_counts,
                                family = "poisson", method = "prepare")
bayes_check
#> Experimental Bayesian landscape-scale model
#> Inference: prepared 
#> Family: poisson | observations: 100 
#> Scale formulas:
#>   land1 : ~1 
#>   land2 : ~1
```

``` r

# A term the Bayesian route rejects fails at once
multiScale_bayes(counts ~ site + land1 + splines::ns(land2, df = 3),
                 inputs = inputs, data = landscape_counts,
                 family = "poisson", method = "prepare")
#> Error:
#> ! Dynamic splines require explicit fixed knots and Boundary.knots, without adaptive df; use the refit backend otherwise.
```

**How to read the output:** The first call succeeds. Its printout lists
the family, the number of observations, and one scale formula per
landscape covariate, where `~1` means a single shared scale. The second
call stops with the same message `profile_glm` gave, because the
Bayesian route also needs a spline basis with written-out knots. Once a
model passes, change `method` to fit it.
[`vignette("v13-bayesian-scales", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v13-bayesian-scales.md)
works through a full fit.

A Bayesian fit standardizes each landscape covariate once, with a fixed
center and spread, while
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
standardizes again at each candidate scale. Regression coefficients from
the two are therefore not on the same scale. Compare their scales of
effect and their predictions.

## Current restrictions

| Route or feature | Restriction |
|:---|:---|
| `profile_glm` | Experimental. In simulations its 95% Wald intervals for sigma held the true value in 76% to 92% of data sets, depending on the response family. Profile-likelihood limits did better |
| `profile_glm`, `joint_nb2` | Need a [`glm()`](https://rdrr.io/r/stats/glm.html) or [`MASS::glm.nb()`](https://rdrr.io/pkg/MASS/man/glm.nb.html) object, and `joint_nb2` needs the negative binomial one. Set `n_cores = 1`. Put offsets in the formula with [`offset()`](https://rdrr.io/r/stats/offset.html) |
| Joint intervals | Withheld when a fit fails its numerical checks, for example when a scale sits at a search bound or when a Poisson model fits negative binomial data equally well. `diagnostics(fit)$profile_glm$reasons` gives the reason, and conditional standard errors remain available |
| `conditional_logit` | Reports sigma with no standard error and no interval. Needs complete choice sets with one used step each, no weights, and `n_cores = 1` |
| [`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html) and `coxph()` starting models | When the formula uses [`scale()`](https://rspatial.github.io/terra/reference/scale.html), [`poly()`](https://rdrr.io/r/stats/poly.html), or a spline, fit with `model = TRUE` and `na.action = stats::na.fail`. Recent versions of `survival` otherwise drop the stored center, spread, or knots, and predictions cannot recover the original units |
| `joint_mixed` | Needs the optional packages `RTMB`, `TMB`, and `glmmTMB`, a Gaussian or exponential kernel, and a `glmmTMB` model. Reports sigma with no standard error, and conditional prediction standard errors only |
| All backends except `refit` | Match model rows to survey points by row name. Build the model data from `kernel_dat` and keep its row names when you subset |
| Estimated radius of a hard-edged metric, and `kernel = "fixed"` | Sigma has no standard error. Use `summary(fit, profile = TRUE)` for an interval. With `kernel = "fixed"`, [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md) screens starting values by default, because a single start cannot move. The estimate is then the best screened candidate |
| [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md) | Experimental, and printed as such. [`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md) and [`bic_tab()`](https://wpeterman.github.io/multiScaleR/reference/bic_tab.md) stop with an error for a Bayesian fit |

## Quick reference

``` r

library(multiScaleR)
library(terra)

# 1. Prepare once. Binned inputs serve every route.
data("landscape_counts"); data("surv_pts")
pts       <- vect(surv_pts)
land_rast <- rast(system.file("extdata", "landscape.tif",
                              package = "multiScaleR"))
inputs <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                      kernel = "gaussian", bin = TRUE)
df     <- data.frame(landscape_counts, inputs$kernel_dat)
start  <- glm(counts ~ site + land1 + land2, family = poisson(), data = df)

# 2. Same model, two routes
fit_refit   <- multiScale_optim(start, inputs)               # any model class
fit_profile <- multiScale_optim(start, inputs, backend = "profile_glm",
                                n_cores = 1)                 # joint intervals
summary(fit_refit)$opt_scale       # sigma, SE, symmetric 95% interval (m)
summary(fit_profile)$opt_scale     # sigma, SE, 95% interval on the log scale

# 3. Prediction uncertainty by route
predict(fit_refit, se.fit = TRUE, uncertainty = "conditional")       # scales fixed
predict(fit_profile, interval = "confidence", uncertainty = "joint") # scales uncertain

# 4. Splines: written-out knots for profile_glm, either recipe for refit
start_knots <- glm(counts ~ site + land1 +
                     splines::ns(land2, knots = c(-0.5, 0.5),
                                 Boundary.knots = c(-4, 4)),
                   family = poisson(), data = df)
start_df    <- glm(counts ~ site + land1 + splines::ns(land2, df = 3),
                   family = poisson(), data = df)
fit_knots <- multiScale_optim(start_knots, inputs, backend = "profile_glm",
                              n_cores = 1)
fit_df    <- multiScale_optim(start_df, inputs)
aic_tab(list(fit_profile, fit_knots, fit_df),
        mod_names = c("linear", "spline, fixed knots", "spline, df = 3"))

# 5. Mixed covariate types: two estimated scales and one fixed-radius metric
features <- msr_vars(
  land1  = kernel_var("land1"),
  land2  = kernel_var("land2"),
  rough3 = surface_var("land3", "sa", radius = 300)
)
inputs_mix <- kernel_prep(pts = pts, raster_stack = land_rast, max_D = 1700,
                          kernel = "gaussian", scale_vars = features,
                          bin = TRUE)
df_mix    <- data.frame(landscape_counts, inputs_mix$kernel_dat)
start_mix <- glm(counts ~ site + land1 + land2 + rough3,
                 family = poisson(), data = df_mix)
fit_mix_refit   <- multiScale_optim(start_mix, inputs_mix)
fit_mix_profile <- multiScale_optim(start_mix, inputs_mix,
                                    backend = "profile_glm", n_cores = 1)

# 6. Check a Bayesian model without fitting it
bayes_check <- multiScale_bayes(counts ~ site + land1 + land2,
                                inputs = inputs, data = landscape_counts,
                                family = "poisson", method = "prepare")
```

### Summary of key functions

| Function | Purpose |
|:---|:---|
| [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md) | Prepare the inputs every route reads. `bin = TRUE` is required by all routes except `refit` |
| [`msr_vars()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md), [`kernel_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md), [`surface_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md), [`landscape_var()`](https://wpeterman.github.io/multiScaleR/reference/msr_vars.md) | Declare kernel-weighted means, surface metrics, and land-cover metrics, each with an estimated or a fixed scale |
| [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md) | Estimate scales by maximum likelihood. `backend` selects the route |
| [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md) | Estimate a posterior for scales and coefficients. `method = "prepare"` checks a model without fitting |
| [`summary()`](https://rspatial.github.io/terra/reference/summary.html) | Report sigma with a Wald or profile-likelihood interval |
| [`predict()`](https://rspatial.github.io/terra/reference/predict.html) | Expected response with the uncertainty the route supports |
| [`diagnostics()`](https://wpeterman.github.io/multiScaleR/reference/diagnostics.md) | Report whether a joint interval is available, and why not |
| [`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md), [`bic_tab()`](https://wpeterman.github.io/multiScaleR/reference/bic_tab.md) | Compare likelihood fits from any [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md) backend |

### See also

- **[`vignette("v01-quickstart", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v01-quickstart.md)**,
  the core workflow with the default backend.
- **[`vignette("v09-estimates-and-intervals", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v09-estimates-and-intervals.md)**,
  what a conditional, joint, or posterior interval includes and leaves
  out.
- **[`vignette("v05-surface-metrics", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v05-surface-metrics.md)**,
  surface roughness and texture metrics.
- **[`vignette("v04-landscape-metrics", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v04-landscape-metrics.md)**,
  land-cover composition and configuration metrics.
- **[`vignette("v10-step-selection", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v10-step-selection.md)**,
  step-selection models with `backend = "conditional_logit"`.
- **[`vignette("v06-year-specific-rasters", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v06-year-specific-rasters.md)**,
  surveys from several years, each with its own map.
- **[`vignette("v13-bayesian-scales", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v13-bayesian-scales.md)**,
  fitting and checking a Bayesian model whose scales vary among years.
- **[`vignette("v08-cross-validation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v08-cross-validation.md)**,
  comparing models by how well they predict held-out observations.
- The help page
  [`?"multiScaleR-glossary"`](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-glossary.md)
  defines the terms used in these vignettes.
