# Estimate scales of effect jointly across independent data sets

**Experimental.** `multiScale_joint()` combines two or more data sets,
each already fitted with
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
and estimates the scale of effect of each landscape covariate jointly.
Each data set keeps its own response, model class, and native parameters
(for example a Poisson count model and an `unmarked` occupancy model).
What is shared is set by `scales` and, if you choose, `share`:

- `scales = "common"`: one scale for each landscape covariate, used by
  every data set that contains it.

- `scales = "separate"`: one scale per data set. This reproduces the
  separate fits and is the reference for testing whether a common scale
  is supported.

- `scales = "footprint"`: one ecological scale, seen through each data
  set's sampling footprint (see Details). The apparent scale in data set
  `d` is `sqrt(sigma^2 + footprint[d]^2)`.

- `scales = "footprint_scaled"`: as `"footprint"`, with the footprint
  SDs multiplied by an estimated factor `kappa`. `kappa = 0` is the
  common scale and `kappa = 1` is the stated footprint.

`share` names landscape covariates whose regression coefficient is also
shared across data sets, so the data sets estimate one landscape effect
but keep their own intercepts and other coefficients.

## Usage

``` r
multiScale_joint(
  fits,
  scales = "common",
  footprint = NULL,
  share = NULL,
  cluster = NULL,
  uncertainty = c("hessian", "jackknife"),
  profile = TRUE,
  level = 0.95,
  groups = 20L,
  seed = NULL,
  verbose = TRUE
)
```

## Arguments

- fits:

  Named list of two or more fitted `multiScaleR` objects from
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
  with the default `backend = "refit"`, one per data set. The names
  identify the data sets in the output. Each fit supplies its starting
  model, its prepared kernel inputs, and starting values. All fits must
  use the same kernel, Gaussian or exponential. Landscape covariates are
  matched across data sets by name, so give the same raster layer the
  same name in every
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  call.

- scales:

  How the scale of each landscape covariate relates across data sets:
  one of `"common"` (default), `"separate"`, `"footprint"`, or
  `"footprint_scaled"` for every covariate, or a named character vector
  with one structure per covariate, such as
  `c(forest = "common", road = "separate")`. A covariate that occurs in
  only one data set always has its own scale.

- footprint:

  Named numeric vector with one sampling footprint SD per data set, in
  map units, required by the footprint structures. Use 0 for a data set
  whose response is measured at a point.
  [`footprint_sd()`](https://wpeterman.github.io/multiScaleR/reference/footprint_sd.md)
  converts common sampling geometries to this SD.

- share:

  Optional character vector of landscape covariates whose coefficient is
  shared across the data sets that contain them. Each such covariate
  must enter every model as one plain linear term, and the models must
  use comparable links (all the same link, or a mix of log and
  complementary log-log, which share the log-intensity scale). Shared
  coefficients currently require models refitted by
  [`update()`](https://rspatial.github.io/terra/reference/update.html):
  `glm`, [`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html),
  lme4, glmmTMB, and mgcv. Default `NULL` shares no coefficient.

- cluster:

  Optional site or block identifiers, used to check for and allow
  dependence between data sets that sample the same places. Either the
  name of a column present in every model's data, or a named list with
  one vector per data set aligned with that model's data rows (or
  `unmarked` sites). Default `NULL`.

- uncertainty:

  `"hessian"` (default) takes standard errors from the Hessian of the
  joint log-likelihood, which assumes the data sets are independent.
  `"jackknife"` deletes one group of clusters at a time and refits,
  which allows dependence within clusters, including between data sets
  that share a site. Its Wald limits use a t distribution with one fewer
  degree of freedom than the number of replicates. It requires `cluster`
  and is not yet available for `unmarked` models.

- profile:

  Logical. If `TRUE` (default), compute profile-likelihood intervals for
  every estimated parameter.

- level:

  Confidence level for the intervals. Default `0.95`.

- groups:

  Number of jackknife groups when `uncertainty = "jackknife"`. Clusters
  are assigned to groups at random. Default `20`, or the number of
  clusters if smaller.

- seed:

  Optional integer seed for the jackknife group assignment.

- verbose:

  Logical. If `TRUE` (default), report progress.

## Value

An object of class `multiScaleR_joint`, a list with

- `parameters`: estimated parameters (scales in map units, `kappa`,
  shared coefficients) with standard errors and Wald and
  profile-likelihood limits.

- `scales`: the apparent scale of each covariate in each data set.

- `fits`: named list of component `multiScaleR` fits.

- `comparison`: log-likelihoods, parameter counts, and the
  likelihood-ratio test against separate scales.

- `residual_correlation`: correlation of residuals between data sets at
  shared clusters, when `cluster` is supplied.

- `logLik`, `df`, `nobs`, `vcov` (parameters on the optimizer scale: log
  scales, `kappa`, coefficients), `structure`, `footprint`, `share`,
  `standardization`, `uncertainty`, `optim`, and `call`.

## Details

**The joint likelihood.** For a given set of scales, each data set's
model is refitted exactly as
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
refits it, and its maximized log-likelihood is returned.
`multiScale_joint()` maximizes the sum of these log-likelihoods over the
scale parameters (and shared coefficients). With `scales = "separate"`
and no shared coefficient the sum separates, and the result equals the
separate fits. Because the joint fit treats the data sets as independent
pieces of evidence, it is valid when they are independent given the
covariates. Data sets that sample the same sites may share unexplained
site effects; supply `cluster` to measure the residual correlation at
shared sites, and use `uncertainty = "jackknife"` when it is not
negligible.

**Comparing structures.** The returned object reports the
likelihood-ratio test of the fitted structure against separate scales.
With two or three data sets this test has little power, so a
nonsignificant result is weak evidence for a common scale. Fit several
structures and compare them with
[`stats::AIC()`](https://rdrr.io/r/stats/AIC.html) (the fits have
[`stats::logLik()`](https://rdrr.io/r/stats/logLik.html) methods) or
[`stats::anova()`](https://rdrr.io/r/stats/anova.html). To test
`kappa = 0`, compare a `"common"` fit with a `"footprint_scaled"` fit in
[`stats::anova()`](https://rdrr.io/r/stats/anova.html). That value lies
on the boundary of the parameter space, so
[`anova()`](https://rdrr.io/r/stats/anova.html) refers the statistic to
a mixture of chi-square distributions, which halves the p-value for one
`kappa`. With a single data set that has a footprint, `kappa` can match
that data set's apparent scale exactly, so the ecological scale then
rests on the point data sets alone.

**Sampling footprints.** A response is often measured over an area, not
at a point: counts within a fixed radius, detections by an acoustic
recorder, a camera viewshed, or animals detected by distance sampling.
If the ecological response follows the landscape smoothed with a
Gaussian kernel of scale `sigma`, then a response averaged over a
footprint with per-axis SD `tau` follows the landscape smoothed by the
convolution of the two, whose per-axis variance is `sigma^2 + tau^2`.
Data sets with larger footprints therefore show larger apparent scales
even when the ecological scale is the same. The relation is exact for a
Gaussian footprint and for a covariate that enters the linear predictor
linearly; for a disk or square footprint it matches the second moment,
and on a log or logit link it is a first-order approximation that holds
while the landscape effect varies little within one footprint. Footprint
structures require the Gaussian kernel.

**Shared coefficients.** For each shared covariate, every data set's
model is refitted with that term replaced by a fixed offset equal to the
shared coefficient times the covariate. The covariate is standardized
with one center and spread pooled across the data sets, evaluated at the
ecological scale, so the coefficient has the same meaning in every data
set.

**Component fits.** `$fits` holds one `multiScaleR` object per data set
at the jointly estimated scales, with standard errors from the joint
fit. Use them with
[`plot.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/plot.multiScaleR.md),
[`plot_marginal_effects()`](https://wpeterman.github.io/multiScaleR/reference/plot_marginal_effects.md),
[`kernel_scale.raster()`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md),
and
[`predict()`](https://rspatial.github.io/terra/reference/predict.html)
as for any other fit. Two cautions apply: `summary(fit, profile = TRUE)`
and
[`profile_sigma()`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md)
on a component profile that data set alone, and a shared coefficient
enters a component's model as an offset named `.msr_joint_offset`, so
effect curves and maps of that covariate are not yet available from the
component.

## See also

[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`footprint_sd()`](https://wpeterman.github.io/multiScaleR/reference/footprint_sd.md),
[`profile_sigma()`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md).
[`vignette("v16-joint-data-sets", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v16-joint-data-sets.md)
works through an example with three surveys, shared sites, and a
sampling footprint.

## Examples

``` r
# \donttest{
# Two independent surveys of one landscape: Poisson counts and binary
# detections. Both respond to `bin1` and `cont2` at the same true scales.
rs <- sim_rast(user_seed = 123, dim = 100, resolution = 30)
rs <- terra::subset(rs, c("bin1", "cont2"))

counts <- sim_dat(alpha = 0.5, beta = c(0.75, -0.75), kernel = "gaussian",
                  sigma = c(75, 150), type = "count", n_points = 60,
                  raster_stack = rs, max_D = 450, user_seed = 1)
occ <- sim_dat(alpha = 0, beta = c(1, -1), kernel = "gaussian",
               sigma = c(75, 150), type = "occ", n_points = 80,
               raster_stack = rs, max_D = 450, user_seed = 2)

in_counts <- kernel_prep(counts$pts, rs, max_D = 450, verbose = FALSE)
in_occ <- kernel_prep(occ$pts, rs, max_D = 450, verbose = FALSE)

fit_counts <- multiScale_optim(
  glm(y ~ bin1 + cont2, family = poisson,
      data = data.frame(y = counts$obs, in_counts$kernel_dat)),
  in_counts, verbose = FALSE)
fit_occ <- multiScale_optim(
  glm(y ~ bin1 + cont2, family = binomial,
      data = data.frame(y = occ$obs, in_occ$kernel_dat)),
  in_occ, verbose = FALSE)

joint <- multiScale_joint(list(counts = fit_counts, occ = fit_occ),
                          scales = "common", verbose = FALSE)
joint
#> 
#> Joint scale-of-effect fit (experimental)
#> Data sets: counts (glm, n = 60); occ (glm, n = 80)
#> Kernel: gaussian
#> 
#> Scale structure
#>  covariate structure   data_sets shared_coefficient
#>       bin1    common counts, occ              FALSE
#>      cont2    common counts, occ              FALSE
#> 
#> Joint parameters (scales in map units)
#>    parameter estimate    SE 95% lower (profile) 95% upper (profile)
#>   sigma:bin1    63.66 12.79               42.97               96.93
#>  sigma:cont2   175.50 51.50               92.27              429.20
#> 
#> Apparent scale in each data set
#>  data_set covariate  sigma    SE 95% lower (profile) 95% upper (profile)
#>    counts      bin1  63.66 12.79               42.97               96.93
#>    counts     cont2 175.50 51.50               92.27              429.20
#>       occ      bin1  63.66 12.79               42.97               96.93
#>       occ     cont2 175.50 51.50               92.27              429.20
#> 
#> logLik = -137.046 (df = 8), AIC = 290.091, n = 140
#> Against separate scales: LR chi-square = 2.322, df = 2, p = 0.313
#> 
#> Note: Standard errors assume the data sets are independent given the covariates.

# Compare with one scale per data set
separate <- multiScale_joint(list(counts = fit_counts, occ = fit_occ),
                             scales = "separate", profile = FALSE,
                             verbose = FALSE)
AIC(joint, separate)
#>          df      AIC
#> joint     8 290.0912
#> separate 10 291.7697
# }
```
