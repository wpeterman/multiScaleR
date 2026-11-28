# Profile and Bootstrap Inference for an NB2 Scale

Examine the full allowed scale range of an NB2 fit with one spatial
predictor, retain separate likelihood regions, and optionally calibrate
likelihood-ratio tests by parametric bootstrap. This diagnostic is
experimental. The likelihood profile and the bootstrap decision retain
their numerical audits, and bootstrap confidence-set calibration is not
established. The diagnostic preserves the fitted object.

## Usage

``` r
scale_inference(
  object,
  sigma_values = NULL,
  level = 0.95,
  grid_points = 41L,
  bootstrap = 0L,
  bootstrap_values = NULL,
  seed = NULL,
  n_starts = 12L
)
```

## Arguments

- object:

  A `multiScaleR` fit from `backend = "profile_glm"` or
  `backend = "joint_nb2"`, with NB2 errors and exactly one optimized
  kernel mean or weighted surface SQ covariate. Weighted SQ requires
  retained raw buffers and a Gaussian or exponential kernel. A Bayesian
  fit from
  [`multiScale_bayes`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  is refused; its posterior already describes scale uncertainty (see
  [`posterior_scales`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md)).

- sigma_values:

  Optional positive scale values in map units to examine. The search
  bounds and the refined optimum are always included.

- level:

  Confidence level between zero and one; default `0.95`.

- grid_points:

  Number of log-spaced scale values when `sigma_values` is omitted;
  default `41`. Increase this to check whether narrow or disconnected
  regions were missed.

- bootstrap:

  Number of parametric bootstrap samples per calibrated scale; zero (the
  default) skips calibration. Small values have substantial Monte Carlo
  uncertainty. At least 999 is advisable for final tail estimation,
  although this can be expensive.

- bootstrap_values:

  Optional subset of scales to calibrate, in map units. They are added
  to the profile grid. Omitting this argument calibrates every grid
  point. A subset supplies pointwise tests, not a complete bootstrap
  confidence set.

- seed:

  Nonnegative integer random seed, required when `bootstrap > 0`. The
  caller's random-number state is restored.

- n_starts:

  Number of additional dispersed starts for inference refits; default
  `12`. Original starts, the selected fit, and both scale bounds are
  also included. Failed starts cannot win.

## Value

A list of class `msr_scale_inference` with:

- estimate:

  Refined scale estimate in map units.

- reference:

  Refined log-parameter vector, negative log likelihood, start audit,
  original negative log likelihood, and improvement. The original fitted
  object remains unchanged.

- grid:

  Scale values, profile negative log likelihood, likelihood-ratio
  statistic, chi-squared acceptance, bootstrap p-value and 95 percent
  Monte Carlo limits, bootstrap acceptance, and whether the Monte Carlo
  limits cross the rejection threshold. Uncalibrated entries are
  missing.

- profile_set:

  Separate accepted scale regions. Columns `lower` and `upper` are
  root-refined endpoints in map units; `lower_outside` and
  `upper_outside` retain adjacent rejected grid points. The two
  `*_bound_limited` columns mark endpoints limited by the optimization
  bounds.

- bootstrap_set:

  Accepted grid components with the same columns, or `NULL` unless every
  grid point was successfully calibrated. Endpoints are accepted grid
  values, not root-refined limits; outside brackets show their
  resolution.

- diagnostics:

  Zero-effect likelihood deviance and compatibility with the profile
  cutoff; profile status (`localized`, `disconnected`, `bound_limited`,
  or `unresolved`); original Wald availability; and bootstrap status,
  failure count, and Monte Carlo ambiguity count.

- bootstrap_records:

  One row per attempted bootstrap refit pair: generating scale,
  replicate, success, likelihood-ratio statistic, fitted scale and
  theta, boundary indicator, both likelihoods, attempted/eligible start
  counts, errors, and warnings. No failed replicate is silently dropped.

- profile_audit:

  Start records at each observed profile grid point.

- level, bootstrap, seed, bounds, call:

  Requested confidence level, bootstrap count and seed, map-unit search
  bounds, and matched call.

## Details

Coefficients and NB2 theta are reestimated at every constrained scale.
Weighted SQ uses the same refit and audit machinery. Existing
kernel-mean coverage studies do not establish interval coverage for
surface SQ or interaction models. Native factor contrasts and product
terms are retained in every constrained, unrestricted, and bootstrap
refit. Exponential-power shape is also reestimated. A stronger
unrestricted refit and a full-range grid can reveal an optimum missed by
the original fit; the returned reference records this without modifying
`object`. If conditional NB2 GLM iterations fail for a full-rank design,
the backend refines coefficients using the likelihood score and observed
information, checks convergence, and rebuilds the conditional GLM
result. For `joint_nb2` fits, these refits use the boundary-aware
conditional likelihood in `rho = log(1 + 1/theta)` coordinates,
including exact Poisson fits. Bootstrap samples use Poisson sampling
when the constrained dispersion estimate is on that boundary. This
preserves the boundary in inference without substituting a large finite
theta. These refits currently profile coefficients rather than using the
joint fitting optimizer.

The profile cutoff is `qchisq(level, 1)` for twice the log-likelihood
difference. Each detected accepted component is retained. Finite grids
can still miss narrow regions. A bound-limited endpoint is not a
resolved confidence limit.

The zero-effect diagnostic removes all response columns involving the
optimized kernel, including its products and factor interactions. If
that model is inside the profile likelihood cutoff, every admissible
scale has a feasible fit inside the cutoff. Scale is then reported as
unresolved. This is a likelihood diagnostic, not a calibrated test for
the presence of a kernel effect.

Bootstrap samples are generated at each constrained scale using its
fitted coefficients and theta, with the landscape, fixed design, and
offsets held fixed. Every sample is refitted both freely and at the
generating scale. The tail estimate is
`(1 + exceedances) / (bootstrap + 1)`. Monte Carlo limits transform
exact binomial limits by the same formula. A test accepts when its tail
estimate exceeds `1 - level`; the ambiguity flag is separate. Any failed
bootstrap pair makes that scale's calibration missing. Only unit
observation weights are currently supported for bootstrap.

A simulation with 32 datasets per scenario and 99 bootstrap draws per
test compared the bootstrap and chi-squared decisions at the known
generating scale. In the interior NB2 scenario, both included the
generating scale in 24 of 32 datasets and agreed on every decision. The
exact 95 percent binomial interval for inclusion was 0.566 to 0.885.
This measures pointwise inclusion in that simulation design. Coverage of
the returned finite-grid sets is a separate target.

In a separate 48-dataset simulation, all full-grid sets were available.
Their coverage changed with the grid even when the bootstrap decisions
at shared nodes were held fixed. A singleton component contains only its
accepted scale. Its outside brackets record adjacent rejected nodes.
[`vignette("v09-estimates-and-intervals", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v09-estimates-and-intervals.md)
summarizes the simulation results.

These are plug-in parametric bootstrap tests for independent
observations, not exact tests with known nuisance parameters. Neither
bootstrap nor positive curvature guarantees calibrated coverage under
weak identification or at parameter boundaries. Full-range,
disconnected, bound-limited, or Monte Carlo-ambiguous results should not
be reduced to a precise scale. Defaults of
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`summary()`](https://rspatial.github.io/terra/reference/summary.html),
and
[`predict()`](https://rspatial.github.io/terra/reference/predict.html)
are unchanged.

## See also

[`multiScale_optim`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`profile_sigma`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md),
[`diagnostics`](https://wpeterman.github.io/multiScaleR/reference/diagnostics.md)

## Examples

``` r
# \donttest{
if (requireNamespace("MASS", quietly = TRUE)) {
  set.seed(42)
  habitat <- terra::rast(nrows = 60, ncols = 60, xmin = 0, xmax = 600,
                         ymin = 0, ymax = 600, crs = "EPSG:32617")
  cells <- seq_len(terra::ncell(habitat))
  terra::values(habitat) <- sin(cells / 41) + cos(cells / 73)
  names(habitat) <- "hab"
  sites <- terra::vect(expand.grid(x = seq(150, 450, length.out = 6),
                                  y = seq(150, 450, length.out = 6)),
                       geom = c("x", "y"), crs = terra::crs(habitat))
  ki <- kernel_prep(sites, habitat, max_D = 120, bin = TRUE, verbose = FALSE)
  dat <- data.frame(ki$kernel_dat)
  dat$y <- rnbinom(nrow(dat), mu = exp(.5 + dat$hab), size = 2)
  initial <- glm(y ~ hab, data = dat, family = MASS::negative.binomial(2))
  fit <- multiScale_optim(initial, ki, backend = "profile_glm",
                           n_cores = 1, verbose = FALSE)
  assessment <- scale_inference(fit, grid_points = 9, n_starts = 2)
  assessment$profile_set   # Separate regions, with explicit search bounds.
  assessment$diagnostics  # A broad region may leave scale unresolved.
}
#> Warning: Joint Wald uncertainty unavailable: Poisson limit within the diagnostic likelihood threshold. Scale SEs are withheld for: hab. Inspect diagnostics and likelihood profiles.
#> $zero_effect_deviance
#> [1] 20.88188
#> 
#> $zero_effect_compatible
#> [1] FALSE
#> 
#> $profile_status
#> [1] "localized"
#> 
#> $original_wald_available
#> [1] FALSE
#> 
#> $bootstrap_status
#> [1] "not_requested"
#> 
#> $bootstrap_failures
#> [1] 0
#> 
#> $monte_carlo_ambiguous
#> [1] 0
#> 
# }
```
