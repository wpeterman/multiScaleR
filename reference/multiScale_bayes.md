# Estimate Scales of Effect with Bayesian Inference

Estimate landscape scales of effect and response coefficients together,
as one posterior distribution that combines what the data say with the
priors you state. Give each landscape predictor one shared scale, let
the scale change with measured covariates, or let it vary among groups
such as regions or years. This interface is experimental. A printed fit
reports "Experimental Bayesian landscape-scale model".

## Usage

``` r
# S3 method for class 'multiScaleR_bayes'
summary(object, ...)

# S3 method for class 'multiScaleR_bayes'
print(x, ...)

multiScale_bayes(
  formula,
  inputs,
  data,
  scales = NULL,
  family = c("nbinom2", "poisson", "binomial", "gaussian"),
  row_id = NULL,
  priors = list(),
  method = c("sample", "prepare", "laplace", "marginal", "slice"),
  chains = 4L,
  iter = 2000L,
  warmup = floor(iter/2),
  seed = NULL,
  control = list(),
  laplace_draws = 2000L,
  cores = 1L
)
```

## Arguments

- object:

  A `multiScaleR_bayes` object from `multiScale_bayes()`, prepared or
  fitted.

- ...:

  Not used. Supplying any further argument is an error.

- x:

  A `multiScaleR_bayes` object to print, prepared or fitted.

- formula:

  Two-sided response formula, such as
  `count ~ habitat + factor(year) + (1 | site)`. Name each landscape
  predictor as it is named in `inputs`. Its kernel-weighted value is
  computed at every candidate scale, so you never supply that column.
  The formula accepts fixed covariates, interactions with landscape
  predictors, offsets, and additive random intercepts. See "Write the
  response formula".

- inputs:

  Prepared landscape inputs from
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  or
  [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md),
  built with a Gaussian or exponential kernel. Kernel means need the
  binned summaries (`bin = TRUE`, the default). Optimized weighted
  surface SQ needs the retained raw cell buffers and source
  specifications. Means and SQ can be combined. Distances, priors, and
  reported scales all use the projected map units of these inputs.
  Declared fixed-radius surface and landscape metrics enter as ordinary
  numeric covariates with their original preparation normalization, so
  keep their prepared columns in `data`. No scale is estimated for them,
  and a constant or undefined metric is an error.

- data:

  Data frame with one row per observation, holding the response, fixed
  covariates, and grouping variables. Landscape predictors whose scales
  are estimated are rebuilt from `inputs` and need not appear here. A
  missing response or covariate value is an error. Rows are never
  dropped.

- scales:

  Named list with one one-sided formula for each landscape predictor in
  `formula`, such as `list(habitat = ~1 + (1 | year))`. Each formula
  models the natural log of that predictor's scale. The default `NULL`
  gives every predictor one shared scale, the same as `~1`. See "Write
  the scale formulas".

- family:

  Response family, one of `"nbinom2"` (default, negative binomial with a
  log link), `"poisson"` (log link), `"binomial"` (logit link), or
  `"gaussian"` (identity link). See "Choose a response family".

- row_id:

  Optional name of a column in `data` that matches each row to its
  prepared landscape, an entry of `inputs$binned$point_ids` (or a row
  name of `inputs$kernel_dat` for surface-only inputs). The default
  `NULL` matches the row names of `data`. IDs may repeat, so repeated
  visits can share one static landscape. With one map per year, use
  site-year IDs so each row reads the map for its own year. The name is
  kept with the fit and reused by
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  and [`simulate()`](https://rdrr.io/r/stats/simulate.html).

- priors:

  Named list of prior settings to change. Settings you leave out take
  their defaults, and the default
  [`list()`](https://rdrr.io/r/base/list.html) uses them all. See "Set
  the priors" for every setting and its default.

- method:

  Computation, one of `"sample"` (default, NUTS sampling), `"prepare"`,
  `"slice"`, `"laplace"`, or `"marginal"`. The method changes how the
  posterior is computed, not the model, likelihood, or priors.
  `"prepare"` checks the data, formulas, priors, and bins without
  fitting. See "Choose a computation method".

- chains:

  Number of independent chains for `"sample"` and `"slice"`; default 4.
  Chain diagnostics need at least two.

- iter:

  Iterations per chain, including warmup; default 2000. Raise it when
  effective sample sizes fall short.

- warmup:

  Iterations at the start of each chain that are discarded while the
  sampler tunes itself; default `floor(iter / 2)`.

- seed:

  Optional positive integer seed for sampling or approximation draws;
  default `NULL`. Slice, Laplace, and marginal fits restore your own
  random-number state after using a supplied seed. Slice sampling draws
  its chain seeds before any work is sent to workers, so a seed gives
  the same chains at every `cores` setting.

- control:

  Named list of method-specific settings; default
  [`list()`](https://rdrr.io/r/base/list.html) uses each method's
  defaults. NUTS accepts `adapt_delta` and `max_treedepth`. Slice
  sampling accepts no settings, so leave the list empty. See "Set
  computation controls" for every setting, including those for
  `"laplace"` and `"marginal"`.

- laplace_draws:

  Number of draws taken from a `"laplace"` or `"marginal"` approximation
  for summaries; default 2000, minimum 2. Without importance correction
  these are independent Gaussian or mixture draws. With it they are
  resampled from the weighted candidates. More draws reduce resampling
  noise in summaries. They do not improve the approximation or raise the
  importance effective sample size.

- cores:

  Positive integer limit on parallel R worker processes. The default of
  1 runs everything in the current session. Parallel fits need the
  optional package processx. See "Use parallel workers".

## Value

An object of class `multiScaleR_bayes`, a list that holds the frozen
model `specification`, the parameter layout in `parameters`, the `call`,
an `inference` label, and `inference_metadata`. Fitted objects add
posterior `draws`, `diagnostics`, and a record of the computation. The
class page
[multiScaleR_bayes](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md)
describes every element, shows where each parameter lives and on which
scale, and explains how to back-transform it. A prepared object
(`method = "prepare"`) has no draws.

[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
returns a list described in "Read the summary output".
[`print()`](https://rdrr.io/r/base/print.html) shows the computation,
family, number of observations, scale formulas, and stored diagnostic
status, and returns `x` invisibly.

## Choose a computation method

|  |  |  |
|----|----|----|
| `method` | What it computes | Optional packages |
| `"prepare"` | Matches rows to landscapes and checks formulas, priors, and bins. Fits nothing | None |
| `"sample"` (default) | The full posterior by NUTS, a gradient-based Markov chain Monte Carlo (MCMC) sampler, run through RTMB and tmbstan with the Laplace option off | RTMB, tmbstan, rstan, posterior |
| `"slice"` | The same full posterior by a slice sampler written in R | None to fit. posterior (\>= 1.6.0) for diagnostics |
| `"laplace"` | A joint Gaussian approximation around the posterior mode, in the same unconstrained, noncentered coordinates the samplers use | RTMB |
| `"marginal"` | Numerical integration over one scale-group SD, with conditional Laplace approximations for every other parameter | RTMB, plus posterior with importance correction |

[`check_bayes_setup()`](https://wpeterman.github.io/multiScaleR/reference/check_bayes_setup.md)
reports whether these packages are installed. Joint Laplace and
uncorrected marginal Laplace are exploratory. In matched simulations
they disagreed with full sampling about annual-scale probabilities even
when their own numerical checks passed. Confirm every quantity you
report with `"sample"` or `"slice"`.
[`vignette("v15-bayesian-computation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v15-bayesian-computation.md)
compares the methods on one model and shows what to do when a fit
struggles.

Preparation, NUTS, slice sampling, and joint Laplace accept several
landscape predictors and additive or crossed random intercepts. Marginal
Laplace, with or without importance correction, needs exactly one
scale-group SD across all landscape predictors. A second random scale
term adds another SD even when it uses the same grouping variable. The
function checks this before the numerical fit begins and never switches
methods for you. A model that the syntax accepts can still fail to
converge, leave a scale weakly identified, or give intervals of unknown
calibration.

## Choose a response family

- `"nbinom2"`: counts with a log link. The variance is \\\mu +
  \mu^2/\theta\\, where theta is the size (dispersion) parameter
  reported by
  [`summary()`](https://rspatial.github.io/terra/reference/summary.html).

- `"poisson"`: counts with a log link.

- `"binomial"`: a logit link. Supply 0 and 1 values, `FALSE` and `TRUE`,
  a two-level factor whose second level is a success, or
  `cbind(successes, failures)` for grouped trials. Response-scale
  predictions are success probabilities.

- `"gaussian"`: a continuous response with an identity link and an
  estimated residual SD.

## Write the response formula

The response formula accepts the following terms.

- Fixed covariates, factors, and offsets.

- Additive random intercepts such as `(1 | site)`. Each grouping
  variable can appear once in a formula and needs at least two levels.
  Use one for repeated observations at a site.

- Products, crossings, nesting, and factor interactions with a landscape
  predictor, such as `habitat * year` or `habitat:forest`.
  `habitat * year` estimates a separate habitat association in each
  year. It changes the response coefficients only and adds no year
  effect to the scale formula. An interaction between two landscape
  predictors still has one scale model per predictor. Every expanded
  coefficient keeps the normal prior on R's native factor and contrast
  coding.

- Transformations of fixed covariates.
  [`scale()`](https://rspatial.github.io/terra/reference/scale.html),
  [`poly()`](https://rdrr.io/r/stats/poly.html),
  [`splines::ns()`](https://rdrr.io/r/splines/ns.html), and
  [`splines::bs()`](https://rdrr.io/r/splines/bs.html) keep their
  training definitions, including inside nested expressions such as
  `I(scale(temp))` and in transformed offsets. Prediction reuses those
  constants for single rows, batches, and effect curves. Give a
  namespace prefix to any custom function that reuses one of these
  names, or compute its predictor before fitting.

- Arithmetic on landscape predictors and literal numbers, such as
  `I(habitat^2)`. The operators `+`, `-`, `*`, `/`, and `^` and the
  functions [`exp()`](https://rdrr.io/r/base/Log.html),
  [`log()`](https://rdrr.io/r/base/Log.html),
  [`sqrt()`](https://rdrr.io/r/base/MathFun.html),
  [`sin()`](https://rdrr.io/r/base/Trig.html),
  [`cos()`](https://rdrr.io/r/base/Trig.html), and
  [`tan()`](https://rdrr.io/r/base/Trig.html) are allowed. These
  expressions act on the standardized kernel value after smoothing, not
  on raster cells, and keep the fixed reference standardization and
  coefficient priors.

- Cubic splines of a landscape predictor, from
  [`splines::ns()`](https://rdrr.io/r/splines/ns.html) or
  [`splines::bs()`](https://rdrr.io/r/splines/bs.html) with literal
  fixed knots. See "Use splines of a landscape predictor".

The following stop with an error: a fixed covariate inside an expression
with a landscape predictor (write the interaction `habitat:forest`
rather than `I(habitat * forest)`), any other function, including
nondifferentiable ones, a landscape predictor in an offset, and adaptive
spline recipes. An expression that is not finite is also an error.
Values are never clipped. Random slopes, spatial random fields,
observation weights, zero inflation, and imperfect detection are not
supported. `method = "prepare"` checks all of this before any optional
package is needed.
[`vignette("v07-fitting-methods", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v07-fitting-methods.md)
lists the covariate types each fitting route accepts.

### Use splines of a landscape predictor

Supply literal `knots` and two increasing `Boundary.knots`, with every
knot distinct and strictly inside the boundaries, for example
`splines::ns(habitat, knots = 0, Boundary.knots = c(-3, 3))`. Knots are
in the units of the standardized predictor, after the fixed reference
transformation. To place knots at raw kernel-weighted values, prepare
the model once with the plain predictor, then subtract
`prepared$specification$center` and divide by
`prepared$specification$spread` for both the knots and the boundaries.
Write the resulting numbers into the formula. Never learn knots from
prediction data.

Basis values and their derivatives update at every candidate scale,
while the knots, the coefficient basis, and the priors stay fixed. All
basis columns share the scale of their source predictor, including in
factor interactions. Natural-spline tails are linear on the link scale.
These are unpenalized bases of fixed complexity, with no
smoothing-variance uncertainty and no monotonicity constraint. Adaptive
`df` recipes, repeated or boundary knots, and degrees other than 3 are
not supported here. The error message points to the refit backend of
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
for those.

## Write the scale formulas

Each formula in `scales` models the natural log of one predictor's
scale. Scale formulas use columns of `data`. They cannot use landscape
predictors, offsets, or a rank-deficient design.

- `~1` (the default): one shared scale for all observations.

- `~ elevation`: the log scale changes linearly with a numeric
  covariate. [`exp()`](https://rdrr.io/r/base/Log.html) of its
  coefficient is the multiplicative change in scale per unit of the
  covariate.

- `~ year`: fixed effects of a factor on the log scale, coded as
  contrasts with the first level.

- `~0 + year`: a separate log scale for each level, each with the same
  distance prior, and no pooling. An intercept-free formula must define
  disjoint factor levels like this. Every other formula needs an
  intercept.

- `~1 + (1 | year)`: partial pooling. Years share a typical scale and an
  estimated between-year SD.

Random intercepts are noncentered, so each group's deviation is a
standard-normal latent value times the group SD. For `~1 + (1 | year)`
the scale of year \\y\\ is \\\log D_y = \alpha + \tau z_y\\, with \\z_y
\sim N(0, 1)\\. Here alpha is the log of the typical scale and tau is
the SD of log scale among years. `exp(tau)` is a multiplier. A year one
SD above typical has a scale `exp(tau)` times the typical one. Years are
exchangeable, and the model imposes no trend or temporal
autocorrelation.

The two formulas describe different parts of the model. A year effect in
the response formula changes baseline abundance or occurrence. A year
effect in the scale formula changes how far the landscape effect
reaches.

### Expand a dot in a scale formula

A dot expands every column of `data`, including the response and the
row-ID column. Remove the response, row keys, and unrelated columns
explicitly, as in `~ . - count - site_year`. Only the active terms
create scale coefficients and priors. A removed column adds no
coefficient and no latent group coordinate. Prediction reuses the fitted
active terms and does not need the removed columns. A landscape
predictor of the response formula cannot also enter a scale formula.

## Set the priors

|  |  |  |
|----|----|----|
| Setting | Default | What it controls |
| `scale_median` | `sqrt(inputs$min_D * inputs$max_D)` | Prior median of each scale, in map units. One positive number, or a named vector with one value per landscape predictor. It is the center of the normal prior on each log-scale intercept, or on each coefficient of a `~0 + year` formula |
| `scale_intercept_sd` | 1 | SD of that normal prior, on the natural log scale |
| `scale_coef_sd` | 0.5 | SD of the zero-centered normal prior on other scale coefficients (slopes and log-distance contrasts) |
| `scale_group_sd` | 0.5 | Scale of the half-normal prior on each between-group SD of log scale |
| `beta_sd` | 2.5 | SD of the zero-centered normal prior on each response coefficient, on the link scale |
| `response_group_sd` | 1 | Scale of the half-normal prior on each response random-intercept SD, on the link scale |
| `dispersion_log_mean`, `dispersion_log_sd` | 0, 1 | Mean and SD of the normal prior on the log of NB2 theta or of the Gaussian residual SD |

The scale priors act on log distances and set no hard limit at `min_D`
or `max_D`. A scale beyond `max_D` weights the extracted neighborhood
ever more evenly. It cannot draw on cells outside the prepared buffers.
[`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md)
reports the posterior mass below `min_D` and above `max_D`. Check what
the priors imply with
[`simulate_bayes_prior()`](https://wpeterman.github.io/multiScaleR/reference/simulate_bayes_prior.md)
before you fit.
[`vignette("v13-bayesian-scales", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v13-bayesian-scales.md)
converts these settings into distances and multipliers.

## Understand the fixed standardization

Kernel means use the finite-cell numerator and denominator stored in the
bins. Missing raster cells are excluded and the remaining weights
renormalized. Each landscape predictor is then centered and scaled once,
across the fitting rows, at its prior median scale (`scale_median`).
That reference transformation is held fixed for every candidate scale,
group, and draw, including when a site has repeated observations or
scales differ among years. It is never estimated again from prediction
rows. The retained values are `specification$center` and
`specification$spread`.

[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
restandardizes at each fitted scale, so the two functions give
coefficients with different meanings. Here a response coefficient is the
change in the linear predictor for one reference SD of the
kernel-weighted predictor at `scale_median`, not at the fitted scale.
Its prior therefore has a fixed meaning, but its value changes when you
change `scale_median`. Report `scale_median` with any coefficient, and
do not compare these coefficients directly with those of
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).

Weighted surface SQ is the square root of the kernel-weighted centered
variance, divided by the sum of the weights, computed over each source's
finite cells. It uses the same pooled reference transformation.

### Plan memory for weighted SQ

Preparation records `specification$surface$memory`, which gives the
number of observation-cell evaluations (`observation_cells`), the
retained buffer bytes (`retained_bytes`), and a planning estimate of
differentiation memory (`estimated_AD_bytes`). Repeated observations
count separately. Before a differentiation tape is recorded for SQ, a
conservative check against installed RAM must pass, and the fit stops
otherwise. The estimate is a planning figure, not a quota or a
guarantee, so leave headroom. Candidate batches are evaluated in bounded
pieces. Other optimized surface metrics are not supported.

## Read the summary output

[`summary()`](https://rspatial.github.io/terra/reference/summary.html)
of a prepared object returns `inference`, `family`, `observations`, and
`scale_formulas`. For a fit with draws it returns a list with these
elements.

- `inference`, `family`, `observations`: the computation label, the
  response family, and the number of observations.

- `coefficients`: a data frame with one row per response `coefficient`
  and columns `mean`, `median`, `sd`, `lower`, and `upper` (a central
  95% credible interval). Values are on the link scale, the log expected
  count for Poisson and NB2, the log odds for binomial, and response
  units for Gaussian models. For a count model,
  [`exp()`](https://rdrr.io/r/base/Log.html) of a landscape coefficient
  is the multiplicative change in the expected count per reference SD of
  the predictor.

- `scales`: the
  [`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md)
  table, with each scale in map units.

- `scale_coefficients`: a list with one data frame per landscape
  predictor. Intercepts and `~0 + year` coefficients are in log map
  units, so [`exp()`](https://rdrr.io/r/base/Log.html) gives a distance.
  Other coefficients are log-distance contrasts, so
  [`exp()`](https://rdrr.io/r/base/Log.html) gives a ratio of distances.

- `group_sd`: a list with `log_scale_group_sd`, the between-group SD of
  log scale (tau) for each scale random term labeled `predictor:group`
  in `term`, and `response_group_sd`, the random-intercept SDs on the
  link scale. Each is `NULL` when the model has no such term. These SDs
  measure variation among groups. They are not the standard error of any
  one group's scale.

- `dispersion`: NB2 theta or the Gaussian residual SD on its natural
  scale, or `NULL` for Poisson and binomial models.

- `diagnostics`: the fit's stored checks, described on the class page
  [multiScaleR_bayes](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md).

Summary intervals are marginal central 95% intervals for individual
quantities, not a joint band. `mean` and `sd` are summaries of the
draws. For long-tailed distances, prefer `median`, `lower`, and `upper`.
For an importance-corrected fit, the `scales` table uses the resampled
draws. `posterior_scales(fit, estimator = "weighted")` summarizes the
weighted candidates directly, as
[`plot()`](https://rspatial.github.io/terra/reference/plot.html) does by
default.

The stored check covers the raw parameters and, for NUTS, sampler
behavior. It does not replace
[`diagnose_bayes_scales()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_scales.md)
for the scales you report, or a posterior predictive check with
[`simulate.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/simulate.multiScaleR_bayes.md).

## Interpret weakly informed scales

A few years or regions say little about variation among them, so the
group-SD prior can drive the answer. Compare plausible scale and
group-SD priors, inspect posterior predictive behavior, and check the
sampling diagnostics. A response coefficient near zero leaves its
landscape scale weakly identified even when the chains mix well. The
half-normal prior on a group SD has no point mass at zero, so a positive
lower interval limit is not a test of equal scales. Compare a
shared-scale model instead. Numerical accuracy depends on the response,
the scale design, and the method. Interval calibration for general
ecological data sets has not been established.

## Use parallel workers

`cores` sets how many local R processes run at once. NUTS and slice
sampling spread chains across workers, so the useful limit is `chains`.
Joint and uncorrected marginal Laplace spread optimizer starts, so the
limit is the number of starts. Marginal integration still visits grid
nodes in order, so each node can start from the previous fits. With
importance correction, workers also evaluate batches of candidate
posterior densities, while proposal generation and diagnostics stay in
the main process. Preparation uses one process.

Each worker holds its own copy of the model data and, where needed, its
own differentiation tape, so memory use rises with `cores`. Temporary
worker files are removed after fitting. Worker start-up and data
transfer can make small fits slower. The effective limit is recorded in
`sampling$cores`, `laplace$cores`, or `marginal$cores`.

## Set computation controls

|  |  |  |  |
|----|----|----|----|
| Setting | Method | Default | What it does |
| `adapt_delta` | `"sample"` | 0.9 | Target acceptance rate during NUTS adaptation, strictly between 0 and 1. Raise it when transitions diverge |
| `max_treedepth` | `"sample"` | 12 | Maximum NUTS tree depth. Raise it when the sampler hits the limit |
| `starts` | `"laplace"`, `"marginal"` | 6 for Laplace, 3 per grid node for marginal | Deterministic, dispersed optimizer starts. Every start adds run time |
| `iter.max` | `"laplace"`, `"marginal"` | 1000 | Optimizer iterations per start |
| `max_grid` | `"marginal"` | 129 | Largest integration grid, 33, 65, 129, or 257 nodes. More nodes cost more time |
| `cdf_tol`, `quantile_tol`, `logz_tol` | `"marginal"` | 0.005, 0.02, 0.01 | Largest change allowed between successive grids in the SD CDF, its 2.5%, 50%, and 97.5% quantiles (in prior SDs), and the log normalizing constant |
| `correction` | `"marginal"` | `"none"` | `"importance"` reweights the approximation against the full posterior |
| `importance_min` | importance | 16384 | Minimum number of production candidates retained |
| `importance_max` | importance | 65536 | Cap on production candidates, 16384, 32768, 65536, or 131072 |
| `importance_df` | importance | 7 | Degrees of freedom of the Student t proposal, 7 or 5 |
| `importance_refinements` | importance | 1 | Further proposal adaptations allowed after a failed stage, 0, 1, or 2 |
| `nested_response` | importance | `TRUE`, or `FALSE` with `conditional_beta = TRUE` | Integrate response random effects when locating conditional modes |
| `conditional_beta` | importance | `FALSE` | Experimental. Propose response coefficients conditional on each candidate's scales and response random effects |

Grid refinement stops when all three changes fall within their
tolerances. It measures numerical integration error, not Laplace error.

`importance_max` is a cap, not a target. A population that passes its
checks at a smaller size stops there. To require more candidates, as you
may for precise annual quantiles, also raise `importance_min`, which
takes the same values and cannot exceed `importance_max`. A separate
training pilot of 16384 candidates always runs first and is discarded.
Larger populations cost more, and they guarantee neither a chosen Monte
Carlo error nor recovery of a missed posterior region.
`importance_df = 5` gives heavier proposal tails, which can help when
importance ratios are unstable, though it can need more candidates. With
`importance_refinements = 1`, a failed first production stage of 16384
candidates trains one further proposal. Its candidates are discarded and
a fresh final stage is drawn. These settings change sampling efficiency,
not the priors or the diagnostic thresholds.

`nested_response = FALSE` keeps the original joint conditional proposal,
which can miss response-SD uncertainty. An explicit
`nested_response = TRUE` requires importance correction.
`conditional_beta = TRUE` also requires importance correction and
selects `nested_response = FALSE` unless you supply it, so an explicit
`nested_response = TRUE` is then an error. Gaussian coefficient
conditionals are exact. Other families use a normalized Student t
proposal with 7 degrees of freedom, and a fixed 20% latent-prior
component guards regions the local geometry misses. Conditional
optimization costs more per candidate and does not guarantee precise
interval endpoints.
[`vignette("v15-bayesian-computation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v15-bayesian-computation.md)
reports the simulation evidence for the corrected proposals and its
limits.

## How each computation works

### Slice sampling

Slice sampling updates blocks of parameters in turn from the full
posterior. For Gaussian responses it draws the regression coefficients
from their exact conditional distribution. It reuses kernel values when
it updates response parameters and recomputes only the changed kernel
during scale updates. Efficiency depends on the data set and on
dependence among scales and response effects, so compare Rhat and
effective sample size per second before choosing a sampler for a larger
fit.

### Joint Laplace approximation

The Laplace method inverts the positive-definite Hessian of the full
negative log posterior at a numerically stationary mode, including the
Jacobians of the log-SD coordinates. It approximates all parameters
jointly. It neither fixes the estimated pooling SDs nor integrates
random effects with a nested approximation. Non-positive curvature, or
failure to find a stationary mode, stops the fit with an error whose
`starts` element keeps every attempted solution. Catch it with
[`tryCatch()`](https://rdrr.io/r/base/conditions.html) to inspect them.
Distinct modes within two log-posterior units trigger a warning.
Skewness, weak scale identification, and dependence between group SDs
and year effects can make the pooling and the intervals inaccurate even
when the checks pass. Compare important results with
`method = "sample"`, especially with few years or weak landscape
effects. `chains`, `iter`, and `warmup` are unused.

### Marginal Laplace approximation

The marginal method needs exactly one scale-group SD, such as the
between-year SD in `list(habitat = ~1 + (1 | year))`. It integrates that
SD numerically on a grid in the coordinates of its half-normal prior
CDF, which leaves the prior unchanged and allows a skewed SD
distribution. Every other parameter, including any response
random-intercept SD, gets a conditional joint Gaussian approximation at
each grid node, and the conditional posterior volume corrects the SD
marginal. This is not full INLA. Linear interpolation between nodes and
triangular mixture sampling give continuous SD draws and carry
conditional dependence into the other parameters. The grid endpoints use
SDs of `1e-8` and 8 times the prior scale. These approximate the prior
CDF endpoints and do not bound the sampled SDs. The checks also fail
when more than 0.001 of the posterior mass lies above prior CDF 0.9999,
or more than 0.01 of the node weight carries competing conditional
modes. A failed node stops the fit rather than being dropped.
Conditional skewness, remaining multimodality, and response-SD
uncertainty still call for comparison with full sampling.

### Importance correction

`control = list(correction = "importance")` targets the same full
posterior as NUTS by reweighting candidates from a proposal built on the
marginal approximation. It does not rescale intervals by an empirical
coverage factor. A training pilot draws from an equal mixture of
multivariate and coordinatewise Student t distributions (7 degrees of
freedom, or 5 with `importance_df = 5`, inflated by a factor of 1.2) and
fits local proposal means and covariances. The pilot is then discarded.
Production candidates come 80% from the learned proposal and 20% from
the original one. Training weights are tempered and learned covariances
are shrunk toward the original, which stabilizes learning without
changing the target. Production weights are neither smoothed nor
truncated, and every mixture density, including the prior-CDF Jacobian,
enters the importance ratio. A 1% uniform component in prior-CDF
coordinates keeps the proposal positive in very small tails. A further
1% normal component for the log SD, with mean `log(2 * prior SD)` and SD
0.6, covers SDs beyond the numerical CDF endpoints. A refinement after a
failed stage fits a new component with 60% of the proposal weight and
keeps 40% for the previous mixture. Training stages never enter the
final estimate, though `history` and `production_candidates` keep their
diagnostics and cost.

With nested response integration, outer curvature comes from central
differences of RTMB's analytic gradient. With exactly one response
random-intercept term, each candidate draws its site effects from
independent Student t proposals with 20 degrees of freedom, centered on
their conditional modes and scaled by their curvatures. This follows the
dependence of site effects on the sampled response SD. Only the smaller
outer covariance is stored at each SD node, which avoids dense
site-by-site matrices. More response terms use an expanded joint
proposal, and its reliability has not been established. Nonfinite site
geometry falls back to the normalized prior proposal for that candidate
and is counted. Candidates are never dropped. Unusable learned
covariances keep their original nodes. The full-posterior ratios still
correct these proposals, and the same checks apply.

The corrected checks require Pareto k at most 0.7, raw weight effective
sample size (ESS) at least 1000, tail-indicator ESS at least 400,
maximum weight at most 0.01, and split-half CDF differences at most
0.05. Indicators are evaluated at the 2.5%, 50%, and 97.5% quantiles of
the unconstrained parameters, annual log scales, and group SDs. These
are importance diagnostics, not MCMC Rhat or ESS. The original proposal
flag stays in `diagnostics$proposal_passed`, but only the corrected
checks decide eligibility. A finite candidate set can still miss modes
or tails. A warning or a failed check calls for NUTS, and passing does
not establish interval calibration.

## See also

The class page
[multiScaleR_bayes](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md)
for the returned object;
[`check_bayes_setup()`](https://wpeterman.github.io/multiScaleR/reference/check_bayes_setup.md),
[`simulate_bayes_prior()`](https://wpeterman.github.io/multiScaleR/reference/simulate_bayes_prior.md),
[`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md),
[`plot.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/plot.multiScaleR_bayes.md),
[`diagnose_bayes_scales()`](https://wpeterman.github.io/multiScaleR/reference/diagnose_bayes_scales.md),
[`predict.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR_bayes.md),
[`predict_scales()`](https://wpeterman.github.io/multiScaleR/reference/predict_scales.md),
[`simulate.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/simulate.multiScaleR_bayes.md),
[`cross_validate_bayes()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_bayes.md),
[`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md),
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
and
[multiScaleR-glossary](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-glossary.md)
for terms;
[`vignette("v13-bayesian-scales", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v13-bayesian-scales.md)
for fitting and reading a model,
[`vignette("v14-bayesian-prediction", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v14-bayesian-prediction.md)
for prediction, and
[`vignette("v15-bayesian-computation", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v15-bayesian-computation.md)
for choosing and checking a computation method.

## Examples

``` r
# Annual survey example: 100 counts at 25 sites over five years, with one
# habitat map per year prepared by kernel_prep_by_group()
example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
  package = "multiScaleR"))

# Check the model before fitting. Preparation matches each row to its
# landscape through row_id and needs no optional packages.
prepared <- multiScale_bayes(count ~ habitat + factor(year) + (1 | site),
  example$inputs, example$data, scales = list(habitat = ~1 + (1 | year)),
  row_id = "site_year", family = "nbinom2",
  priors = list(scale_median = 250), method = "prepare")
prepared                        # family, observations, and scale formula
#> Experimental Bayesian landscape-scale model
#> Inference: prepared 
#> Family: nbinom2 | observations: 100 
#> Scale formulas:
#>   habitat : ~1 + (1 | year) 
prepared$inference_metadata$uncertainty$intervals  # "not_provided" yet
#> [1] "not_provided"

# The same model fitted with NUTS, stored as 2,000 thinned draws
fit <- example$fits$pooled
fit_summary <- summary(fit)
fit_summary$diagnostics$passed    # TRUE: the sampler checks passed
#> [1] TRUE
# Response coefficients on the log expected-count scale
fit_summary$coefficients[c("coefficient", "median", "lower", "upper")]
#>        coefficient    median       lower     upper
#> 1      (Intercept) 2.5875680  2.34823418 2.8156719
#> 2 factor(year)2022 0.3441289  0.09189045 0.5941000
#> 3 factor(year)2023 0.1302952 -0.21349941 0.4921239
#> 4 factor(year)2024 0.1284989 -0.12032150 0.3825631
#> 5 factor(year)2025 0.4055208  0.11122988 0.6711479
#> 6          habitat 0.9068686  0.74879486 1.2338224
# Sigma for each year, in meters, after partial pooling
fit_summary$scales[c("group", "median", "lower", "upper")]
#>       group   median    lower    upper
#> 1 year=2021 173.0828 39.08530 431.8495
#> 2 year=2022 173.7422 41.71146 423.5986
#> 3 year=2023 143.6104 34.01918 376.3533
#> 4 year=2024 175.0161 41.81391 441.1051
#> 5 year=2025 155.9366 49.68400 373.9298
# Tau, the SD of log sigma among years, and the multiplier exp(tau)
fit_summary$group_sd$log_scale_group_sd[c("term", "median", "lower", "upper")]
#>           term    median      lower     upper
#> 1 habitat:year 0.2351976 0.01032759 0.8547712
exp(fit_summary$group_sd$log_scale_group_sd$median)
#> [1] 1.265159

# A slice-sampled fit with one shared sigma per layer. Its call shows the
# sampling arguments, and its scale intercepts are in log meters.
shared <- readRDS(system.file("extdata", "cache-bayesian-scales.rds",
  package = "multiScaleR"))$fit
shared$call
#> multiScale_bayes(formula = counts ~ site + land1 + land2, inputs = inputs, 
#>     data = landscape_counts, family = "poisson", priors = list(scale_median = c(land1 = 300, 
#>         land2 = 300)), method = "slice", seed = 5211)
shared_coef <- summary(shared)$scale_coefficients
exp(sapply(shared_coef, function(x) x$median))   # median sigma in meters
#>    land1    land2 
#> 246.6455 501.1973 
# See vignette("v13-bayesian-scales", package = "multiScaleR") for fitting
# and interpretation.
```
