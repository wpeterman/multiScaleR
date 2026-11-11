# Multiscale optimization

Identifies the kernel scale of effect (sigma) for each raster covariate
by maximizing the log-likelihood of a fitted statistical model.
Repeatedly replaces kernel-weighted covariate values at different scales
and refits the model, using `optim` (or `optimParallel` for parallel
execution) with the L-BFGS-B algorithm.

## Usage

``` r
multiScale_optim(
  fitted_mod,
  kernel_inputs,
  join_by = NULL,
  par = NULL,
  start_strategy = c("single", "screen"),
  screen_n_sigma = 5,
  screen_n_jitter = 6,
  screen_maxit = 8,
  screen_jitter_sd = 0.5,
  n_cores = NULL,
  PSOCK = FALSE,
  verbose = TRUE,
  refit_fn = NULL
)
```

## Arguments

- fitted_mod:

  A fitted model object whose covariates include the kernel-weighted
  variables defined in `kernel_inputs`. Supported classes include `lm`,
  `glm`, `gls` (nlme), and `unmarkedFit` (unmarked). Many other model
  classes are also supported via
  [`stats::update()`](https://rdrr.io/r/stats/update.html) or a custom
  `refit_fn`. A predictor cannot be named `"x"`, which is reserved for
  internal coordinate data.

- kernel_inputs:

  A list of class `"multiScaleR_data"` created by
  [`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md).
  Must contain elements `raw_cov`, `d_list`, `min_D`, `max_D`,
  `unit_conv`, and `kernel`.

- join_by:

  Default: `NULL`. For `unmarked` models, a data frame with one unique
  key row per spatial site in `kernel_inputs`, in the same order as
  those sites. Its named key columns must also occur in the fitted
  model's `siteCovs`. Use it when model rows repeat sites, such as one
  row per site and year. See Details.

- par:

  Optional numeric vector of starting values for the optimizer. Values
  must be divided by `max_D` to match the internal scaled parameter
  space. Length must equal the number of optimized covariates (or twice
  that for `kernel = "expow"`, where shape parameters follow sigma
  parameters). Default: `NULL`; starting values are chosen
  automatically.

- start_strategy:

  Character. How to choose starting values when `par = NULL`: `"single"`
  (default) uses one shared default start, while `"screen"` performs a
  low-cost prescreen of sigma values and short screening optimizations
  before launching one full optimization.

- screen_n_sigma:

  Positive integer. Number of log-spaced sigma values evaluated per
  covariate during the prescreen when `start_strategy = "screen"`.
  Default: `5`.

- screen_n_jitter:

  Non-negative integer. Number of jittered candidate starts evaluated
  with short screening optimizations after the marginal prescreen.
  Default: `6`. Set to `0` to skip the short jittered runs.

- screen_maxit:

  Positive integer. Maximum iterations used for each short screening
  optimization. Default: `8`.

- screen_jitter_sd:

  Non-negative numeric scalar. Standard deviation of multiplicative
  sigma jitter on the log scale during screening. Default: `0.5`.

- n_cores:

  Positive integer. Number of cores for parallel optimization via
  `optimParallel`. Default: `NULL` (single-threaded). When
  `start_strategy = "screen"`, the same `n_cores` setting is also used
  for short screening attempts before the full optimization. Parallel
  optimization is beneficial for models with many covariates or slow
  log-likelihood evaluations. When `n_cores > 1`, `multiScale_optim()`
  checks
  [`estimate_multiscale_ram`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md)
  and stops before screened starts or worker creation if the requested
  workers exceed the conservative RAM-based limit for the current
  inputs. Use `n_cores = 1` for serial execution, or rebuild inputs with
  `kernel_prep(..., bin = TRUE, store_cell_data = FALSE)` for lean mode
  when every covariate is kernel-weighted. Landscape and surface metrics
  cannot use lean mode. An unavailable RAM estimate does not block the
  run; a passing estimate does not guarantee sufficient memory. No
  worker count is reduced automatically.

- PSOCK:

  Logical. On Windows, a PSOCK cluster is always used. On Unix, a FORK
  cluster is used by default (faster). Set `TRUE` to force a PSOCK
  cluster on Unix. Default: `FALSE`.

- verbose:

  Logical. Print optimization status and warnings to the console.
  Default: `TRUE`.

- refit_fn:

  Optional function for refitting `fitted_mod` during optimization. When
  `NULL` (default),
  [`stats::update()`](https://rdrr.io/r/stats/update.html) (or the
  unmarked equivalent) is used. Provide this only when the default refit
  path is insufficient. See Details.

## Value

A list of class `"multiScaleR"` containing:

- `scale_est`:

  Data frame with one row per optimized covariate and two columns:
  `Mean` (optimized sigma on the original projection scale) and `SE`
  (Hessian-based standard error). `SE` is `NA` for hard-radius landscape
  and unweighted surface metrics, whose likelihood changes in steps. Row
  names are covariate names.

- `shape_est`:

  Data frame of the same structure as `scale_est` for the shape
  parameter, or `NULL` when `kernel != "expow"`.

- `optim_results`:

  The raw list returned by
  [`stats::optim`](https://rdrr.io/r/stats/optim.html) or
  [`optimParallel::optimParallel`](https://rdrr.io/pkg/optimParallel/man/optimParallel.html),
  including `par`, `value` (negative log-likelihood), `hessian`, and
  convergence codes. Also contains `par_unscale` (sigma values back on
  the projection scale) and `hessian_unscale`.

- `opt_mod`:

  The refitted model at the optimized sigma values. This is the model
  object to use for inference, prediction, and
  [`plot_marginal_effects`](https://wpeterman.github.io/multiScaleR/reference/plot_marginal_effects.md).

- `fitted_mod_original`:

  The original `fitted_mod` passed by the user, stored for reference.

- `min_D`:

  Numeric. Lower bound for sigma during optimization.

- `max_D`:

  Numeric. Upper bound for sigma during optimization.

- `kernel_inputs`:

  The `kernel_inputs` list (minus `min_D` and `max_D`, which are stored
  separately).

- `scl_params`:

  Named list with `mean` and `sd` vectors for each covariate: the
  centering and scaling parameters from the optimized kernel data. Used
  by
  [`kernel_scale.raster`](https://wpeterman.github.io/multiScaleR/reference/kernel_scale.raster.md)
  when `scale_center = TRUE`.

- `join_by`:

  The `join_by` data frame, or `NULL`.

- `opt_context`:

  Internal optimization context object storing the model-class-specific
  refit logic. Retained for use by
  [`profile_sigma`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md).

- `profile_scale_est`:

  Profile-likelihood CI data frame, or `NULL` until
  [`profile_sigma`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md)
  has been run.

- `diagnostics`:

  List of diagnostic objects; see
  [`diagnostics`](https://wpeterman.github.io/multiScaleR/reference/diagnostics.md).

- `next_run`:

  List of recommended values for a follow-up fit. Includes `max_D` for
  the next
  [`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  call, optimized `start_sigma` values in map units, optional
  `start_shape` values for `kernel = "expow"`, and `start_par` on the
  internal scaled parameter space expected by `multiScale_optim`.
  `n_cores` is `NULL`; a conservative parallel worker-count suggestion
  is available on demand from
  [`estimate_multiscale_ram`](https://wpeterman.github.io/multiScaleR/reference/estimate_multiscale_ram.md)
  (it is no longer computed automatically, as that query can be slow and
  shells out to an external process).

- `warn_message`:

  Integer vector of triggered warning codes: 1 = max-distance, 2 = sigma
  precision, 3 = shape precision.

- `call`:

  The matched call.

## Details

**Optimization approach**

The optimizer uses the L-BFGS-B algorithm (bounded quasi-Newton) to
maximize the log-likelihood over sigma (and shape for `"expow"`) while
holding all regression coefficients at their fitted values for each
candidate scale. Standard errors are Hessian-based approximations from
the outer optimization. Summary methods report profile-likelihood
confidence intervals for sigma when
[`profile_sigma`](https://wpeterman.github.io/multiScaleR/reference/profile_sigma.md)
has been run on the object; otherwise they fall back to Hessian-based
intervals. Candidate scales that make a covariate undefined at a modeled
observation are rejected, so the optimizer compares likelihoods for the
same sample. The final `sample_size` diagnostic reports prepared and
fitted row counts; for repeated-site `unmarked` data it also reports the
number of distinct prepared spatial sites. For hard-radius landscape and
unweighted surface metrics, Hessian uncertainty is not meaningful; use
`summary(x, profile = TRUE)` for profile-likelihood intervals.

**Binned acceleration**

When `kernel_inputs` was created by
[`kernel_prep`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
with `bin = TRUE` (the default), kernel-weighted covariates are
evaluated from precomputed distance-binned summaries rather than by
iterating every buffer cell on each optimizer evaluation. This makes the
per-evaluation cost independent of `max_D` and typically speeds up
optimization by an order of magnitude or more for large buffers, with a
negligible binning approximation. Inputs created with
`store_cell_data = FALSE` ("lean" mode) carry only the binned summaries;
these still optimize, profile, and refit normally. Landscape-metric
covariates always use the exact per-cell path and therefore require
cell-level data.

**Parallel optimization**

To ensure that fitted model function calls are properly serialized for
parallel workers, fit models using fully namespace-qualified functions.
For example, fit a negative binomial model as
`fitted_mod <- MASS::glm.nb(y ~ x, data = df)` rather than loading the
namespace implicitly.

**Joining unmarked multi-year data**

When `unmarked` model rows repeat spatial sites, provide a `join_by`
data frame with one row per prepared site. For example, if the spatial
sites have IDs `c("A", "B")` and model rows have site IDs
`c("A", "B", "A", "B")` for two years, use
`join_by = data.frame(site = c("A", "B"))` and include the four-row
`site` column in the model's `siteCovs`. Multiple key columns are
allowed. Prepared keys must be unique and every model row must match.
The model's response and site-covariate row order is preserved. Without
`join_by`, the number and order of model rows must equal the prepared
spatial sites. When both objects retain the same row IDs, their order is
checked. Otherwise, the order must be supplied correctly by the caller.
Optimized covariates must be static `siteCovs`; `obsCovs` and
`yearlySiteCovs` cannot be optimized through this site-level workflow.

**Custom refit functions**

By default, `multiScale_optim()` refits the model using
[`stats::update()`](https://rdrr.io/r/stats/update.html) (standard
models) or the `unmarked` equivalent. When the default path is
insufficient, supply `refit_fn`. The function must accept arguments
`model`, `data`, and `context`, and return a fitted model whose
log-likelihood can be extracted by
[`stats::logLik()`](https://rdrr.io/r/stats/logLik.html) or
[`insight::get_loglikelihood()`](https://easystats.github.io/insight/reference/get_loglikelihood.html).

A minimal custom refit function:


    refit_fn <- function(model, data, context) {
      stats::update(model, data = data)
    }

For models requiring explicit call reconstruction:


    refit_fn <- function(model, data, context) {
      call <- model$call
      call$data <- quote(data)
      eval(call, envir = list(data = data), enclos = parent.frame())
    }

With PSOCK clusters, `refit_fn` must serialize cleanly. Prefer
namespace-qualified calls (e.g.,
[`stats::update()`](https://rdrr.io/r/stats/update.html)) and avoid
closures over large local objects.

**Wrapper model objects**

Some packages return wrapper objects around another fitted model (e.g.,
`amt::fit_clogit()` wraps a
[`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html) fit
in its `model` component). When possible, `multiScale_optim()` unwraps
these automatically. For `amt::fit_clogit()`, pass `model = TRUE` when
fitting so the nested `clogit` model retains its model frame.

When \`start_strategy = "screen"\` and \`par\` is left as \`NULL\`,
\`multiScale_optim()\` first scouts the sigma space with one-dimensional
log-spaced scans, then optionally tests a few jittered starts using
short, Hessian-free optimizations. The marginal sigma scans are serial,
but the short screening attempts use \`n_cores\` when parallel
optimization is requested. These screening steps are used only to choose
one starting vector for the single full optimization. On Windows, those
parallel screening attempts use PSOCK workers, so the same PSOCK
serialization guidance applies as in the full optimization. For
reproducible screened starts, call \`set.seed()\` before
\`multiScale_optim()\`.

## See also

[`kernel_dist`](https://wpeterman.github.io/multiScaleR/reference/kernel_dist.md)

## Examples

``` r
# \donttest{
set.seed(555)

points <- vect(cbind(c(5,7,9,11,13),
                     c(13,11,9,7,5)))

mat_list <- list(r1 = rast(matrix(rnorm(20^2),
                                  nrow = 20)),
                 r2 = rast(matrix(rnorm(20^2),
                                  nrow = 20)))
rast_stack <- rast(mat_list)
kernel_inputs <- kernel_prep(pts = points,
                             raster_stack = rast_stack,
                             max_D = 5,
                             kernel = 'gaussian',
                             sigma = NULL)
#> 
#> No sigma values provided...
#> Creating necessary elements to optimize sigma
#> 
#> Extracting values from raster layers...
## Example response data
y <- rnorm(5)

## Create data frame with raster variables
dat <- data.frame(y = y,
                  kernel_inputs$kernel_dat)
mod1 <- glm(y ~ r1 + r2,
            data = dat)

## NOTE: This code is only for demonstration
## Optimization results will have no meaning
opt_mod <- multiScale_optim(fitted_mod = mod1,
                            kernel_inputs = kernel_inputs,
                            par = NULL,
                            n_cores = NULL)
#> 
#> 
#> Optimization complete
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent.
#>  Consider increasing max_D in `kernel_prep` to >=16.44 to ensure accurate estimation of scale.
#> 
#> 

## Using package data
data('pts')
data('count_data')
hab <- terra::rast(system.file('extdata',
                   'hab.tif', package = 'multiScaleR'))

kernel_inputs <- kernel_prep(pts = pts,
                             raster_stack = hab,
                             max_D = 250,
                             kernel = 'gaussian')
#> 
#> No sigma values provided...
#> Creating necessary elements to optimize sigma
#> 
#> Extracting values from raster layers...

mod <- glm(y ~ hab,
           family = poisson,
           data = count_data)

## Optimize scale
opt <- multiScale_optim(fitted_mod = mod,
                        kernel_inputs = kernel_inputs)
#> 
#> 
#> Optimization complete
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent.
#>  Consider increasing max_D in `kernel_prep` to >=250.56 to ensure accurate estimation of scale.
#> 
#> 

## Summary of fitted model
summary(opt)
#> 
#> Call:
#> multiScale_optim(fitted_mod = mod, kernel_inputs = kernel_inputs)
#> 
#> 
#> Kernel used:
#> gaussian
#> 
#> ***** Optimized Scale of Effect *****
#> 
#>         Mean       SE     2.5%    97.5%
#> hab 76.16637 9.567451 57.09847 95.23428
#> 
#> 
#>   ==================================== 
#> 
#> ***** Optimized Scale of Effect -- Distance *****
#> Scale distances (map units) at 90%, two definitions.
#> 
#> 2D, radius of the circle for the ideal untruncated kernel (selected):
#>       Mean   2.5%  97.5%
#> hab 163.45 122.53 204.37
#> 
#> 1D, along a line through the point (comparison):
#>       Mean  2.5%  97.5%
#> hab 125.28 93.92 156.65
#> 
#> 
#>   ==================================== 
#> 
#>  *****     Fitted Model Summary     *****
#> 
#> 
#> Call:
#> glm(formula = y ~ hab, family = poisson, data = data)
#> 
#> Coefficients:
#>             Estimate Std. Error z value Pr(>|z|)    
#> (Intercept)   0.4807     0.1021   4.707 2.51e-06 ***
#> hab           0.8691     0.0855  10.166  < 2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> (Dispersion parameter for poisson family taken to be 1)
#> 
#>     Null deviance: 186.585  on 74  degrees of freedom
#> Residual deviance:  70.185  on 73  degrees of freedom
#> AIC: 232.29
#> 
#> Number of Fisher Scoring iterations: 5
#> 
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent. 
#>  Consider increasing max_D in `kernel_prep` to >= 250.56 to ensure accurate estimation of scale.
#> 
#> 
#> 
#> ***** Follow-up recommendation *****
#> Recommended `max_D`: 250.56
#> Recommended starting `sigma` values:
#>      hab 
#> 76.16637 
#> 
#> Action:
#> Refit from the optimized parameters and address the flagged diagnostics. 

## 'True' parameter values data were simulated from:
# hab scale = 75
# Intercept = 0.5,
# hab slope estimate = 0.75

## Plot and visualize kernel density
plot(opt)



## Apply optimized kernel to the environmental raster
opt_hab <- kernel_scale.raster(hab, multiScaleR = opt)
#> 
#> Smoothing spatRaster variable 'hab' from layer 'hab' at sigma = 76

plot(c(hab, opt_hab))


## Project model; scale & center
opt_hab.s_c <- kernel_scale.raster(raster_stack = hab,
                                   multiScaleR = opt,
                                   scale_center = TRUE)
#> 
#> Smoothing spatRaster variable 'hab' from layer 'hab' at sigma = 76

mod_pred <- predict(opt_hab.s_c, opt$opt_mod, type = 'response')
plot(mod_pred)


## Custom refit hook for model classes that need explicit control.
## This example still uses glm(), but the same pattern can be used for
## classes whose default update path is not sufficient.
refit_glm <- function(model, data, context) {
  stats::update(model, data = data)
}

opt_custom <- multiScale_optim(fitted_mod = mod,
                               kernel_inputs = kernel_inputs,
                               refit_fn = refit_glm)
#> 
#> 
#> Optimization complete
#> 
#>  WARNING!!!
#>  The estimated scale of effect approaches or exceeds the current search extent.
#>  Consider increasing max_D in `kernel_prep` to >=250.56 to ensure accurate estimation of scale.
#> 
#> 
# }
```
