# Create an Audited Refit Function for Distance-Removal Models

Build a refit function for `multiScale_optim(refit_fn = ...)` that
refits an `unmarked` distance-removal (GDR) model at each candidate
scale of effect from its original estimates and from extra starting
values you declare. It keeps the best clean fit and writes a record of
every attempt to an audit directory. A distance-removal likelihood can
have more than one local optimum, and several starts make it less likely
that a poor optimum is accepted unnoticed.

## Usage

``` r
make_gdr_refitter(
  model,
  cold_starts,
  audit_dir,
  provenance = list(),
  tmb_settings = NULL
)
```

## Arguments

- model:

  The starting model. One of:

  - an `unmarkedFitGDR` fitted by
    [`unmarked::gdistremoval()`](https://ecoverseR.github.io/unmarked/reference/gdistremoval.html)
    with `engine = "C"`;

  - an `unmarkedFitGDR` fitted with `engine = "TMB"`, including models
    with random intercepts;

  - an
    [unmarkedFitGDRField](https://wpeterman.github.io/multiScaleR/reference/unmarkedFitGDRField-class.md)
    from
    [`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md).

  The refit function is tied to this model. It reuses the model's data,
  formulas, cutoff `K`, mixture, distance function, output units and
  optimizer settings, and it never switches engine. The model needs
  `unmarked` 1.5.2, plus `TMB` 1.9.25 for the TMB engine. Formulas with
  [`offset()`](https://rdrr.io/r/stats/offset.html) terms and
  zero-inflated TMB models (`mixture = "ZIP"` with `engine = "TMB"`) are
  not supported.

- cold_starts:

  Named list of one or more extra starting vectors on the link scales of
  the model, for example `list(cold = cold)`. For a C-engine model, name
  each vector exactly like `unmarked::coef(model)`, such as
  `"lambda(Int)"` and `"dist(Int)"`. For a TMB-engine or field model,
  supply one value for every fixed coefficient and log SD, either
  unnamed in the order of `gdr_start_values(model)$starts`, or named by
  parameter block and position, such as `"beta_lambda[1]"` or
  `"lsigma_lambda[1]"`. Random effects are not included. Each vector
  must differ from the model's estimates and from the other vectors, and
  the names `"original"` and `"final"` are reserved. These values are
  model coefficients, not kernel scales of effect and not the cutoff
  `K`.

- audit_dir:

  Character string naming the directory for audit records. It must be
  new or empty, and it is created if needed. (For a field model, an
  existing directory with records is also accepted.) Keep it, because
  the records are the only trace of failed attempts. See "Read the audit
  directory".

- provenance:

  Optional named list of your own identifiers, such as a data source, a
  script version or an archive DOI. Use simple values: atomic vectors of
  at most 64 elements, or nested lists of them, with at most 64 KiB in
  total. It is copied into every record as supplied and is not checked.

- tmb_settings:

  For a TMB-engine or field model, a named list of the optimizer
  settings to use for every refit. Leave it `NULL` (the default) for a
  C-engine model, whose settings are read from the fitted model.
  Required elements:

  - `method`: one of `"BFGS"`, `"Nelder-Mead"`, `"CG"` or `"L-BFGS-B"`.

  - `se`: logical. Whether the final refit of each call computes
    standard errors.

  - `threads`: must be 1.

  - `control`: named list of finite numeric
    [`optim()`](https://rdrr.io/r/stats/optim.html) controls with a
    positive integer `maxit`, such as
    `list(maxit = 2000L, reltol = 1e-10)`.

  Optional elements:

  - `lower`, `upper`: parameter bounds for `method = "L-BFGS-B"`, either
    one number or one value per coordinate in the fitted optimizer order
    (`names(gdr_start_values(model)$fitted_values)`), which can differ
    from the start order.

  - `random.start`: how TMB starts each search for the random-effect
    values. `"cached_best"` (the TMB default) starts from the best
    values found so far. `"zero"` starts from zero.

  - `inner.control`: a list with positive `maxit`, `tol`, `grad.tol` and
    `step.tol` for that inner search.

  - `record_memory`: `TRUE` records the memory used by the R process at
    each phase of a fit. It requires the `ps` package. The default is
    `FALSE`.

  For a field model, pass the identical list used in
  [`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md).

## Value

A function with arguments `model`, `data` and `context`, for the
`refit_fn` argument of
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md).
The optimizer calls it with the starting model, a copy of its frame in
which the kernel-weighted covariates have new values, and a context
whose `covs` element names those columns. Each call returns one fitted
model of the same class as `model`, with an attribute
`msr_gdr_refit_audit` that summarizes the call:

- `status`: `"returned"` for a successful call.

- `attempts`: one record per start (`"original"`, each cold start, and
  `"final"`), each with the `start`, `converged`, the `logLik`,
  `coefficients` on the link scale, `warnings`, any `error`, the
  `covariance` status, whether it was `eligible`, and the optimizer
  output (`par`, `value`, `counts`, `convergence`, `message`). TMB
  attempts add `marginal_score` (largest absolute gradient, its 0.001
  threshold and whether it passed) and `optimizer_point_file`.

- `selected`: the index of the start with the highest log-likelihood.

- `objective_only`: `TRUE` when the call skipped the final refit.

- `contract_sha256` and `input_context_sha256`: fingerprints that tie
  the call to `contract.rds` and to the inputs it received. `path` is
  the summary record.

A field model returns the record described in
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md)
instead. A failed call stops with an error of class
`msr_gdr_refit_error`, whose `audit` and `audit_file` elements hold the
record and its path.

## How each refit works

At every call the function fits the model from its original estimates
(the start named `"original"`) and from each cold start, then keeps the
start with the highest finite log-likelihood. That start must be
eligible: converged (code 0) with no warnings and, for TMB models,
usable covariance and a largest absolute gradient below 0.001. If it is
not eligible, the call fails rather than returning a lower-likelihood
fit in its place.

When
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
needs the fitted model rather than only its likelihood, a final fit
starts from the selected coefficients with the original standard-error
request. It must also be eligible and may lose no more than 1e-5
log-likelihood units. C-engine intermediate fits skip coefficient
Hessians, which saves time.

Only the columns named in `context$covs` may change between calls, and
only static site covariates qualify. A covariate whose name also appears
in the yearly or observation covariates is rejected, so give it a
distinct name. Changes to the counts, the other covariates, the survey
protocol or the row order stop the call.

Each TMB attempt builds a new TMB objective and starts its random
effects at zero. It always runs
[`TMB::sdreport()`](https://rdrr.io/pkg/TMB/man/sdreport.html), even
with `se = FALSE`, so TMB refits are not faster without standard errors.
Usable covariance requires a positive-definite Hessian, fixed covariance
and joint precision. The random effects must sit at a mode with
positive-definite curvature and a largest gradient and Newton decrement
below 0.001. The objective, recomputed at the optimizer's coordinates,
must also agree with the reported value within 1e-7. These are local
checks. They do not show that the mode is unique or that TMB's Laplace
approximation is accurate.

## Declare TMB settings

A C-engine fit keeps the optimizer settings it used, so they are read
automatically. A TMB fit in `unmarked` 1.5.2 does not keep them, and its
constructor does not pass optimizer controls through. Declare them in
`tmb_settings` instead. The refit function passes them to
[`optim()`](https://rdrr.io/r/stats/optim.html) directly, without
changing `unmarked`. A typical list is:

    list(method = "BFGS", se = TRUE, threads = 1L,
         control = list(maxit = 2500L, reltol = 1e-10))

`random.start` and `inner.control` change how TMB searches for the
random effects. A different setting can reach a different numerical
solution, so treat a change as a new analysis and check its results. No
setting is changed automatically after a failure. With
`record_memory = TRUE`, each phase record stores the resident memory of
the whole R process and its lifetime peak, in bytes. These numbers
include everything else loaded in the session, so they are not the
memory used by one phase.

## Read the audit directory

|  |  |  |
|----|----|----|
| File | Written | Contents |
| `contract.rds` | When the refit function is created | Fingerprints of the original data and design, the settings, the starts, the original coefficients and optimizer summary, and the identity of `unmarked`, TMB and `multiScaleR`, with a SHA-256 of the whole contract |
| `refit-*-request.rds` | At the start of each call | The call's opening record |
| `refit-*-input.rds` | After the inputs pass their checks | Fingerprints of the new data and the changed columns |
| `refit-*-attempt-*-request.rds` | Before each start enters `unmarked` | The start and its name |
| `refit-*-attempt-*-point.rds` | TMB only, after optimization and before covariance | The optimizer result |
| `refit-*-attempt-*.rds` | After each start | The completed attempt record |
| `refit-*.rds` | At the end of each call | The summary, also attached to the returned model |

A request without a matching completed record documents an attempt that
never finished, for example after an interruption. Keep it. Records
contain numbers, settings and SHA-256 fingerprints. They contain no raw
data, design matrices, Hessians or fitted models, so they stay small.

The optimizer point file (class `msr_gdr_optimizer_point`) is written
before the covariance step, so a later covariance failure still leaves
evidence that the optimization finished. It holds the optimizer's `par`,
`value`, `counts`, `convergence` and `message`, the coordinate labels,
data and layout fingerprints, the cutoff `K`, the declared settings and
the elapsed time of each phase. It is always marked
`inference_eligible = FALSE`. It is not a fitted model, it gives no
standard errors, and a fit cannot be resumed from it.

At each call the refit function also confirms that `unmarked`, TMB and
this function's own code have not changed since it was created. A
change, including a different installation in a new worker process,
stops the call before fitting.

## Refit a response field

For a model from
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md),
each call refits the field with all of that function's checks, with the
range held fixed. Use the same `cold_starts` and the identical
`tmb_settings`, and run
[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md)
with `n_cores = 1`. The ordinary
[`update()`](https://rspatial.github.io/terra/reference/update.html)
method stops for field models, so this refit function is the way to
refit one.

## What the checks do not establish

Coefficient standard errors from the returned model treat the selected
scales of effect as known. The function adds no priors, bootstrap
intervals or kernel-scale uncertainty. Agreement among the declared
starts does not guarantee that a better optimum does not exist
elsewhere.

The cutoff `K` (the largest abundance the likelihood considers at a
point) is unrelated to the `K` column of
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md),
which counts parameters. Reusing one cutoff does not replace a check
that results are insensitive to it.

See the "Refit from several starting values" section of
[`vignette("v12-distance-removal", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v12-distance-removal.md)
for a worked example.

## See also

[`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
[`gdr_start_values()`](https://wpeterman.github.io/multiScaleR/reference/gdr_start_values.md),
[`fit_gdr_response_field()`](https://wpeterman.github.io/multiScaleR/reference/fit_gdr_response_field.md),
[`validate_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/validate_unmarked_data.md),
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md)

## Examples

``` r
if (as.character(utils::packageVersion("unmarked")) == "1.5.2") {
  # Five point counts: birds by distance band and by time interval
  yd <- rbind(c(0, 0), c(1, 0), c(0, 2), c(2, 1), c(0, 0))
  yr <- rbind(c(0, 0, 0), c(0, 1, 0), c(1, 0, 1), c(1, 1, 1), c(0, 0, 0))
  frame <- unmarked::unmarkedFrameGDR(yd, yr, dist.breaks = c(0, 25, 75),
    unitsIn = "m", period.lengths = c(1, 2, 1),
    siteCovs = data.frame(hab = seq(-1, 1, length.out = 5)))
  # Fit the starting model (C engine) before declaring extra starts
  model <- unmarked::gdistremoval(lambdaformula = ~hab, data = frame,
    output = "abund", K = 40, se = FALSE, engine = "C",
    starts = c(log(8), 0, log(40), qlogis(0.3)), control = list(maxit = 2000))
  # One extra start, named like the model coefficients (link scales)
  cold <- unmarked::coef(model)
  cold[] <- 0
  cold["dist(Int)"] <- log(30)
  audit_dir <- tempfile("msr-gdr-audit-")
  refit <- make_gdr_refitter(model, list(cold = cold), audit_dir = audit_dir)

  # multiScale_optim() calls refit() at each candidate scale. Calling it
  # once here, on the original data, shows what each call records.
  fit <- refit(model, frame, context = list(covs = "hab"))
  audit <- attr(fit, "msr_gdr_refit_audit")
  # Log-likelihood and eligibility of each start, then the selected start
  print(sapply(audit$attempts, function(a)
    c(logLik = a$logLik, eligible = a$eligible)))
  print(audit$attempts[[audit$selected]]$name)
  # contract.rds, plus request, input, attempt and summary records
  print(list.files(audit_dir))
}
#>               [,1]      [,2]      [,3]
#> logLik   -14.77773 -14.77773 -14.77773
#> eligible   1.00000   1.00000   1.00000
#> [1] "original"
#>  [1] "contract.rds"                                           
#>  [2] "refit-7970-20261117221256687209-1-attempt-1-request.rds"
#>  [3] "refit-7970-20261117221256687209-1-attempt-1.rds"        
#>  [4] "refit-7970-20261117221256687209-1-attempt-2-request.rds"
#>  [5] "refit-7970-20261117221256687209-1-attempt-2.rds"        
#>  [6] "refit-7970-20261117221256687209-1-attempt-3-request.rds"
#>  [7] "refit-7970-20261117221256687209-1-attempt-3.rds"        
#>  [8] "refit-7970-20261117221256687209-1-input.rds"            
#>  [9] "refit-7970-20261117221256687209-1-request.rds"          
#> [10] "refit-7970-20261117221256687209-1.rds"                  
```
