# Prepare a fixed set of predictions for repeated use

Compute the predictions for one fixed set of rows once and keep them in
a compact snapshot. For the analytic and conditional logistic backends,
the snapshot holds the prepared spatial design, so
[`predict()`](https://rspatial.github.io/terra/reference/predict.html)
on it can switch between scales and uncertainty choices without
rebuilding the kernels. For joint mixed, `unmarked`, and other native
refit models, it holds the native numerical outputs for the controls you
bind. A snapshot does not add uncertainty, update the fit, or predict
new locations.

## Usage

``` r
prepare_prediction(
  object,
  newdata = NULL,
  kernel_inputs = NULL,
  max_bytes = 64 * 1024^2,
  component = NULL,
  native_args = list(),
  join_by = NULL,
  include_fixed = FALSE
)

# S3 method for class 'multiScaleR_prediction'
predict(
  object,
  type = c("response", "link"),
  se.fit = FALSE,
  uncertainty = c("joint", "conditional"),
  interval = c("none", "confidence"),
  level = 0.95,
  ...
)

# S3 method for class 'multiScaleR_prediction'
print(x, ...)
```

## Arguments

- object:

  A fitted `multiScaleR` object from
  [`multiScale_optim()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_optim.md),
  using `backend = "profile_glm"`, `"joint_nb2"`, `"conditional_logit"`,
  or `"joint_mixed"`, or the default refit backend with a supported
  native model in `object$opt_mod`. For
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html),
  a snapshot returned by `prepare_prediction()`.

- newdata, kernel_inputs:

  Prediction rows and matching spatial inputs from
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md),
  with the same requirements as in
  [`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md).
  The default `NULL` for both uses the rows retained from fitting.

- max_bytes:

  Positive number. The largest allowed size of the returned snapshot, in
  bytes. The default, `64 * 1024^2`, is 64 MiB. A larger snapshot raises
  an error, so use fewer rows or raise the limit. It limits the stored
  size, not the temporary memory used while preparing.

- component:

  Character scalar naming the native prediction type to keep.

  - `unmarked` models: required. Name one component, such as `"state"`,
    `"det"`, `"psi"`, `"col"`, or `"ext"`, as the fitted model supports.

  - `"joint_mixed"`: optional. Name one native glmmTMB type, such as
    `"zprob"`. The default `NULL` prepares both `"response"` and
    `"link"`.

  - Other refits: the default `NULL` prepares `"response"` and `"link"`
    for `lm`, `glm`, glmmTMB, lme4, `gls`, and `lme` models, and `"lp"`
    and `"risk"` for Cox models. Other native providers need an explicit
    component.

  - Analytic and conditional logistic backends: leave `NULL`. Their
    snapshots keep every type (`"response"`, `"link"`, and, for
    conditional logistic fits, `"risk"`).

- native_args:

  A named list of native prediction controls to bind to the snapshot.
  Default [`list()`](https://rdrr.io/r/base/list.html), which keeps the
  native defaults. Changing any control later needs a new snapshot.
  Supported controls:

  - `"joint_mixed"`: `re.form`, `allow.new.levels`, and `na.action`.

  - `unmarked`: `backTransform`, `na.rm`, `appendData`, `level`, and
    `re.form`.

  - Other refits: `re.form`, `allow.new.levels`, `na.action`,
    `reference`, and `level` (such as nlme `level = 0`), when the native
    method accepts them. Ordinary `lm` refits also accept `interval`
    (`"confidence"` for mean limits, `"prediction"` for
    future-observation limits). spaMM refits accept `variances` and
    `intervals`, but cross-row covariance needs point predictions
    instead of a snapshot. mgcv GAMs accept `exclude`, `terms`,
    `unconditional`, `iterms.type`, and `block.size`, and BAMs also
    accept `discrete`. `unconditional = TRUE` requests the native
    smoothing-parameter covariance, which adds no kernel-scale
    uncertainty.

  Native argument abbreviations are resolved to their full names before
  the list is checked and are recorded by those names. Selected lme4 and
  spaMM effects and GAM terms are recorded in the prediction contracts.

- join_by:

  For native refit and `unmarked` spatial predictions, a data frame of
  site keys with one row per prepared site, in `kernel_inputs` order. It
  needs both `newdata` and `kernel_inputs`. `unmarked` spatial snapshots
  require it (see
  [`prepare_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/prepare_unmarked_data.md)).
  The default `NULL` aligns general refit rows by row name, and
  `unmarked` data without spatial inputs pass unchanged. Not used by the
  analytic, conditional logistic, and joint mixed backends.

- include_fixed:

  Logical. The default `FALSE` keeps fixed-radius predictors as explicit
  inputs in `newdata`. Set `TRUE` to recompute them from the spatial
  inputs after verifying their training preparation, as described in
  [`prepare_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/prepare_unmarked_data.md)
  and
  [`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md).
  It needs spatial inputs and verified row alignment, with explicit site
  keys for repeated sites and for `unmarked` models.

- type, se.fit, uncertainty, interval, level:

  Output and uncertainty controls for
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  on a snapshot, with the defaults and checks of
  [`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md).
  For analytic and conditional logistic snapshots they choose the scale
  and the uncertainty summary of the fixed design. For joint mixed and
  refit snapshots, `type` selects one prepared component (the first one
  prepared by default), `se.fit` selects the cached point or
  standard-error output, and `uncertainty` must be `"conditional"` when
  standard errors or native interval attributes are returned. `unmarked`
  snapshots accept only `type`.

- ...:

  For
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html),
  any additional argument, including new prediction data, raises an
  error. Prepare a new snapshot for new rows. Ignored by
  [`print()`](https://rdrr.io/r/base/print.html).

- x:

  A `multiScaleR_prediction` snapshot from `prepare_prediction()`.

## Value

`prepare_prediction()` returns an object of class
`multiScaleR_prediction`, printed as a one-line summary of the backend,
the number of fixed predictions, and the stored size. Its elements are
internal, except `prediction_contract`, a named list with one
[prediction_contract](https://wpeterman.github.io/multiScaleR/reference/prediction_contract.md)
per prepared output type. Each records the quantity, its scale and
conditioning, the random-effect policy, and the uncertainty
restrictions. This metadata does not establish fit eligibility or
interval calibration.

[`predict()`](https://rspatial.github.io/terra/reference/predict.html)
on a snapshot returns the same output as
[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md)
would return for the original model and fixed rows. `unmarked` snapshots
return the native table with `Predicted`, `SE`, `lower`, and `upper`
columns. Joint mixed snapshots return native glmmTMB vectors or
standard-error lists for the bound type and controls. Numerical outputs
keep their usual structure, including native attributes.

[`print()`](https://rdrr.io/r/base/print.html) returns the snapshot
invisibly.

## What a snapshot keeps

A snapshot binds the coefficients, covariance, selected scales, training
standardization, row order, choice sets, and native controls of one
call. It keeps neither raw spatial buffers nor the fitted native model,
so it can be saved with
[`saveRDS()`](https://rspatial.github.io/terra/reference/serialize.html)
and restored with
[`readRDS()`](https://rspatial.github.io/terra/reference/serialize.html).
Later changes to the original model or data do not update it. Call
`prepare_prediction()` again after changing the model, locations, choice
sets, metric definitions, distance units, or native controls. New data
cannot be passed to
[`predict()`](https://rspatial.github.io/terra/reference/predict.html)
on a snapshot.

Supported native and analytic models keep their trained standard
scaling, polynomial, and spline calls, including those inside nested
expressions and formula offsets. Recovery needs one unambiguous
learner's attributes per evaluated predictor. Missing or ambiguous
metadata raises an error rather than learning a transformation from the
prediction rows. Ordinary native `lm`, `glm`, and time-fixed unpenalized
Cox snapshots use the original fitted terms for dot shortcuts and active
variables, so removed variables add no prediction inputs and new data
columns never enter the design.

## Analytic and conditional logistic snapshots

These snapshots keep the spatial design, so
[`predict()`](https://rspatial.github.io/terra/reference/predict.html)
on them accepts the same `type`, `se.fit`, `uncertainty`, `interval`,
and `level` choices as
[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md),
with the same restrictions. Reuse avoids rebuilding the spatial engine
when you need several summaries of the same rows. Conditional logistic
snapshots without valid choice-set IDs still give linear predictors and
relative scores, but they reject response probabilities and keep the
reason in their response contract. Their standard errors need
`uncertainty = "conditional"`, and no scale uncertainty is added.

## Native snapshots

Joint mixed, `unmarked`, and other refit snapshots cache the native
numerical outputs for the bound controls, rather than keeping a native
model or rebuilding its prediction machinery. Native warnings are
recorded and repeated each time the cached output is requested. A native
error in the point predictions stops preparation.

- **Joint mixed models.** Point and standard-error outputs are both
  cached. Standard errors stay conditional on the selected scales. Zero
  random-effect controls use formula syntax independently of its
  environment, and grouping columns may be omitted when the provider
  permits zero-effect predictions. glmmTMB `re.form = NA` sets random
  effects to zero and does not integrate over their distribution.

- **General refits.** If the provider cannot supply standard errors, the
  point predictions remain available and a later standard-error request
  repeats the recorded error. Native prediction intervals and response
  variances keep the provider's observation noise. Cached native
  intervals or variance attributes need `uncertainty = "conditional"`,
  even with `se.fit = FALSE`. Leave the wrapper's `interval = "none"`
  when you retrieve those native outputs. GAM rows keep every
  formula-input requirement, including inputs to excluded terms, and the
  unchecked `newdata.guaranteed` shortcut is not supported. spaMM refits
  use their retained native frames to verify fixed features. An
  unsupported native standard-error request stays an explicit error,
  recorded separately from the point output.

- **`unmarked` models.**
  [`predict()`](https://rspatial.github.io/terra/reference/predict.html)
  returns the cached native table. Its confidence level and
  back-transformation are fixed by `native_args` and cannot be changed
  at prediction time.

## Unmarked snapshots

Without `kernel_inputs`, supply native prediction rows with spatial
covariates already on their fitted scale, including raw units when the
native refit used them. With `kernel_inputs`, supply explicit `join_by`
keys and native prediction rows.
[`prepare_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/prepare_unmarked_data.md)
then computes the spatial columns and applies the training
standardization without changing native observation or year ordering or
non-spatial factor and offset inputs. Spatial snapshots keep the
prediction row names and site keys in their provenance.

Ordinary occupancy and `pcount` state and detection terms use the
original raw training transformations, fitted factor contrasts, and
native component covariance. Missing contrasts can be recovered from
unambiguous original main-factor columns. Other missing coding needs
explicit training contrasts. Preparation verifies the fixed designs and
offsets against the retained native fit. Explicit data-frame rows
prepare only the requested component, while native frame inputs or an
omitted `newdata` keep both components for native missing-row removal.
ZIP outputs keep their native zero-inflation estimate and uncertainty,
and every output keeps the full model's parameter count and diagnostics.
Keep the original covariates and fitted designs. Frame subsets keep the
original known-group columns.

`gdistremoval()` snapshots likewise freeze the active `"lambda"`,
`"phi"`, `"dist"`, or `"rem"` component's original training terms and
native covariance, with distinct site, primary-period, and
removal-interval recipes. Data-frame rows limit preparation to the
requested component, and native GDR frames keep every active design and
their original missing-input checks. No model or executable TMB tape is
saved in the snapshot. GDR models with active formula offsets or TMB
zero inflation are rejected before fitting, model selection, or
prediction, as
[`aic_tab()`](https://wpeterman.github.io/multiScaleR/reference/aic_tab.md)
explains.

Conditional learned random-effect terms need precomputed training
predictors. Zero-effect contexts need only their fixed design and
offsets, and they reject a later conditional-effect request. Any
non-`NULL` native `re.form` sets the effects to zero. It neither selects
some effects nor integrates over them. For dynamic occupancy, `"psi"`
describes initial occupancy, and `"col"` and `"ext"` describe
transitions. These component predictions are distinct from smoothed or
projected occupancy trajectories.

## See also

[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md),
[prediction_contract](https://wpeterman.github.io/multiScaleR/reference/prediction_contract.md),
[multiScaleR-prediction](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-prediction.md),
[`prepare_unmarked_data()`](https://wpeterman.github.io/multiScaleR/reference/prepare_unmarked_data.md)

## Examples

``` r
## Analytic fit from the running example: counts and two landscape layers
data("landscape_counts")
data("surv_pts")
pts <- terra::vect(surv_pts)
land <- terra::rast(system.file("extdata", "landscape.tif",
                                package = "multiScaleR"))[[c("land1", "land2")]]
inputs <- kernel_prep(pts, land, max_D = 1700, kernel = "gaussian",
                      bin = TRUE, verbose = FALSE)
dat <- data.frame(landscape_counts, inputs$kernel_dat)
start <- glm(counts ~ site + land1 + land2, family = poisson(), data = dat)
fit <- multiScale_optim(start, inputs, backend = "profile_glm",
                        n_cores = 1, verbose = FALSE)

## Prepare three new survey points once. Row names of new_sites match the
## prepared point IDs.
new_inputs <- kernel_prep(pts[1:3], land, max_D = 1700, kernel = "gaussian",
                          bin = TRUE, verbose = FALSE)
new_sites <- data.frame(site = c(-1, 0, 1))
prepared <- prepare_prediction(fit, new_sites, new_inputs)
prepared
#> Prepared multiscale predictions (profile_glm): 3 fixed predictions, 87.5 Kb.

## Several summaries of the same rows, without rebuilding the kernels
predict(prepared)                                  # expected counts
#>         1         2         3 
#> 0.6399094 1.8208897 5.2991735 
predict(prepared, type = "link", se.fit = TRUE)    # log scale, joint SEs
#> $fit
#>          1          2          3 
#> -0.4464287  0.5993252  1.6675509 
#> 
#> $se.fit
#>         1         2         3 
#> 0.1998652 0.1898499 0.1614594 
#> 
#> $type
#> [1] "link"
#> 
#> $uncertainty
#> [1] "joint"
#> 
predict(prepared, interval = "confidence",
        uncertainty = "conditional")$upper         # upper 95% limits
#>         1         2         3 
#> 0.9398905 2.1443647 6.8517965 

## A native unmarked snapshot binds one component. This stored N-mixture
## fit is used in vignette("v02-user-guide").
nmix <- readRDS(system.file("extdata", "cache-user-guide-nmix.rds",
                            package = "multiScaleR"))
detection <- prepare_prediction(nmix, component = "det")
head(predict(detection))   # detection probability per survey, 95% limits
#>   Predicted         SE     lower    upper
#> 1 0.5106106 0.03446574 0.4432754 0.577563
#> 2 0.5106106 0.03446574 0.4432754 0.577563
#> 3 0.5106106 0.03446574 0.4432754 0.577563
#> 4 0.5106106 0.03446574 0.4432754 0.577563
#> 5 0.5106106 0.03446574 0.4432754 0.577563
#> 6 0.5106106 0.03446574 0.4432754 0.577563
```
