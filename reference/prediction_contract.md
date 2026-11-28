# Read prediction contracts

A prediction contract is a list of metadata that travels with a
prediction. It names the quantity that was predicted, its scale, how
random effects were treated, and which uncertainty the standard errors
or intervals include and leave out. Read it when you are unsure what a
result describes. A contract describes a calculation. It does not
certify that a fit is eligible for prediction or that an interval is
calibrated.
[multiScaleR-prediction](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-prediction.md)
explains the quantities and uncertainties that the fields name.

## Find the contract

|  |  |
|----|----|
| Result | Where the contract is |
| Snapshot from [`prepare_prediction()`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md) | `snapshot$prediction_contract`, a named list with one contract per prepared output type, such as `$response` and `$link` |
| Native refit prediction from [`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md) | `attr(result, "prediction_contract")` |
| Raster prediction from [`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md) | `attr(map, "prediction_contract")`, beside the `projection_contract` attribute |
| Effect curves from [`plot_marginal_effects()`](https://wpeterman.github.io/multiScaleR/reference/plot_marginal_effects.md) for `unmarked`, glmmTMB, lme4, multi-component GAM, and Bayesian fits | `attr(curves, "prediction_contract")` |
| Bayesian summaries and simulations from [`predict.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR_bayes.md), [`simulate.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/simulate.multiScaleR_bayes.md), and [`predict_scales()`](https://wpeterman.github.io/multiScaleR/reference/predict_scales.md) | `attr(result, "prediction_contract")` |

Point predictions from the analytic backends return plain vectors or
lists without a contract. Prepare a snapshot to get one. Raw posterior
draw matrices and lists keep their existing attributes, and native
numerical outputs inside a snapshot keep their native attributes.

## Read the fields

- `schema_version`: version of the metadata structure.

- `inference_target`, `method`: the posterior, likelihood, conditional
  likelihood, or native target, and the label of the computation that
  produced it.

- `quantity`, `output_scale`: what the numbers are, on the stated scale.
  The quantity is an expected response, a linear predictor, a kernel
  scale, a contrast between kernel scales, a replicate response, a
  choice probability, or a named native component. Native component
  meanings remain model specific.

- `target_population`: the supplied prediction rows, or the specified
  location pairs for private spatial scale contrasts.

- `effects`: the policy for known and unseen groups. Zero deviations
  (`"zero_deviation_reference"`) are reference values, not integration
  over the group distribution. Sampled new effects are shared by each
  group's rows within a complete joint vector.

- `conditioning`, `integration`: the quantities held fixed and the
  recorded treatment of uncertainty. Missing historical solver metadata
  is marked as unrecorded rather than inferred.

- `uncertainty`: the treatment of the population and kernel parameters,
  whether observation noise is included (`observation_noise`), what the
  intervals mean, what is omitted (`omitted`), and `calibration`, which
  is always `"not_established"`. Analytic snapshots also record whether
  joint covariance was available at preparation (`joint_available`), and
  `joint_conditional_on` names any kernel parameters at a search bound
  that it holds fixed. This does not bypass the covariance checks in
  [`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md).

- `original_coordinates`, `joint_parameter_map_complete`: the parameter
  representation used by the adapter. For Bayesian predictions,
  completeness checks the stored block dimensions against the original
  prepared parameter map. This structural check does not validate a
  coordinate transformation, reconstruct integrated coefficients, or
  establish posterior accuracy.

- `diagnostics_status`: the existing fit status when available,
  `"passed"`, `"failed"`, or `"unavailable"`. Joint covariance
  availability is separate from this status.

- `kernel_support`: for native predictions, any recoverable table of
  optimization bounds and the existing extraction-window diagnostic.
  These flags qualify the selected kernel. They do not add kernel
  uncertainty, alter conditional coefficient standard errors, or
  establish prediction eligibility.

- `metadata_establishes_eligibility`: always `FALSE`. Read the original
  fit diagnostics before you interpret predictions.

- `fixed_inputs`, `component`: whether the rows and controls are bound
  to a snapshot, and the named component when one applies.

## Read the uncertainty record

Expected-response intervals exclude observation noise, and replicates
include it. New-group predictions also carry the group variation you
requested. Native confidence intervals and standard errors condition on
the selected kernel parameters. Likelihood summaries conditional on one
population parameter vector omit its estimation uncertainty.

- Native point and raster contracts record `uncertainty$standard_errors`
  as `"provided"` or `"not_requested"`, and separately keep
  `se.fit_requested`. Provider controls such as glmmTMB's `cov.fit` can
  return standard errors even when the wrapper's `se.fit` is `FALSE`.
  Returned native covariance is conditional on the selected kernel
  parameters.

- spaMM variance attributes keep their native meaning and are listed in
  `native_variances`, with their numerical attributes unchanged.
  `fixefVar` and `predVar` are on the linear-predictor scale. `residVar`
  and `respVar` are on the response scale, and the response variance
  includes observation noise, which is distinct from uncertainty in the
  mean.

- Native interval attributes are listed in `native_intervals` and keep
  the provider's target. Native `lm` interval columns record
  `"native_mean_confidence"` or `"native_prediction_interval"` and the
  bound `interval_level`. Only the prediction interval includes
  future-observation noise. Both include native parameter uncertainty
  and hold the selected kernel fixed.

- Native [`mgcv::ziplss()`](https://rdrr.io/pkg/mgcv/man/ziplss.html)
  response means combine both linear predictors. Their standard errors
  propagate the full native coefficient covariance through that mean,
  keeping signed cross-covariance. The contract records this as
  `response_standard_error_method`, distinct from the provider's own
  response standard-error formula. Kernel uncertainty and count noise
  remain omitted, and native smoothing-covariance requests keep their
  meaning.

- Refit snapshot contracts describe their cached point output, and
  `conditional_se_available` records whether a cached standard-error
  request succeeded. `unmarked` snapshots describe the retained native
  confidence tables. Mixed native snapshots keep both point and standard
  error choices, and `"selection_by_request"` identifies that
  capability. None of these fields certifies finite standard errors,
  covariance eligibility, or calibration.

## Random effects in unmarked contracts

`unmarked` treats any non-`NULL` `re.form` as zero random effects. The
contract records zero deviations and no integration for that request,
including formula controls that other providers use to select some
effects. Native argument abbreviations are resolved before this policy
is recorded. Zero-effect contexts do not evaluate unused random
transformations or add latent variation. Conditional requests still need
their original random-design metadata.

## Choice sets and selection scores

Conditional logistic linear predictors and exponentiated relative
selection scores need no prediction choice sets. Choice probabilities
need complete sets. A snapshot can keep the reason probabilities are
unavailable while its link and risk outputs remain available. Relative
scores are not occurrence probabilities or utilization distributions.

## Spatial response fields

Native glmmTMB response and link predictions with recognized
exponential, Gaussian, or Matern response covariance add
`response_spatial_field`. It records the provider and layout status, the
conditional-response component, and each field's native term,
covariance, coordinate and grouping variable names, number of locations,
independent fields, and covariance-parameter count. Conditioning records
either the native spatial prediction or all response random effects set
to zero. New-location controls keep their native meaning. Coordinate
units are not inferred from factor labels. The field range is distinct
from the ecological kernel scale, and kernel uncertainty is not
included. Identification checks the provider's registered covariance
codes and the fitted block layout. Coordinate-aware field maps also
record `geometry`: the projected coordinate reference system and units,
native column names and shared group, reference and context hashes,
serial batch size, and memory budget. Cell centers supply the
coordinates, and the working-memory estimate is not an operating-system
memory cap.

Identified native `gdistremoval()` response fields record a fixed
correlation range, conditional physical abundance effects, and a
distinct residual variance at new locations. Their standard errors keep
the full fitted fixed and latent covariance. Density point predictions
are inverse-link conditional predictors, not lognormal field averages.
These native field calculations do not enable the private RTMB fitting
route, field simulation, spatial cross-validation, or ecological-scale
uncertainty.

Spatial-field cross-validation records a joint assessment-response
component. It uses exact Gaussian conditioning for Gaussian responses
and the Laplace marginal likelihood ratio for other responses. Both
include observation noise and conditional field variation at the fitted
training parameters. Log-link response moments use exponentially tilted
Laplace ratios, and Bernoulli moments use a Gaussian approximation. The
predictive standard deviation describes observation variation. For
non-Gaussian responses, the contract identifies field variance on the
linear-predictor scale and observation variance on the response scale.

## Kernel-scale contrasts

Private spatial contrasts use `quantity = "kernel_scale_contrast"` and
`output_scale = "log_ratio"` or `"ratio"`. Both locations use the same
complete joint vector. Their intervals summarize the contrasts rather
than differences between marginal interval endpoints. Observation noise
is omitted, and calibration remains unestablished.

## Subsets keep the original contract

A contract describes the original prediction call. Ordinary data-frame
subsetting can keep the attribute, but it does not recompute the
contract for a smaller batch. Use the existing row-level new-group flags
for individual rows. No contract establishes interval calibration or
upgrades a rejected fit.

## See also

[multiScaleR-prediction](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-prediction.md),
[`prepare_prediction()`](https://wpeterman.github.io/multiScaleR/reference/prepare_prediction.md),
[`predict.multiScaleR()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR.md),
[`predict.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR_bayes.md),
[`predict_scales()`](https://wpeterman.github.io/multiScaleR/reference/predict_scales.md),
[`simulate.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/simulate.multiScaleR_bayes.md),
and
[`vignette("v09-estimates-and-intervals", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v09-estimates-and-intervals.md).

## Examples

``` r
## A stored N-mixture fit used in vignette("v02-user-guide")
nmix <- readRDS(system.file("extdata", "cache-user-guide-nmix.rds",
                            package = "multiScaleR"))

## A snapshot stores one contract per prepared component
snap <- prepare_prediction(nmix, component = "state")
contract <- snap$prediction_contract$state
contract$quantity                       # a named native component
#> [1] "native_component"
contract$uncertainty$kernel_parameters  # "fixed_at_estimate": sigma held fixed
#> [1] "fixed_at_estimate"
contract$uncertainty$omitted            # what the limits leave out
#> [1] "kernel_parameter_uncertainty" "observation_noise"           

## A native refit prediction carries its contract as an attribute
abundance <- predict(nmix, type = "state", se.fit = TRUE,
                     uncertainty = "conditional")
head(abundance)                         # expected abundance with 95% limits
#>   Predicted         SE     lower     upper
#> 1 0.8899691 0.19558543 0.5785079 1.3691169
#> 2 0.6609388 0.10951077 0.4776678 0.9145272
#> 3 2.5001106 0.27095126 2.0216671 3.0917815
#> 4 3.4768612 0.43972555 2.7135273 4.4549261
#> 5 0.9695994 0.14269679 0.7266416 1.2937918
#> 6 0.3752776 0.07535654 0.2531801 0.5562573
attr(abundance, "prediction_contract")$uncertainty$intervals
#> [1] "native_conditional_confidence"
```
