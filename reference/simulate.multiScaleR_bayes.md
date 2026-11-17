# Simulate Posterior Predictive Responses

Generate replicate responses for fitted or new observations. Each
replicate includes joint parameter uncertainty, any requested new-group
effects, and the scatter of individual observations around their means,
so replicates show the values you might record.

## Usage

``` r
# S3 method for class 'multiScaleR_bayes'
simulate(
  object,
  nsim = 1L,
  seed = NULL,
  newdata = NULL,
  kernel_inputs = NULL,
  row_id = NULL,
  re.form = NULL,
  scale_re.form = NULL,
  new_scale_levels = c("error", "sample", "zero"),
  new_response_levels = c("error", "sample", "zero"),
  trials = NULL,
  include_fixed = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  object with posterior draws, from sampling or either approximation.

- nsim:

  Positive integer number of replicate data sets; default 1.

- seed:

  Optional positive integer simulation seed. Your own random-number
  state is restored after simulation, even when it fails.

- newdata:

  `NULL` (default) uses the fitted observations. Otherwise a nonempty
  data frame of new observations that holds the fixed covariates, the
  offset variables, and the grouping variables of every component you
  include. Response values are neither needed nor used.

- kernel_inputs:

  Prepared landscape inputs for `newdata`, from
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  or
  [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md).
  Required with new data and `NULL` otherwise. Use the fitted kernel,
  predictor definitions, and physical map units. See "Match new rows to
  their landscapes" in
  [`predict.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR_bayes.md).

- row_id:

  Name of the column in `newdata` that matches each row to an entry of
  `kernel_inputs$binned$point_ids`, or to a row name of
  `kernel_inputs$kernel_dat` for surface-only inputs. IDs may repeat.
  The default `NULL` uses the fitted row-ID column when the fit kept it,
  and the row names of `newdata` otherwise. Supply it for fits saved by
  older versions of the package. Used only with new data.

- re.form:

  `NULL` (default) includes the response random effects. `NA` sets them
  all to zero while keeping group-specific scales. This is a reference
  prediction, not an average over the random effects.

- scale_re.form:

  `NULL` (default) includes the fitted scale random effects and any
  sampled new ones. `NA` sets every scale deviation to zero. It is
  separate from `re.form`, which controls response effects.

- new_scale_levels, new_response_levels:

  What to do with a group level that was not in the fitted data, set
  separately for scale and response random effects. `"error"` (default)
  stops, `"sample"` draws a new effect, and `"zero"` sets the new effect
  to zero. A sampled effect is one standard-normal value per new group
  and posterior draw, multiplied by that draw's group SD, and shared by
  every row in the group. Zeroing gives a typical-group reference, not a
  population average. Known groups keep their fitted effects. Neither
  choice creates a coefficient for a new level of a fixed factor.

- trials:

  Positive integer binomial trial totals, one value or one per
  observation row. The default `NULL` keeps the fitted totals when you
  simulate the fitted rows. For new observations it uses one trial for a
  binary fit, and it is an error for a grouped binomial fit, which needs
  explicit totals. Use only with binomial models. Response columns in
  `newdata` never supply the totals.

- include_fixed:

  Logical. With new data, `TRUE` recomputes declared fixed-radius
  landscape predictors from `kernel_inputs`, using their original
  preparation means and SDs. The default `FALSE` requires those columns
  in `newdata`. Source, metric, radius, kernel, and resolution must
  match the fitted preparation, including its `max_D` buffer extent.
  This adds no scale parameters or uncertainty.

- ...:

  Not used. Supplying any further argument is an error.

## Value

A data frame with one row per fitted or new observation and one column
per replicate data set, named `sim_1`, `sim_2`, and so on. Values are
counts for Poisson and NB2 models, success counts out of the fitted or
supplied trial totals for binomial models (not proportions), and
response units for Gaussian models. Each column uses one complete
posterior vector, drawn at random with replacement, together with the
requested random effects.

For new data, a `point_ids` attribute gives the prepared landscape ID of
each row, and a `new_groups` attribute lists the unseen levels that were
sampled or set to zero, as in
[`predict()`](https://rspatial.github.io/terra/reference/predict.html).
A `prediction_contract` attribute (see
[prediction_contract](https://wpeterman.github.io/multiScaleR/reference/prediction_contract.md))
marks the result as posterior predictive replicates that include
observation noise and records the policies for included effects. It does
not establish diagnostic eligibility or calibration.

## Check the model against the data

Compare summaries of the replicate data sets with the same summaries of
the observed response. An observed value far outside the replicated
range means the model cannot reproduce that feature. Good sampling
diagnostics do not show that the observation model fits the data, and
agreement on one statistic does not establish adequacy.
[`vignette("v13-bayesian-scales", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v13-bayesian-scales.md)
works through such a check. To use bayesplot, convert the replicates to
a matrix with one simulated data set per row,
`t(as.matrix(replicated))`. For a count model, `ppc_stat()` with
`stat = function(y) mean(y == 0)` checks the fraction of zeros.

## Simulate new groups jointly

Within one call, new-group effects are drawn afresh for each replicate
data set and shared by all of its rows in that group. Different new
groups are independent given the population parameters, whose joint
posterior dependence is kept. Observation noise is then drawn separately
for each row. Separate calls do not share a realization across
components or batches, so put every row that enters a joint question
into one call. Reusing a seed across
[`predict_scales()`](https://wpeterman.github.io/multiScaleR/reference/predict_scales.md),
[`predict()`](https://rspatial.github.io/terra/reference/predict.html),
and [`simulate()`](https://rdrr.io/r/stats/simulate.html) does not
promise matching latent effects, because the three functions use their
random numbers differently.

To get an interval for a total, sum each replicate column first with
[`colSums()`](https://rspatial.github.io/terra/reference/rowSums.html)
and then take quantiles of the sums. Adding the limits of individual
intervals overstates the range of the total.

The three prediction functions return different orientations, which
matters for joint contrasts and totals.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) returns
observations by replicates. `predict(summary = FALSE)` returns posterior
draws by observations, and `predict_scales(summary = FALSE)` returns a
named list of draws-by-observations matrices, one per landscape
predictor.

## Supply trials for binomial data

For a binary fit, new observations default to one trial. For a fit to
`cbind(successes, failures)`, supply the actual future totals, for
example `trials = rep(10L, nrow(newdata))`.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) then returns
success counts out of those totals, while
[`predict()`](https://rspatial.github.io/terra/reference/predict.html)
returns success probabilities whatever the number of trials.

## See also

[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
[`simulate_bayes_prior()`](https://wpeterman.github.io/multiScaleR/reference/simulate_bayes_prior.md),
[`predict.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/predict.multiScaleR_bayes.md),
[`predict_scales()`](https://wpeterman.github.io/multiScaleR/reference/predict_scales.md),
[prediction_contract](https://wpeterman.github.io/multiScaleR/reference/prediction_contract.md),
[multiScaleR_bayes](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md);
[`vignette("v14-bayesian-prediction", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v14-bayesian-prediction.md).

## Examples

``` r
example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
  package = "multiScaleR"))
replicated <- simulate(example$fits$pooled, nsim = 20, seed = 94)
dim(replicated)  # one row per original observation, one column per replicate
#> [1] 100  20

# Posterior predictive check of the fraction of zero counts
range(vapply(replicated, function(y) mean(y == 0), numeric(1)))
#> [1] 0.00 0.04
mean(example$data$count == 0)  # observed fraction of zero counts
#> [1] 0.04
# These replicate counts include observation noise. Fitted means do not.

# Counts at eight future surveys, and an interval for their total
forecast <- readRDS(system.file("extdata", "bayesian-prediction-example.rds",
  package = "multiScaleR"))
future <- simulate(forecast$fit, nsim = 200, seed = 803,
  newdata = forecast$newdata, kernel_inputs = forecast$kernel_inputs,
  row_id = "site_year", new_scale_levels = "sample",
  new_response_levels = "sample")
quantile(colSums(future), c(0.025, 0.5, 0.975))  # total of the eight surveys
#>    2.5%     50%   97.5% 
#>  45.925  84.500 161.025 
```
