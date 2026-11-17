# Predict Means from a Bayesian Scale Model

Predict expected responses, or linear predictors, at fitted or new
observations. Each prediction carries the joint posterior uncertainty in
the scales, the coefficients, and the random effects. A site or year
that was not in the fitted data needs an explicit choice, either drawing
its effect from the estimated spread among groups or setting it to zero.

## Usage

``` r
# S3 method for class 'multiScaleR_bayes'
predict(
  object,
  newdata = NULL,
  type = c("response", "link"),
  re.form = NULL,
  summary = TRUE,
  level = 0.95,
  draws = NULL,
  kernel_inputs = NULL,
  row_id = NULL,
  scale_re.form = NULL,
  new_scale_levels = c("error", "sample", "zero"),
  new_response_levels = c("error", "sample", "zero"),
  seed = NULL,
  include_fixed = FALSE,
  ...
)
```

## Arguments

- object:

  A fitted
  [`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md)
  object with posterior draws, from sampling or either approximation.

- newdata:

  `NULL` (default) uses the fitted observations. Otherwise a nonempty
  data frame of new observations that holds the fixed covariates, the
  offset variables, and the grouping variables of every component you
  include. Response values are neither needed nor used.

- type:

  `"response"` (default) returns expected responses, which are expected
  counts for Poisson and NB2 models, success probabilities for binomial
  models (also for grouped binomial data), and means for Gaussian
  models. `"link"` returns the linear predictor.

- re.form:

  `NULL` (default) includes the response random effects. `NA` sets them
  all to zero while keeping group-specific scales. This is a reference
  prediction, not an average over the random effects.

- summary:

  Logical. `TRUE` (default) returns posterior summaries. `FALSE` returns
  the matrix of draws by observations.

- level:

  Probability of the central interval; default 0.95.

- draws:

  Optional positive integer number of evenly spaced posterior draws to
  use. The default `NULL` uses every draw. Fewer draws save time and
  memory at a cost in simulation precision.

- kernel_inputs:

  Prepared landscape inputs for `newdata`, from
  [`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
  or
  [`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md).
  Required with new data and `NULL` otherwise. Use the fitted kernel,
  predictor definitions, and physical map units. See "Match new rows to
  their landscapes" in `predict.multiScaleR_bayes()`.

- row_id:

  Name of the column in `newdata` that matches each row to an entry of
  `kernel_inputs$binned$point_ids`, or to a row name of
  `kernel_inputs$kernel_dat` for surface-only inputs. IDs may repeat.
  The default `NULL` uses the fitted row-ID column when the fit kept it,
  and the row names of `newdata` otherwise. Supply it for fits saved by
  older versions of the package. Used only with new data.

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

- seed:

  Optional positive integer seed for sampled new effects. Your own
  random-number state is restored afterward, even when prediction fails.

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

With `summary = TRUE`, a data frame with one row per observation, in the
supplied order.

- `row`: the position of the observation.

- `point_id`: its prepared landscape ID.

- `mean`, `median`, `sd`, `lower`, `upper`: summaries of the predicted
  quantity on the scale set by `type`. `lower` and `upper` bound the
  central `level` credible interval.

- `new_scale_level`, `new_response_level` (new data only): whether the
  row belongs to an unseen group of an included scale or response
  component.

With `summary = FALSE`, a numeric matrix with one row per posterior draw
and one column per observation.

Results for new data carry a `new_groups` attribute, a data frame with
`component`, `covariate`, `group`, `level`, and `action` that lists each
unseen level and whether it was sampled or set to zero. Summaries also
carry a `prediction_contract` attribute that describes the target, the
included effects, the uncertainty, and the stored parameter map (see
[prediction_contract](https://wpeterman.github.io/multiScaleR/reference/prediction_contract.md)).
Its diagnostic status repeats the fit's own status. Prediction does not
make a failed fit eligible.

## What the intervals include

Intervals describe uncertainty in the expected response or the linear
predictor. They include uncertainty in the scales, the coefficients, and
the fitted group effects, and, when new effects are sampled, the
variation among new groups. They leave out the scatter of a single new
observation around its mean. Use
[`simulate.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/simulate.multiScaleR_bayes.md)
for that. Each interval is marginal for one observation, not a
simultaneous band. `mean` and `sd` summarize a finite set of draws and
do not show that the corresponding moment exists, so prefer the median
and the interval limits for long-tailed predictions.

## Match new rows to their landscapes

`kernel_inputs` supplies the prepared bins for the prediction locations
and maps. Each prepared location has one unique entry in
`binned$point_ids`, and `row_id` matches every row of `newdata` to one
of them. Rows can come in any order and can repeat, as for repeated
visits to one location. Build the inputs with
[`kernel_prep()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep.md)
or
[`kernel_prep_by_group()`](https://wpeterman.github.io/multiScaleR/reference/kernel_prep_by_group.md),
as in
[`vignette("v06-year-specific-rasters", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v06-year-specific-rasters.md).
Keep the fitted kernel, predictor definitions, and projected map units.
The numerical distance normalizer `unit_conv` may differ and is
converted internally, but that does not convert coordinate systems or
physical units. Kernel means need bins. Weighted SQ needs matching
source and metric specifications and the raw cell buffers.

Prediction reuses the fitted reference center and spread, the contrasts,
and the formula bases. It never standardizes the new rows on their own,
never learns a new scale reference, and never refits the model or uses
responses in `newdata`. A new kernel uses only the neighborhood you
prepared, so a predicted scale larger than `max_D` cannot draw on cells
beyond it.

## Choose how to treat new groups

[`vignette("v14-bayesian-prediction", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v14-bayesian-prediction.md)
compares the choices in a table. A few rules matter for every call. An
unseen level of a fixed year factor, such as `factor(year)`, is always
an error, under every policy. To predict response years you have not
surveyed, use a random intercept for year or a scientifically justified
numeric trend. Exchangeable annual effects impose no trend or
autocorrelation.
[`predict_scales()`](https://wpeterman.github.io/multiScaleR/reference/predict_scales.md)
can still predict a new year's scale when the response formula has fixed
year effects. Setting new effects to zero is not integration over group
variation, and with a nonlinear link or kernel the two predictions
generally differ. Grouping columns used only by a component you leave
out with `re.form = NA` or `scale_re.form = NA` may be absent from
`newdata`.

## Predict jointly across rows

Put every row whose joint prediction matters into one call. Within a
call, all new rows in the same grouping term and level share one new
effect for each posterior draw. Distinct terms and groups follow the
model's conditional independence, and the joint posterior dependence
among parameters is kept. Separate calls draw separate new effects, so
do not combine independently generated batches to infer joint spatial or
annual uncertainty or a total. For a total of expected values, set
`summary = FALSE` and sum each row of the matrix with
[`rowSums()`](https://rspatial.github.io/terra/reference/rowSums.html).

The full matrix of draws by observations is built even when
`summary = TRUE`, at 8 bytes per value. When you need only an interval
for each location of a large map, predict smaller batches and keep their
summaries, or lower `draws` at a cost in simulation precision.

## Use formula transformations safely

[`scale()`](https://rspatial.github.io/terra/reference/scale.html),
[`poly()`](https://rdrr.io/r/stats/poly.html),
[`splines::ns()`](https://rdrr.io/r/splines/ns.html), and
[`splines::bs()`](https://rdrr.io/r/splines/bs.html) keep their training
definitions inside ordinary expressions such as `I(scale(temp))` and in
transformed offsets, so single rows, batches, and effect curves use the
same constants. An expression such as `I(x - mean(x))` instead takes its
mean from the new batch, by ordinary R formula rules. Compute such a
covariate with training constants before fitting when you need the
original centering. Transforms must come from base R or stats, or carry
an explicit package namespace. Custom functions defined only in the
fitting workspace are not saved with the model, so compute them as
covariate columns instead. New levels of a fixed factor are always
rejected. For fits saved without design metadata, the rebuilt formula
must reproduce the stored training matrix from the original data, or
prediction stops. Give a namespace prefix to any custom function that
reuses a learner name such as `scale`.

## Know what prediction cannot fix

Prediction passes the fitted posterior forward. It does not repair
failed sampling, an inaccurate approximation, or a misspecified
observation model, and it does not check the fit's diagnostics, so
inspect `summary(object)$diagnostics` first. Every method uses its
stored joint draws. For an importance-corrected fit these are resampled
candidates, so prediction has the finite precision of those resamples
rather than a fresh weighted integration.

Prediction stops on unmatched IDs, new fixed-factor levels, unavailable
formula transforms, incompatible kernel inputs, and nonfinite results.
It never clips scales to the input limits or replaces a problematic
draw.

To test predictions, hold out data at the level you intend to predict. A
held-out visit at a known site and year is easier to predict than a
whole new site or year, and repeated sites or years are not independent
validation replicates.
[`cross_validate_bayes()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_bayes.md)
refits the model without the held-out rows and scores them.

## See also

[`multiScale_bayes()`](https://wpeterman.github.io/multiScaleR/reference/multiScale_bayes.md),
[multiScaleR_bayes](https://wpeterman.github.io/multiScaleR/reference/multiScaleR_bayes-class.md),
[`posterior_scales()`](https://wpeterman.github.io/multiScaleR/reference/posterior_scales.md),
[prediction_contract](https://wpeterman.github.io/multiScaleR/reference/prediction_contract.md),
[`predict_scales()`](https://wpeterman.github.io/multiScaleR/reference/predict_scales.md),
[`simulate.multiScaleR_bayes()`](https://wpeterman.github.io/multiScaleR/reference/simulate.multiScaleR_bayes.md),
[`cross_validate_bayes()`](https://wpeterman.github.io/multiScaleR/reference/cross_validate_bayes.md),
[multiScaleR-prediction](https://wpeterman.github.io/multiScaleR/reference/multiScaleR-prediction.md)
for prediction across the package;
[`vignette("v14-bayesian-prediction", package = "multiScaleR")`](https://wpeterman.github.io/multiScaleR/articles/v14-bayesian-prediction.md).

## Examples

``` r
example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
  package = "multiScaleR"))
fit <- example$fits$pooled
fitted_means <- predict(fit, draws = 100)
head(fitted_means[c("point_id", "median", "lower", "upper")])
#>   point_id    median    lower     upper
#> 1  s1_2021  3.248331 1.679910  4.551568
#> 2  s2_2021  4.039417 2.317309  7.140053
#> 3  s3_2021 15.054402 9.627522 20.920631
#> 4  s4_2021  4.283389 2.053277  6.393943
#> 5  s5_2021  5.032599 3.267581  7.769735
#> 6  s6_2021  4.818284 2.761564  6.980679
# Expected counts at the fitted rows, including fitted site effects.
# The draw limit keeps this demonstration short. Use all draws for
# final summaries.

# New surveys at known site-years reuse their annual maps and group effects
new_surveys <- example$data[1:3, c("site", "year", "site_year")]
predict(fit, new_surveys, kernel_inputs = example$inputs,
  row_id = "site_year", draws = 100)
#>   row point_id      mean    median        sd    lower     upper new_scale_level
#> 1   1  s1_2021  3.207937  3.248331 0.8470514 1.679910  4.551568           FALSE
#> 2   2  s2_2021  4.176684  4.039417 1.1973646 2.317309  7.140053           FALSE
#> 3   3  s3_2021 15.239878 15.054402 3.2055779 9.627522 20.920631           FALSE
#>   new_response_level
#> 1              FALSE
#> 2              FALSE
#> 3              FALSE
# simulate() would add the scatter of individual counts.

# A fit with random year intercepts can predict new sites and years
forecast <- readRDS(system.file("extdata", "bayesian-prediction-example.rds",
  package = "multiScaleR"))
future <- predict(forecast$fit, forecast$newdata,
  kernel_inputs = forecast$kernel_inputs, row_id = "site_year",
  new_scale_levels = "sample", new_response_levels = "sample",
  draws = 200, seed = 802)
future[c("point_id", "median", "lower", "upper", "new_response_level")]
#>          point_id    median     lower     upper new_response_level
#> 1         s1_2026  2.947675  1.347987  6.703234               TRUE
#> 2         s2_2026 12.196475  6.777054 22.035955               TRUE
#> 3 new_site_A_2026 19.774712  8.500104 46.423083               TRUE
#> 4 new_site_B_2026  4.269301  1.776990  8.776078               TRUE
#> 5         s1_2027  3.298215  1.528912  6.098285               TRUE
#> 6         s2_2027 13.633387  7.175980 29.062856               TRUE
#> 7 new_site_A_2027 22.960803 10.462792 47.634440               TRUE
#> 8 new_site_B_2027  4.550484  1.923561 11.587764               TRUE
attr(future, "new_groups")   # the unseen levels and what was done to them
#>   component covariate group      level action
#> 1     scale   habitat  year       2026 sample
#> 2     scale   habitat  year       2027 sample
#> 3  response            site new_site_A sample
#> 4  response            site new_site_B sample
#> 5  response            year       2026 sample
#> 6  response            year       2027 sample
```
