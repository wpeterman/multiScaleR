#' Simulate Posterior Predictive Responses
#'
#' Generate replicate responses for fitted or new observations. Each
#' replicate includes joint parameter uncertainty, any requested new-group
#' effects, and the scatter of individual observations around their means, so
#' replicates show the values you might record.
#'
#' @param object A fitted [multiScale_bayes()] object with posterior draws,
#'   from sampling or either approximation.
#' @param nsim Positive integer number of replicate data sets; default 1.
#' @param seed Optional positive integer simulation seed. Your own
#'   random-number state is restored after simulation, even when it fails.
#' @inheritParams predict.multiScaleR_bayes
#' @param trials Positive integer binomial trial totals, one value or one per
#'   observation row. The default `NULL` keeps the fitted totals when you
#'   simulate the fitted rows. For new observations it uses one trial for a
#'   binary fit, and it is an error for a grouped binomial fit, which needs
#'   explicit totals. Use only with binomial models. Response columns in
#'   `newdata` never supply the totals.
#' @param ... Not used. Supplying any further argument is an error.
#'
#' @return A data frame with one row per fitted or new observation and one
#'   column per replicate data set, named `sim_1`, `sim_2`, and so on. Values
#'   are counts for Poisson and NB2 models, success counts out of the fitted
#'   or supplied trial totals for binomial models (not proportions), and
#'   response units for Gaussian models. Each column uses
#'   one complete posterior vector, drawn at random with replacement, together
#'   with the requested random effects.
#'
#'   For new data, a `point_ids` attribute gives the prepared landscape ID of
#'   each row, and a `new_groups` attribute lists the unseen levels that were
#'   sampled or set to zero, as in `predict()`. A `prediction_contract`
#'   attribute (see [prediction_contract]) marks the result as posterior
#'   predictive replicates that include observation noise and records the
#'   policies for included effects. It does not establish diagnostic
#'   eligibility or calibration.
#'
#' @details
#' # Check the model against the data
#'
#' Compare summaries of the replicate data sets with the same summaries of
#' the observed response. An observed value far outside the replicated range
#' means the model cannot reproduce that feature. Good sampling diagnostics do
#' not show that the observation model fits the data, and agreement on one
#' statistic does not establish adequacy.
#' `vignette("v13-bayesian-scales", package = "multiScaleR")` works through
#' such a check. To use bayesplot, convert the replicates to a matrix with one
#' simulated data set per row, `t(as.matrix(replicated))`. For a count model,
#' `ppc_stat()` with `stat = function(y) mean(y == 0)` checks the fraction of
#' zeros.
#'
#' # Simulate new groups jointly
#'
#' Within one call, new-group effects are drawn afresh for each replicate
#' data set and shared by all of its rows in that group. Different new groups
#' are independent given the population parameters, whose joint posterior
#' dependence is kept. Observation noise is then drawn separately for each
#' row. Separate calls do not share a realization across components or
#' batches, so put every row that enters a joint question into one call.
#' Reusing a seed across [predict_scales()], `predict()`, and `simulate()`
#' does not promise matching latent effects, because the three functions use
#' their random numbers differently.
#'
#' To get an interval for a total, sum each replicate column first with
#' `colSums()` and then take quantiles of the sums. Adding the limits of
#' individual intervals overstates the range of the total.
#'
#' The three prediction functions return different orientations, which
#' matters for joint contrasts and totals. `simulate()` returns observations
#' by replicates. `predict(summary = FALSE)` returns posterior draws by
#' observations, and `predict_scales(summary = FALSE)` returns a named list of
#' draws-by-observations matrices, one per landscape predictor.
#'
#' # Supply trials for binomial data
#'
#' For a binary fit, new observations default to one trial. For a fit to
#' `cbind(successes, failures)`, supply the actual future totals, for example
#' `trials = rep(10L, nrow(newdata))`. `simulate()` then returns success
#' counts out of those totals, while `predict()` returns success
#' probabilities whatever the number of trials.
#'
#' @seealso [multiScale_bayes()], [simulate_bayes_prior()],
#'   [predict.multiScaleR_bayes()], [predict_scales()],
#'   [prediction_contract], [multiScaleR_bayes-class];
#'   `vignette("v14-bayesian-prediction", package = "multiScaleR")`.
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' replicated <- simulate(example$fits$pooled, nsim = 20, seed = 94)
#' dim(replicated)  # one row per original observation, one column per replicate
#'
#' # Posterior predictive check of the fraction of zero counts
#' range(vapply(replicated, function(y) mean(y == 0), numeric(1)))
#' mean(example$data$count == 0)  # observed fraction of zero counts
#' # These replicate counts include observation noise. Fitted means do not.
#'
#' # Counts at eight future surveys, and an interval for their total
#' forecast <- readRDS(system.file("extdata", "bayesian-prediction-example.rds",
#'   package = "multiScaleR"))
#' future <- simulate(forecast$fit, nsim = 200, seed = 803,
#'   newdata = forecast$newdata, kernel_inputs = forecast$kernel_inputs,
#'   row_id = "site_year", new_scale_levels = "sample",
#'   new_response_levels = "sample")
#' quantile(colSums(future), c(0.025, 0.5, 0.975))  # total of the eight surveys
#' @exportS3Method stats::simulate
simulate.multiScaleR_bayes <- function(object, nsim = 1L, seed = NULL,
    newdata = NULL, kernel_inputs = NULL, row_id = NULL, re.form = NULL,
    scale_re.form = NULL, new_scale_levels = c("error", "sample", "zero"),
    new_response_levels = c("error", "sample", "zero"), trials = NULL, include_fixed = FALSE, ...) {
  if (length(list(...))) stop("Unused simulation arguments.", call. = FALSE)
  validate_scalar_numeric(nsim, "nsim", integerish = TRUE, lower = 1)
  ids <- .msr_bayes_draw_ids(object, NULL)
  include_response <- .msr_bayes_prediction_re_form(re.form, "re.form")
  context <- .msr_bayes_prediction_seed(seed,
    .msr_bayes_prediction_prepare(object, newdata, kernel_inputs, row_id,
      match.arg(new_scale_levels), match.arg(new_response_levels),
      .msr_bayes_prediction_re_form(scale_re.form, "scale_re.form"), include_response,
      include_fixed))
  spec <- context$spec
  if (!is.null(trials)) {
    if (spec$family != "binomial") stop("trials applies only to binomial models.", call. = FALSE)
    if (!is.numeric(trials) || !is.null(dim(trials)) ||
        !length(trials) %in% c(1L, length(spec$y)) ||
        any(!is.finite(trials) | trials < 1 | trials != floor(trials)))
      stop("trials must be positive integer totals, one or one per prediction row.", call. = FALSE)
    spec$trials <- rep(trials, length.out = length(spec$y))
  } else if (!is.null(newdata) && spec$family == "binomial" &&
      any(object$specification$trials != 1))
    stop("Supply trials explicitly for new grouped-binomial observations.", call. = FALSE)
  if (!is.null(seed)) {
    validate_scalar_numeric(seed, "seed", integerish = TRUE, lower = 1,
                             upper = .Machine$integer.max)
    had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
    on.exit({
      if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  ans <- matrix(NA_real_, length(spec$y), nsim)
  selected <- sample(ids, nsim, replace = TRUE)
  for (j in seq_len(nsim)) {
    p <- .msr_bayes_prediction_parameters(.msr_bayes_draw(object, selected[j]), context)
    eta <- .msr_bayes_components(spec, p, include_response = include_response)$eta
    ans[, j] <- switch(spec$family,
      gaussian = stats::rnorm(length(eta), eta, exp(p$log_dispersion)),
      poisson = stats::rpois(length(eta), exp(eta)),
      binomial = stats::rbinom(length(eta), spec$trials, stats::plogis(eta)),
      nbinom2 = stats::rnbinom(length(eta), mu = exp(eta), size = exp(p$log_dispersion)))
    if (!is.null(newdata) && any(!is.finite(ans[, j])))
      stop("Nonfinite new-observation replicate; inspect the design, priors and parameter draws.", call. = FALSE)
  }
  colnames(ans) <- paste0("sim_", seq_len(nsim))
  ans <- as.data.frame(ans)
  if (!is.null(newdata)) {
    attr(ans, "point_ids") <- spec$observation_ids
    attr(ans, "new_groups") <- context$new_groups
  }
  attr(ans, "prediction_contract") <- .msr_bayes_prediction_contract(
    object, context, "replicate_response", "response",
    .msr_bayes_prediction_re_form(scale_re.form, "scale_re.form"), include_response)
  ans
}
