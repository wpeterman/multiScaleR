#' Simulate Posterior Predictive Responses
#'
#' Generate replicate responses for fitted or new observations, including
#' joint parameter uncertainty, explicit new-group effects and observation
#' sampling variation.
#'
#' @param object A \code{\link{multiScale_bayes}} fit with posterior draws
#'   from sampling or either approximation.
#' @param nsim Positive number of replicate datasets; default 1.
#' @param seed Optional positive integer simulation seed. When supplied, the
#'   caller's random-number state is restored after simulation, including on failure.
#' @inheritParams predict.multiScaleR_bayes
#' @param trials Positive integer binomial trial totals, one scalar or one per
#'   observation row. \code{NULL} keeps fitted totals for fitted-row simulation
#'   and uses one trial for new observations from a binary fit. New observations
#'   from a grouped-binomial fit require explicit totals. This argument applies
#'   only to binomial models. New response values never supply the trial totals.
#' @param ... Reserved for extensions; unsupported arguments cause an error.
#' @return A data frame with fitted or new observation rows and one column per
#'   replicate dataset. Binomial values are success counts out of the fitted
#'   or supplied trial totals, not proportions. Each column uses one sampled
#'   complete posterior vector and the requested included random effects.
#'   For new data, \code{point_ids} and \code{new_groups} attributes identify
#'   original row order and sampled/zeroed levels, as in \code{predict()}.
#'   The \code{prediction_contract} attribute identifies posterior predictive
#'   replicates, which include observation noise, and records the included
#'   effect policies. It does not establish diagnostic eligibility or calibration.
#' @details Compare summaries of these replicate datasets with the observed
#'   response to assess model adequacy. Good sampling diagnostics do not
#'   establish that the observation model fits the data.
#'   For new groups, effects are sampled afresh for each replicate dataset,
#'   shared by its rows in that group and independent across replicate datasets
#'   conditional on the sampled posterior vector. Observation noise is then
#'   generated separately for each row. Use a single batch to retain dependence.
#'   Summaries of individual predictive intervals cannot be summed to obtain
#'   an interval for a total; sum each complete replicate instead.
#' @seealso \code{\link{multiScale_bayes}}, \code{\link{simulate_bayes_prior}},
#'   \code{\link{prediction_contract}},
#'   \code{\link{predict.multiScaleR_bayes}};
#'   \code{vignette("bayesian-prediction", package = "multiScaleR")}.
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' replicated <- simulate(example$fits$pooled, nsim = 20, seed = 94)
#' dim(replicated)  # one row per original observation, one column per replicate
#' range(vapply(replicated, function(y) mean(y == 0), numeric(1)))
#' mean(example$data$count == 0)  # observed fraction of zero counts
#' # These replicate counts include observation variation; fitted means do not.
#' # Agreement on one statistic alone does not establish model adequacy.
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
