#' Simulate Prior Predictive Responses and Landscape Scales
#'
#' Check the implications of the complete model priors before fitting. Draw
#' parameters from their specified priors and generate responses for the
#' prepared observations, including observation sampling variation.
#'
#' @param object A \code{\link{multiScale_bayes}} object. An object from
#'   \code{method = "prepare"} is sufficient. A fitted object's posterior
#'   draws are ignored; its specification supplies the same prior model.
#' @param nsim Positive integer number of prior predictive datasets; default
#'   100. Larger values explore the priors more thoroughly and use more memory.
#' @param seed Optional positive integer simulation seed. When supplied, the
#'   caller's random-number state is restored, including when simulation fails.
#' @return A list with:
#'   \describe{
#'     \item{responses}{A data frame with the original observation rows and
#'       one column per simulated dataset, named \code{sim_1}, \code{sim_2},
#'       and so on. Gaussian responses use response units. Poisson and NB2
#'       responses are counts. Binomial responses are success counts out of
#'       the prepared trial totals.}
#'     \item{means}{A data frame of matching expected responses. Binomial
#'       values are success probabilities, not expected success counts;
#'       multiply by the prepared trial totals for expected counts. Other
#'       families use the response scale. Each mean conditions on its replicate's
#'       prior parameter vector and excludes observation variation; it is not
#'       an expectation integrated over the entire prior.}
#'     \item{log_scales}{A named list, one matrix per landscape predictor.
#'       Rows follow its distinct scale-formula combinations, labeled as in
#'       \code{\link{posterior_scales}}; columns match the simulated datasets.
#'       Entries are natural log distances in projected map units. Exponentiate
#'       medians and interval limits for distance summaries.}
#'     \item{observation_ids}{The prepared observation identities, in the same
#'       order as the rows of \code{responses} and \code{means}.}
#'     \item{priors}{The complete prior settings used.}
#'   }
#' @details Each dataset uses one complete, independent prior parameter vector.
#'   All observations in a group share that vector's group effect, retaining
#'   the dependence specified by the scale and response formulas. Scale effects
#'   and response effects are drawn separately, as specified by the model.
#'
#'   Simulation conditions on the prepared maps, fixed covariates, offsets,
#'   group identities, binomial trial totals and frozen predictor
#'   standardization. It does not use observed response values, posterior
#'   draws or likelihood weights. It does not generate new sites or years.
#'   No optional Bayesian packages are required.
#'
#'   Compare replicated magnitudes, zero fractions and scale quantiles with
#'   scientifically plausible ranges before fitting. Extreme predictions can
#'   reveal a prior that is implausible for your design. They do not demonstrate
#'   poor posterior fit or prescribe a prior automatically. Broad distance
#'   priors can have infinite arithmetic moments; prefer log-scale summaries
#'   and distance quantiles. Broad log-link response group-SD priors can also
#'   permit infinite prior predictive count moments. A finite simulation does
#'   not establish that an arithmetic mean over prior replicates exists.
#'   Nonfinite calculations cause an informative error
#'   identifying the replicate instead of clipping values or replacing a draw.
#' @seealso \code{\link{multiScale_bayes}},
#'   \code{\link{simulate.multiScaleR_bayes}},
#'   \code{vignette("bayesian-scales", package = "multiScaleR")}.
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' prepared <- multiScale_bayes(count ~ habitat + factor(year) + (1 | site),
#'   example$inputs, example$data, scales = list(habitat = ~1 + (1 | year)),
#'   row_id = "site_year", priors = list(scale_median = 250), method = "prepare")
#' prior <- simulate_bayes_prior(prepared, nsim = 20, seed = 63)
#' range(vapply(prior$responses, function(y) mean(y == 0), numeric(1)))
#' # Prior responses can be much more variable than the observed counts.
#' exp(apply(prior$log_scales$habitat, 1, stats::quantile, c(.025, .5, .975)))
#' # These are prior scale intervals, not posterior estimates or evidence of fit.
#' @export
simulate_bayes_prior <- function(object, nsim = 100L, seed = NULL) {
  if (!inherits(object, "multiScaleR_bayes"))
    stop("object must be a multiScaleR_bayes object.", call. = FALSE)
  validate_scalar_numeric(nsim, "nsim", integerish = TRUE, lower = 1)
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
  spec <- object$specification
  priors <- spec$priors
  responses <- means <- matrix(NA_real_, length(spec$y), nsim)
  log_scales <- lapply(spec$scales, function(sm)
    matrix(NA_real_, length(sm$representatives), nsim,
      dimnames = list(sm$labels, paste0("sim_", seq_len(nsim)))))
  for (j in seq_len(nsim)) {
    # Draw from the actual model priors. Importance-proposal helpers also
    # contain proposal-specific distributions and must not be used here.
    p <- spec$parameters
    p$beta <- stats::rnorm(length(p$beta), sd = priors$beta_sd)
    for (nm in names(spec$scales)) {
      sm <- spec$scales[[nm]]
      location <- if (is.null(sm$prior_mean))
        c(log(priors$scale_median[[nm]] / spec$unit_conv),
          rep(0, length(sm$beta_index) - 1L)) else sm$prior_mean
      spread <- if (is.null(sm$prior_sd))
        c(priors$scale_intercept_sd, rep(priors$scale_coef_sd,
          length(sm$beta_index) - 1L)) else sm$prior_sd
      p$scale_beta[sm$beta_index] <- stats::rnorm(length(sm$beta_index), location, spread)
    }
    p$scale_z <- stats::rnorm(length(p$scale_z))
    p$response_z <- stats::rnorm(length(p$response_z))
    p$log_scale_sd <- log(abs(stats::rnorm(length(p$log_scale_sd), sd = priors$scale_group_sd)))
    p$log_response_sd <- log(abs(stats::rnorm(length(p$log_response_sd), sd = priors$response_group_sd)))
    p$log_dispersion <- stats::rnorm(length(p$log_dispersion),
      priors$dispersion_log_mean, priors$dispersion_log_sd)
    if (any(!is.finite(unlist(p, use.names = FALSE))))
      stop("Nonfinite prior parameter in replicate ", j,
        "; inspect the prior settings. No draw was replaced.", call. = FALSE)
    eta <- .msr_bayes_components(spec, p)$eta
    mu <- switch(spec$family, gaussian = eta, binomial = stats::plogis(eta), exp(eta))
    if (any(!is.finite(eta)) || any(!is.finite(mu)))
      stop("Nonfinite prior predictive calculation in replicate ", j,
        "; inspect the prior settings and prepared design. No draw was replaced.", call. = FALSE)
    means[, j] <- mu
    responses[, j] <- switch(spec$family,
      gaussian = stats::rnorm(length(eta), mu, exp(p$log_dispersion)),
      poisson = stats::rpois(length(eta), mu),
      binomial = stats::rbinom(length(eta), spec$trials, mu),
      nbinom2 = stats::rnbinom(length(eta), mu = mu, size = exp(p$log_dispersion)))
    if (any(!is.finite(responses[, j])))
      stop("Nonfinite prior predictive response in replicate ", j,
        "; inspect the prior settings. No draw was replaced.", call. = FALSE)
    for (nm in names(spec$scales)) {
      sm <- spec$scales[[nm]];rows <- sm$representatives
      value <- drop(sm$X[rows, , drop = FALSE] %*% p$scale_beta[sm$beta_index])
      for (term in sm$random)
        value <- value + exp(p$log_scale_sd[term$sd_index]) *
          p$scale_z[term$z_index][term$index[rows]]
      if (any(!is.finite(value)))
        stop("Nonfinite prior log scale in replicate ", j,
          "; inspect the prior settings. No draw was replaced.", call. = FALSE)
      log_scales[[nm]][, j] <- value + log(spec$unit_conv)
    }
  }
  colnames(responses) <- colnames(means) <- paste0("sim_", seq_len(nsim))
  list(responses = as.data.frame(responses), means = as.data.frame(means),
    log_scales = log_scales, observation_ids = spec$observation_ids, priors = priors)
}
