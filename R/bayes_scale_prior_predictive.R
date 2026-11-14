#' Simulate Prior Predictive Responses and Landscape Scales
#'
#' Check what the complete set of priors implies before you fit. Each
#' simulated data set draws every parameter from its prior and then generates
#' responses for the prepared observations, including the scatter of
#' individual observations around their means.
#'
#' @param object A [multiScale_bayes()] object. An object from
#'   `method = "prepare"` is enough. For a fitted object, the posterior draws
#'   are ignored and its specification supplies the same prior model.
#' @param nsim Positive integer number of prior predictive data sets; default
#'   100. More data sets explore the priors more thoroughly and use more
#'   memory.
#' @param seed Optional positive integer simulation seed. Your own
#'   random-number state is restored afterward, even when simulation fails.
#'
#' @return A list with these elements.
#'   * `responses`: a data frame with one row per prepared observation and
#'     one column per simulated data set, named `sim_1`, `sim_2`, and so on.
#'     Values are counts for Poisson and NB2 models, success counts out of the
#'     prepared trial totals for binomial models, and response units for
#'     Gaussian models.
#'   * `means`: a data frame of the matching expected responses. Binomial
#'     values are success probabilities, so multiply by the trial totals for
#'     expected counts. Other families use the response scale. Each mean is
#'     conditional on its own prior draw and excludes observation noise. It is
#'     not an expectation over the whole prior.
#'   * `log_scales`: a named list with one matrix per landscape predictor. Rows
#'     are its distinct combinations of scale covariates and groups, labeled as
#'     in [posterior_scales()], and columns match the simulated data sets.
#'     Entries are natural log distances in projected map units, so
#'     exponentiate medians and interval limits to read them as distances.
#'   * `observation_ids`: the prepared landscape ID of each observation, in
#'     the row order of `responses` and `means`.
#'   * `priors`: the complete prior settings used, with defaults filled in.
#'
#' @details
#' # Read the simulations
#'
#' Compare the simulated magnitudes, fractions of zeros, and scale quantiles
#' with ranges that are scientifically plausible for your design. Extreme
#' simulations show that a prior is loose for your design and may call for a
#' tighter setting, such as a smaller `beta_sd`. They are not evidence of a
#' poor posterior fit, and they do not choose a prior for you.
#' `vignette("v13-bayesian-scales", package = "multiScaleR")` works through
#' this check.
#'
#' Broad distance priors can have infinite arithmetic moments, so summarize
#' scales with log-scale summaries and distance quantiles. Broad response
#' group-SD priors with a log link can also give infinite prior predictive
#' count moments. A finite simulation does not show that a mean over prior
#' replicates exists.
#'
#' # How the simulation works
#'
#' Each data set uses one complete, independent prior parameter vector. All
#' observations in a group share that vector's group effect, which keeps the
#' dependence that the scale and response formulas specify. Scale effects and
#' response effects are drawn separately, as the model specifies.
#'
#' Simulation conditions on the prepared maps, fixed covariates, offsets,
#' group identities, binomial trial totals, and the frozen predictor
#' standardization. It uses no observed responses and no posterior draws, and
#' it creates no new sites or years. No optional Bayesian package is needed. A
#' nonfinite value stops with an error that names the replicate. Values are
#' never clipped and draws are never replaced.
#'
#' @seealso [multiScale_bayes()], [simulate.multiScaleR_bayes()],
#'   [posterior_scales()];
#'   `vignette("v13-bayesian-scales", package = "multiScaleR")`.
#' @examples
#' example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'   package = "multiScaleR"))
#' prepared <- multiScale_bayes(count ~ habitat + factor(year) + (1 | site),
#'   example$inputs, example$data, scales = list(habitat = ~1 + (1 | year)),
#'   row_id = "site_year", priors = list(scale_median = 250), method = "prepare")
#' prior <- simulate_bayes_prior(prepared, nsim = 20, seed = 63)
#'
#' # Fraction of zero counts in each prior data set, against the observed 0.04
#' range(vapply(prior$responses, function(y) mean(y == 0), numeric(1)))
#' # The default priors are loose and allow far more zeros than were observed.
#'
#' # Prior range of sigma for each year, in meters
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
