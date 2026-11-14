# Proposal calculations use unconstrained parameter coordinates. Remove beta
# from the outer covariance before drawing its normalized conditional density.
.msr_laplace_audit_joint_proposal <- function(object, keep) {
  original <- object$laplace
  beta <- setdiff(seq_along(original$mode), keep)
  covariance <- original$covariance[keep, keep, drop = FALSE]
  hessian <- chol2inv(chol(covariance))
  geometry <- .msr_bayes_laplace_geometry(hessian)
  list(node = list(mode = original$mode[keep], hessian = hessian, root = geometry$root),
    start = original$mode[beta])
}

.msr_laplace_audit_draws <- function(object, n) {
  spec <- object$specification
  sizes <- lengths(object$parameters)
  position <- seq_len(sum(sizes))
  keep <- setdiff(position, seq_len(sizes[["beta"]]))
  is_joint <- !is.null(object$laplace)
  if (is_joint) {
    project <- .msr_laplace_audit_joint_proposal(object, keep)
    prior_object <- object
    # Reuse the latent-prior density and its local dispersion t component.
    prior_object$marginal <- list(hyper_index = integer(), weights = 1,
      nodes = list(object$laplace))
    draw_local <- function(n) .msr_bayes_importance_tdraw(project$node, integer(), n)
    log_local <- function(x) .msr_bayes_importance_density(x, project$node)
  } else {
    original <- object$marginal
    map <- list(par = stats::setNames(numeric(length(position)), rep(names(sizes), sizes)))
    project <- .msr_bayes_beta_project(original, map)
    grid <- project$grid
    mixture <- .msr_bayes_marginal_weights(grid$x, log(.99 * grid$density + .01))
    grid[names(mixture)] <- mixture
    grid$tail_fraction <- .01
    prior_object <- object
    draw_local <- function(n) .msr_bayes_importance_candidates(grid, project$hyper, n)
    log_local <- function(x) .msr_bayes_importance_logq(x, grid, project$hyper)
  }
  prior <- .msr_bayes_latent_prior_proposal(prior_object, keep)
  defensive <- stats::runif(n) < .2
  outer <- matrix(0, n, length(keep))
  if (any(defensive))
    outer[defensive, ] <- .msr_bayes_latent_prior_draws(prior, sum(defensive))
  if (any(!defensive))
    outer[!defensive, ] <- draw_local(sum(!defensive))
  logq <- .msr_bayes_logadd(log(.8) + log_local(outer),
    log(.2) + .msr_bayes_latent_prior_logq(outer, prior))
  beta <- .msr_bayes_beta_candidates(outer, keep, spec, project$start)
  list(samples = beta$samples, logq = logq + beta$logq,
    beta = c(fallbacks = beta$fallbacks, max_step = beta$max_step),
    defensive_draws = sum(defensive))
}

.msr_laplace_audit_log_scale <- function(draws, spec, covariate, row) {
  sm <- spec$scales[[covariate]]
  value <- drop(draws$scale_beta[, sm$beta_index, drop = FALSE] %*% sm$X[row, ])
  for (term in sm$random)
    value <- value + exp(draws$log_scale_sd[, term$sd_index]) *
      draws$scale_z[, term$z_index[term$index[row]]]
  value + log(spec$unit_conv)
}

.msr_laplace_audit_event <- function(approx, reference, weight, tolerance, overlap) {
  if (anyNA(approx) || anyNA(reference))
    stop("Nonfinite fitted scale values in approximation assessment.", call. = FALSE)
  p <- mean(approx)
  r <- sum(weight * reference)
  mcse <- sqrt(sum(weight^2 * (reference - r)^2))
  half <- seq_along(weight) <= length(weight) / 2
  split <- abs(sum(weight[half] * reference[half]) / sum(weight[half]) -
    sum(weight[!half] * reference[!half]) / sum(weight[!half]))
  # Wilson bounds avoid claiming zero simulation error for an all-zero or
  # all-one approximate sample. The weighted reference uses a delta MCSE.
  z <- 3; n <- length(approx)
  center <- (p + z^2 / (2 * n)) / (1 + z^2 / n)
  radius <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n)
  left <- center - radius - r - z * mcse
  right <- center + radius - r + z * mcse
  event_ess <- if (r > 0 && r < 1 && mcse > 0) r * (1 - r) / mcse^2 else NA_real_
  status <- if (!overlap || !is.finite(event_ess) || event_ess < 400 ||
      !is.finite(split) || split > min(.05, tolerance)) "unresolved" else
    if (left > tolerance || right < -tolerance) "discrepancy" else
      if (left >= -tolerance && right <= tolerance) "within_tolerance" else "unresolved"
  data.frame(approximation = p, reference = r, reference_mcse = mcse,
    difference = p - r, difference_lower = left, difference_upper = right,
    indicator_ess = event_ess, split_reference_difference = split, status = status)
}

#' Assess Fitted Scale Probabilities from an Uncorrected Laplace Fit
#'
#' Check the probability that each fitted scale lies at or below a distance
#' you name, as reported by a joint Laplace or uncorrected marginal Laplace
#' fit, against a fresh importance-weighted estimate of the same posterior.
#' Read this check alongside the fit's own optimizer and integration
#' diagnostics.
#'
#' This function is experimental. Its arguments and results may change in a
#' later release without a deprecation period.
#'
#' @param object A fitted [multiScale_bayes()] object from
#'   `method = "laplace"`, or from `method = "marginal"` without importance
#'   correction, with its approximation draws and geometry intact. The check
#'   needs the RTMB and posterior (version 1.6.0 or later) packages.
#' @param thresholds Named list of scale thresholds, as in
#'   [diagnose_bayes_probabilities()]. Each name selects a fitted landscape
#'   predictor, and each vector holds distinct finite thresholds, such as
#'   `list(habitat = c(150, 350))`. Every threshold is checked for every
#'   distinct fitted combination of scale covariates and groups, such as each
#'   year.
#' @param scale Units of `thresholds`. `"distance"` (default) uses map units
#'   and needs positive values. `"log"` uses natural log map units. The event
#'   is a fitted kernel scale less than or equal to the threshold.
#' @param tolerance Absolute difference in probability that matters for your
#'   inference, strictly between 0 and 1; default 0.05. Choose it, and the
#'   thresholds, before looking at the results.
#' @param candidates Number of fresh importance candidates, at least 2000;
#'   default 16384. The validation study used 65536. More candidates can
#'   reduce simulation error when the proposal overlaps the posterior well.
#'   Every candidate is kept for the weight and event diagnostics, so memory
#'   use grows with this number and with the number of parameters.
#' @param seed Optional positive integer seed for the assessment. Your own
#'   random-number state is restored afterward.
#'
#' @return A list with these elements.
#'   * `available`: `TRUE` when importance ratios and reference summaries
#'     could be calculated. A reference that fails its overlap checks is still
#'     available, but its statuses stay `"unresolved"`.
#'   * `reason`: why the reference is unavailable, why it failed its overlap
#'     checks, or that the original fit failed its own diagnostics. An empty
#'     string otherwise.
#'   * `raw_diagnostics`: the fit's original `diagnostics`, unchanged.
#'   * `overlap`: the importance checks of the reference (`passed`,
#'     `pareto_k`, `ess`, `min_tail_ess`, `max_weight`, `split_cdf_error`), or
#'     `NULL` when the reference is unavailable.
#'   * `target_evaluation`: counts of posterior density evaluations and of
#'     rechecks of invalid values.
#'   * `probabilities`: a data frame with one row per fitted scale combination
#'     and threshold. `covariate`, `group`, `row`, `threshold`, and
#'     `threshold_scale` identify the event. `approximation` is the
#'     probability from the fit's draws and `reference` the importance-weighted
#'     probability, with its Monte Carlo standard error `reference_mcse`.
#'     `difference` is the approximation minus the reference, and
#'     `difference_lower` and `difference_upper` bound it numerically.
#'     `indicator_ess` is the effective sample size of the reference event
#'     indicator, `split_reference_difference` the absolute difference between
#'     the two halves of the reference, and `status` the verdict described in
#'     "Read the status".
#'   * `candidates`, `seed`, `tolerance`: the settings used.
#'   * `coefficient_proposal`: the number of fallback coefficient proposals
#'     (`fallbacks`) and the largest remaining Newton step (`max_step`).
#'   * `defensive_draws`: the number of candidates drawn from the defensive
#'     prior component.
#'   * `metadata_establishes_eligibility`: always `FALSE`.
#'
#' @details
#' # Read the status
#'
#' | `status` | Meaning | What to do |
#' |:---------|:--------|:-----------|
#' | `"discrepancy"` | The numerical bounds on the difference lie wholly beyond `tolerance` on one side | Do not report the approximate probability. Sample the model |
#' | `"within_tolerance"` | Both bounds lie inside `tolerance` | You can report the approximate probability, for this threshold only |
#' | `"unresolved"` | Every other case, including a reference too imprecise or overlapping the posterior too poorly to decide | Treat the probability as unchecked. Raise `candidates`, or sample |
#'
#' An event can be resolved only when the reference passes the package's
#' importance overlap checks, the original fit passed its own diagnostics, the
#' event has an indicator ESS of at least 400, and the two independent halves
#' of the reference differ by no more than the smaller of 0.05 and
#' `tolerance`. The bounds combine three standard errors of the weighted
#' reference with Wilson bounds for the original approximate draws. A constant
#' reference indicator also leaves the event unresolved. Expect
#' `"unresolved"` often, and read it as "not checked", not as "close enough".
#' `vignette("v15-bayesian-computation", package = "multiScaleR")` reports how
#' often each status occurred in simulations.
#'
#' # How the reference is built
#'
#' Candidates come 80% from heavy-tailed local proposals around the
#' approximation and 20% from a defensive component based on the latent-effect
#' and scale priors. Response coefficients are drawn from their exact Gaussian
#' conditional, or from a normalized heavy-tailed conditional proposal for
#' other families. The full posterior density and the actual mixture density
#' determine every weight. Neither the approximation nor its draws are
#' replaced.
#'
#' This exploratory check targets only the requested scale probabilities. A
#' heavy-tailed proposal can still miss posterior regions. The report does not
#' certify approximation accuracy, general calibration, held-out predictions,
#' or the fit's original diagnostics. Compare with a full posterior fit
#' whenever your conclusions depend on these probabilities.
#'
#' @seealso [multiScale_bayes()], [diagnose_bayes_probabilities()],
#'   [posterior_scales()], [multiScaleR_bayes-class];
#'   `vignette("v15-bayesian-computation", package = "multiScaleR")`.
#' @examples
#' \donttest{
#' if (requireNamespace("RTMB", quietly = TRUE) &&
#'     requireNamespace("posterior", quietly = TRUE)) {
#'   example <- readRDS(system.file("extdata", "bayesian-annual-example.rds",
#'     package = "multiScaleR"))
#'   # A joint Laplace fit of the annual model, which takes about a second
#'   fast_fit <- multiScale_bayes(count ~ habitat + factor(year) + (1 | site),
#'     example$inputs, example$data, scales = list(habitat = ~1 + (1 | year)),
#'     row_id = "site_year", family = "nbinom2",
#'     priors = list(scale_median = 250, scale_group_sd = 0.5),
#'     method = "laplace", seed = 527001)
#'
#'   # Probability that each year's sigma is 150 m or less. The smallest
#'   # candidate budget keeps this example short. The validation study used
#'   # 65536 candidates.
#'   check <- diagnose_bayes_laplace(fast_fit,
#'     thresholds = list(habitat = 150), candidates = 2000, seed = 127)
#'   print(check$reason)   # empty when the reference and the fit passed
#'   print(check$probabilities[c("group", "approximation", "reference",
#'                               "status")])
#'   # "unresolved" means the probability was not checked, not that it is
#'   # close enough.
#' }
#' }
#' @export
diagnose_bayes_laplace <- function(object, thresholds, scale = c("distance", "log"),
    tolerance = .05, candidates = 16384L, seed = NULL) {
  if (!inherits(object, "multiScaleR_bayes") ||
      is.null(object$laplace) == is.null(object$marginal) ||
      !is.null(object$importance) || is.null(object$draws))
    stop("Supply an uncorrected joint or marginal Laplace fit with retained draws.", call. = FALSE)
  scale <- match.arg(scale)
  spec <- object$specification
  if (.msr_has_gaussian_blocks(spec))
    stop("Private Gaussian projection blocks are outside this assessment.", call. = FALSE)
  .msr_bayes_validate_probability_thresholds(thresholds, spec, scale)
  validate_scalar_numeric(tolerance, "tolerance", lower = 0, upper = 1,
    inclusive_lower = FALSE, inclusive_upper = FALSE)
  validate_scalar_numeric(candidates, "candidates", integerish = TRUE, lower = 2000,
    upper = .Machine$integer.max)
  if (!requireNamespace("RTMB", quietly = TRUE) ||
      !.msr_bayes_package_status("posterior", "1.6.0")$passed)
    stop("Install RTMB and posterior (>= 1.6.0) to assess a Laplace fit.", call. = FALSE)
  out <- list(available = FALSE, reason = "", raw_diagnostics = object$diagnostics,
    overlap = NULL, target_evaluation = NULL, candidates = candidates,
    seed = seed, tolerance = tolerance, coefficient_proposal = NULL,
    defensive_draws = NA_integer_, metadata_establishes_eligibility = FALSE)
  audit <- tryCatch({
    population <- .msr_bayes_prediction_seed(seed, .msr_laplace_audit_draws(object, candidates))
    out$coefficient_proposal <- population$beta
    out$defensive_draws <- population$defensive_draws
    obj <- .msr_bayes_ad_object(spec)
    target <- .msr_bayes_importance_target(population$samples, obj, spec)
    out$target_evaluation <- target$diagnostics
    log_ratio <- target$log_density - population$logq
    if (!all(is.finite(population$logq)) || anyNA(log_ratio) ||
        any(log_ratio == Inf) || !any(is.finite(log_ratio)))
      stop("Invalid importance ratios. No candidate was dropped.", call. = FALSE)
    checks <- .msr_bayes_importance_checks(population$samples, log_ratio, spec)
    list(checks = checks, fresh = .msr_bayes_parameter_blocks(population$samples, spec$parameters))
  }, error = identity)
  valid <- !inherits(audit, "error")
  if (!valid) out$reason <- conditionMessage(audit) else {
    checks <- audit$checks
    out$available <- TRUE
    out$overlap <- checks[setdiff(names(checks), "weights")]
    if (!isTRUE(checks$passed))
      out$reason <- "The reference did not pass the posterior importance overlap checks."
    if (!isTRUE(object$diagnostics$passed))
      out$reason <- paste(out$reason, "The original approximation did not pass its fit diagnostics.")
  }
  probabilities <- list()
  for (nm in names(thresholds)) {
    sm <- spec$scales[[nm]]
    for (j in seq_along(sm$representatives)) {
      row <- sm$representatives[j]
      approx <- .msr_laplace_audit_log_scale(object$draws, spec, nm, row)
      if (any(!is.finite(approx)))
        stop("Nonfinite fitted scale values in the original approximation.", call. = FALSE)
      reference <- if (valid) .msr_laplace_audit_log_scale(audit$fresh, spec, nm, row)
      reference_valid <- valid && all(is.finite(reference))
      for (threshold in thresholds[[nm]]) {
        limit <- if (scale == "log") threshold else log(threshold)
        a <- approx <= limit
        result <- if (reference_valid) .msr_laplace_audit_event(a, reference <= limit,
          checks$weights, tolerance, isTRUE(checks$passed) && isTRUE(object$diagnostics$passed)) else
            data.frame(approximation = mean(a), reference = NA_real_, reference_mcse = NA_real_,
              difference = NA_real_, difference_lower = NA_real_, difference_upper = NA_real_,
              indicator_ess = NA_real_, split_reference_difference = NA_real_, status = "unresolved")
        id <- data.frame(covariate = nm, group = sm$labels[j], row = as.integer(row),
          threshold = unname(threshold), threshold_scale = scale)
        probabilities[[length(probabilities) + 1L]] <- cbind(id, result)
      }
    }
  }
  out$probabilities <- do.call(rbind, probabilities)
  rownames(out$probabilities) <- NULL
  out
}
