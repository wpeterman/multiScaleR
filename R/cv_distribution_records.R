.msr_cv_distribution_model <- function(fit) .analysis_model(if (inherits(fit, "multiScaleR")) fit$opt_mod else fit)

.msr_cv_distribution_parameters <- function(fit, control = NULL) {
  m <- .msr_cv_distribution_model(fit)
  if (!inherits(m, "glm") || inherits(m, "gam") || !is.data.frame(m$model) ||
      !isTRUE(m$converged) || !is.null(m$na.action) || !is.null(m$call$subset) || !is.null(m$call$offset) ||
      !is.null(m$th.warn) || !identical(m$family$link, "log"))
    stop("Require a converged, complete, native log-link NB/Gamma training GLM without theta warnings or external offsets.", call. = FALSE)
  exact_poisson <- inherits(fit, "multiScaleR") && identical(fit$backend, "joint_nb2") &&
    isTRUE(fit$dispersion_est$exact_poisson) && isTRUE(fit$dispersion_est$poisson_limit) &&
    identical(as.numeric(fit$dispersion_est$alpha), 0)
  family <- if (grepl("^Negative Binomial", m$family$family) ||
      (exact_poisson && identical(m$family$family, "poisson"))) "nbinom2" else m$family$family
  if (!.msr_cv_distribution_family(family)) stop("Native family is not an NB2/Gamma CV target.", call. = FALSE)
  y <- stats::model.response(m$model); w <- stats::model.weights(m$model)
  if (is.null(w)) w <- rep(1, nrow(m$model))
  dat <- data.frame(outcome = y, row.names = rownames(m$model))
  observed <- .msr_cv_distribution_outcome(outcome ~ 1, dat, family)$observed
  if (!is.numeric(w) || any(!is.finite(w) | w < 0 | w != floor(w)) || !any(w > 0) || !is.finite(sum(w)))
    stop("Native likelihood weights must be finite nonnegative integer frequencies with positive total.", call. = FALSE)
  mu <- stats::fitted(m)
  if (length(mu) != length(y) || any(!is.finite(mu) | mu <= 0) ||
      !is.numeric(m$rank) || length(m$rank) != 1L || !is.finite(m$rank) || m$rank < 0 ||
      any(is.nan(stats::coef(m))) || any(!is.finite(stats::coef(m)[!is.na(stats::coef(m))])))
    stop("Native training mean, estimable coefficients or rank are invalid.", call. = FALSE)
  theta <- dispersion <- NA_real_; poisson_limit <- FALSE
  if (family == "nbinom2") {
    estimated <- inherits(m, "negbin")
    env <- environment(m$family$aic)
    if (exact_poisson && identical(m$family$family, "poisson")) native_theta <- Inf else {
      if (!is.environment(env) || !exists(".Theta", env, inherits = FALSE))
        stop("Cannot recover the native negative-binomial size parameter.", call. = FALSE)
      native_theta <- get(".Theta", env, inherits = FALSE)
    }
    theta <- if (estimated) m$theta else native_theta
    poisson_limit <- exact_poisson
    if (!is.numeric(theta) || length(theta) != 1L || is.na(theta) || theta <= 0 ||
        (!is.finite(theta) && !(is.infinite(theta) && poisson_limit)) ||
        !is.numeric(native_theta) || length(native_theta) != 1L || is.na(native_theta) || native_theta <= 0 ||
        (!is.finite(native_theta) && !poisson_limit) ||
        (!estimated && !identical(as.numeric(theta), as.numeric(native_theta))) ||
        (poisson_limit && (is.finite(theta) || !identical(as.numeric(theta), as.numeric(native_theta)))))
      stop("NB2 theta must match its native family; infinity requires the original verified joint-NB2 Poisson-limit flag.", call. = FALSE)
    # MASS can retain the last IRLS family's size while its final theta and
    # twologlik use the following alternation. Preserve that native fit/SE;
    # the normalized density below independently verifies its final likelihood.
    convention <- if (estimated) "estimated_nb2_theta" else "fixed_nb2_theta"
  } else {
    fixed <- m$family$msr_cv_fixed_dispersion
    estimated <- is.null(fixed)
    if (estimated) {
      if (!is.null(m$family$dispersion) && !is.na(m$family$dispersion))
        stop("An unverified fixed Gamma dispersion field does not establish a normalized likelihood.", call. = FALSE)
      dispersion <- m$deviance/sum(w)
      convention <- "native_gamma_deviance_scale"
    } else {
      dispersion <- fixed
      if (!identical(as.numeric(m$family$dispersion), as.numeric(fixed)))
        stop("Known Gamma dispersion metadata disagree.", call. = FALSE)
      convention <- "fixed_gamma_dispersion"
    }
    if (!is.numeric(dispersion) || length(dispersion) != 1L || !is.finite(dispersion) || dispersion <= 0)
      stop("Gamma training dispersion must be finite and strictly positive.", call. = FALSE)
  }
  if (!is.null(control)) {
    if (!identical(control$family, family) ||
        (family == "nbinom2" && (estimated != is.null(control$theta) ||
          (!estimated && !identical(as.numeric(theta), control$theta)))) ||
        (family == "Gamma" && (estimated != is.null(control$dispersion) ||
          (!estimated && !identical(as.numeric(dispersion), control$dispersion)))))
      stop("Native fixed/estimated nuisance target disagrees with CV controls.", call. = FALSE)
  }
  used <- w > 0
  density <- if (family == "Gamma") stats::dgamma(observed[used], shape = 1/dispersion,
    scale = mu[used] * dispersion, log = TRUE) else if (poisson_limit) stats::dpois(observed[used], mu[used], log = TRUE) else
      stats::dnbinom(observed[used], mu = mu[used], size = theta, log = TRUE)
  direct <- sum(w[used] * density); ll <- .msr_model_loglik(m); k <- .msr_parameter_count(m)
  if (!is.finite(direct) || !is.finite(ll) ||
      !isTRUE(all.equal(as.numeric(ll), direct, tolerance = 1e-8, scale = 1)) || k != m$rank + as.integer(estimated))
    stop("Native NB/Gamma likelihood or nuisance parameter count disagrees with its normalized training density.", call. = FALSE)
  list(family = family, theta = as.numeric(theta), dispersion = as.numeric(dispersion), poisson_limit = poisson_limit,
    estimated = estimated, convention = convention, rank = as.integer(m$rank), native_K = k,
    logLik = as.numeric(ll), nobs = .msr_selection_nobs(m), training_deviance = m$deviance, sum_weights = sum(w))
}

.msr_cv_distribution_training_identity <- function(fit, split, targets, family) {
  m <- .msr_cv_distribution_model(fit)
  if (!is.data.frame(m$model) || !identical(rownames(m$model), split$analysis))
    stop("Retained distribution GLMs must match every analysis row ID; full-data fits cannot replace training fits.", call. = FALSE)
  p <- .msr_cv_distribution_parameters(fit)
  if (!identical(p$family, family)) stop("Native distribution family disagrees with the held-out target.", call. = FALSE)
  dat <- data.frame(outcome = stats::model.response(m$model), row.names = rownames(m$model))
  y <- .msr_cv_distribution_outcome(outcome ~ 1, dat, family)
  w <- stats::model.weights(m$model); if (is.null(w)) w <- rep(1, nrow(dat))
  expected <- targets[match(split$analysis, targets$row_id), , drop = FALSE]
  if (!identical(y$observed, expected$observed) || !identical(y$trials, expected$trials) ||
      !identical(as.numeric(w), expected$weight)) stop("Native training outcomes or frequency weights disagree with CV targets.", call. = FALSE)
  invisible(TRUE)
}

.msr_cv_distribution_record <- function(fit, diagnostics, controls) {
  out <- list(version = 3L, controls = controls, diagnostics = diagnostics); messages <- character()
  out$data <- tryCatch(withCallingHandlers({
    m <- .msr_cv_distribution_model(fit)
    nuisance <- .msr_cv_distribution_parameters(fit, controls$distribution_control)
    w <- stats::model.weights(m$model); if (is.null(w)) w <- rep(1, nrow(m$model))
    list(logLik = nuisance$logLik, native_K = nuisance$native_K,
      ecological_K = .msr_selection_parameter_count(fit, 0L), K = .msr_selection_parameter_count(fit, nuisance$native_K),
      nobs = nuisance$nobs, observation_ids = .msr_observation_ids(m), nuisance = nuisance,
      targets = data.frame(row_id = rownames(m$model), observed = as.numeric(stats::model.response(m$model)),
        trials = 1, weight = as.numeric(w), stringsAsFactors = FALSE),
      signature = .msr_selection_signature(m), formula = paste(deparse(stats::formula(m)), collapse = " "),
      family = nuisance$family, link = "log")
  }, warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = function(e) { out$error <<- conditionMessage(e); NULL })
  out$warnings <- messages; out
}

.msr_cv_distribution_record_identity <- function(record, x, name, targets) {
  s <- .msr_cv_fold_settings(x, name)
  if (!is.list(record) || !identical(record$version, 3L) || !is.character(record$warnings) ||
      !identical(record$controls, .msr_cv_selection_controls(s)) || !identical(record$diagnostics, x$diagnostics[[name]]))
    stop("Invalid distribution selection-record schema, controls or diagnostics.", call. = FALSE)
  if (!is.null(record$error)) {
    if (!is.character(record$error) || length(record$error) != 1L || is.na(record$error) || !is.null(record$data))
      stop("Malformed distribution metadata failure.", call. = FALSE)
    return(invisible(TRUE))
  }
  d <- record$data; p <- d$nuisance; expected <- targets[match(x$folds[[name]]$analysis, targets$row_id), , drop = FALSE]
  rownames(expected) <- NULL
  scalar <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v)
  if (!is.list(d) || !identical(d$family, s$family) || !identical(d$link, "log") || !identical(d$targets, expected) ||
      !is.list(p) || !identical(p$family, s$family) || !is.logical(p$estimated) || length(p$estimated) != 1L || is.na(p$estimated))
    stop("Distribution record family, nuisance metadata or training targets disagree.", call. = FALSE)
  for (nm in c("K", "native_K", "ecological_K", "nobs"))
    if (!scalar(d[[nm]]) || d[[nm]] < 0 || d[[nm]] != floor(d[[nm]])) stop("Invalid native distribution parameter count or nobs.", call. = FALSE)
  if (!scalar(d$logLik) || !scalar(p$rank) || p$rank < 0 || p$rank != floor(p$rank) ||
      d$native_K != p$rank + as.integer(p$estimated) || d$K != d$native_K + d$ecological_K ||
      !identical(p$logLik, d$logLik) || p$native_K != d$native_K || p$nobs != d$nobs ||
      !is.character(d$formula) || length(d$formula) != 1L || is.na(d$formula) || !nzchar(d$formula))
    stop("Distribution record native likelihood or nuisance count is invalid.", call. = FALSE)
  .msr_cv_distribution_metadata_check(p, s$distribution_control)
  if (p$poisson_limit && !identical(s$backend, "joint_nb2"))
    stop("Recorded Poisson limit requires the verified joint-NB2 backend.", call. = FALSE)
  if (!scalar(p$sum_weights) || p$sum_weights != sum(expected$weight) ||
      (p$family == "Gamma" && p$estimated && (!scalar(p$training_deviance) ||
        !identical(p$dispersion, p$training_deviance/p$sum_weights))))
    stop("Distribution training nuisance convention disagrees with native deviance or frequency total.", call. = FALSE)
  used <- expected$weight > 0
  if (d$nobs != sum(used) || !identical(d$observation_ids, expected$row_id[used]))
    stop("Distribution native effective rows or observation count disagree with targets.", call. = FALSE)
  e <- expected[used, , drop = FALSE]; e <- e[order(e$row_id), , drop = FALSE]
  signature <- list(response = e$observed, weights = e$weight, measure = if (s$family == "Gamma") "continuous" else "discrete",
    trials = NULL, reml = FALSE, design = NULL, partial_method = NULL, strata = NULL)
  if (!isTRUE(all.equal(d$signature, signature, tolerance = 0, check.attributes = FALSE)))
    stop("Distribution likelihood signature disagrees with training outcomes and weights.", call. = FALSE)
  invisible(TRUE)
}

.msr_cv_distribution_metadata_check <- function(p, control) {
  scalar <- function(v) is.numeric(v) && length(v) == 1L && is.finite(v) && v > 0
  missing_scalar <- function(v) is.numeric(v) && length(v) == 1L && is.na(v) && !is.nan(v)
  if (!is.list(p) || !identical(p$family, control$family) || !is.logical(p$estimated) || length(p$estimated) != 1L || is.na(p$estimated) ||
      !is.logical(p$poisson_limit) || length(p$poisson_limit) != 1L || is.na(p$poisson_limit))
    stop("Invalid distribution nuisance target metadata.", call. = FALSE)
  if (p$family == "Gamma") {
    if (!scalar(p$dispersion) || !missing_scalar(p$theta) || p$poisson_limit || p$estimated != is.null(control$dispersion) ||
        !identical(p$convention, if (p$estimated) "native_gamma_deviance_scale" else "fixed_gamma_dispersion") ||
        (!p$estimated && !identical(p$dispersion, control$dispersion))) stop("Gamma fixed/estimated nuisance metadata disagree.", call. = FALSE)
  } else if ((!scalar(p$theta) && !(is.numeric(p$theta) && length(p$theta) == 1L && is.infinite(p$theta) && p$theta > 0 && p$poisson_limit)) ||
      !missing_scalar(p$dispersion) || (p$poisson_limit && (!p$estimated || !is.infinite(p$theta))) ||
      p$estimated != is.null(control$theta) ||
      !identical(p$convention, if (p$estimated) "estimated_nb2_theta" else "fixed_nb2_theta") ||
      (!p$estimated && !identical(p$theta, control$theta))) stop("NB2 fixed/estimated nuisance metadata disagree.", call. = FALSE)
  invisible(TRUE)
}
