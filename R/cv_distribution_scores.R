.msr_cv_distribution_score <- function(p, metrics) {
  required <- c("distribution_family", "theta", "dispersion", "poisson_limit")
  if (!all(required %in% names(p)) || anyNA(p$distribution_family) ||
      any(!p$distribution_family %in% c("nbinom2", "Gamma", "poisson", "gaussian")) ||
      !is.logical(p$poisson_limit) || anyNA(p$poisson_limit) || any(!is.finite(p$prediction)) ||
      any(!is.finite(p$observed)) || any(!is.finite(p$weight) | p$weight < 0 | p$weight != floor(p$weight)) ||
      any(p$trials != 1)) stop("Invalid recorded distribution predictions, observations or case frequencies.", call. = FALSE)
  nb <- p$distribution_family == "nbinom2"; ga <- p$distribution_family == "Gamma"
  po <- p$distribution_family == "poisson"; normal <- p$distribution_family == "gaussian"
  if (any((nb | po) & (p$observed < 0 | p$observed != floor(p$observed) | p$prediction < 0)) ||
      any(ga & (p$observed <= 0 | p$prediction <= 0)) || any(normal & p$weight != 1) ||
      any(ga & (!is.finite(p$dispersion) | p$dispersion <= 0)) ||
      any(nb & (is.na(p$theta) | p$theta <= 0 | (!is.finite(p$theta) & !(is.infinite(p$theta) & p$poisson_limit)))) ||
      any(p$poisson_limit & (!nb | !is.infinite(p$theta))) ||
      any(!nb & !is.na(p$theta)) || any(!ga & !is.na(p$dispersion)))
    stop("Distribution means, responses or training-only nuisance parameters are outside their declared support.", call. = FALSE)
  eligible <- p$weight > 0
  if (!any(eligible)) stop("The assessment set has no positive-weight observations.", call. = FALSE)
  p$eligible <- eligible
  for (metric in metrics) {
    value <- switch(metric, rmse = (p$observed - p$prediction)^2, mae = abs(p$observed - p$prediction),
      neg_log_score = {
        if (any(normal)) stop("Gaussian density scoring requires a separately qualified training dispersion.", call. = FALSE)
        logp <- rep(NA_real_, nrow(p))
        if (any(po | p$poisson_limit)) {
          take <- po | p$poisson_limit
          logp[take] <- stats::dpois(p$observed[take], p$prediction[take], log = TRUE)
        }
        take <- nb & !p$poisson_limit
        if (any(take)) logp[take] <- stats::dnbinom(p$observed[take], size = p$theta[take], mu = p$prediction[take], log = TRUE)
        if (any(ga)) logp[ga] <- stats::dgamma(p$observed[ga], shape = 1/p$dispersion[ga],
          scale = p$prediction[ga] * p$dispersion[ga], log = TRUE)
        -logp
      }, stop("Distribution CV supports rmse, mae and neg_log_score.", call. = FALSE))
    if (any(is.nan(value) | (!is.na(value) & value == -Inf))) stop("Invalid held-out distribution loss.", call. = FALSE)
    value[!eligible] <- NA_real_; p[[metric]] <- value
  }
  p
}

# Canonical reference rows acquire metadata only inside a distribution comparison
# or selection result. Their original returned CV objects remain untouched.
.msr_cv_distribution_rows <- function(x, name) {
  s <- .msr_cv_fold_settings(x, name)
  p <- x$predictions[x$predictions$fold == name, , drop = FALSE]
  fields <- c("distribution_family", "theta", "dispersion", "poisson_limit")
  if (!any(fields %in% names(p)) && s$family %in% c("poisson", "gaussian")) {
    p$distribution_family <- s$family; p$theta <- p$dispersion <- NA_real_; p$poisson_limit <- FALSE
  }
  if (!all(fields %in% names(p)) || !is.character(p$distribution_family) ||
      anyNA(p$distribution_family) || any(p$distribution_family != s$family))
    stop("Recorded prediction families disagree with the actual training-fold settings.", call. = FALSE)
  if (.msr_cv_distribution_family(s$family)) {
    d <- x$distribution[[name]]
    .msr_cv_distribution_metadata_check(d, s$distribution_control)
    if (!identical(as.numeric(p$theta), rep(d$theta, nrow(p))) ||
        !identical(as.numeric(p$dispersion), rep(d$dispersion, nrow(p))) ||
        !identical(p$poisson_limit, rep(d$poisson_limit, nrow(p))) ||
        (d$poisson_limit && !identical(s$backend, "joint_nb2")))
      stop("Assessment nuisance parameters disagree with the frozen native training distribution.", call. = FALSE)
  } else if (any(!is.na(p$theta)) || any(!is.na(p$dispersion)) || any(p$poisson_limit))
    stop("Canonical reference rows carry unqualified distribution parameters.", call. = FALSE)
  p
}

.msr_cv_distribution_metadata <- function(x, name) {
  s <- .msr_cv_fold_settings(x, name)
  if (.msr_cv_distribution_family(s$family)) return(x$distribution[[name]])
  list(family = s$family, theta = NA_real_, dispersion = NA_real_, poisson_limit = FALSE,
    estimated = FALSE, convention = "not_required_for_mean_or_poisson_density")
}

.msr_cv_distribution_comparison <- function(mod_list, metric = NULL) {
  signatures <- lapply(mod_list, function(x) {
    if (!inherits(x, "multiScaleR_cv") || inherits(x, "multiScaleR_choice_cv") || !isTRUE(x$complete) ||
        !is.data.frame(x$fold_results) || any(x$fold_results$status != "scored") ||
        !identical(x$settings$prediction_target, "held_out_GLM_expected_response") ||
        !identical(x$settings$parameter_uncertainty, "not_integrated") ||
        !identical(x$settings$observation_intervals, "not_provided"))
      stop("Compare complete qualified response-distribution CV targets; failed folds cannot be dropped.", call. = FALSE)
    if (!is.null(metric) && !metric %in% x$settings$metrics) stop("Every candidate must record the requested metric.", call. = FALSE)
    used <- sort(unique(unlist(lapply(x$folds, function(f) c(f$analysis, f$assessment)), use.names = FALSE)))
    .msr_cv_folds(x$folds, used)
    t <- x$targets; columns <- c("row_id", "observed", "trials", "weight")
    if (!is.data.frame(t) || !all(columns %in% names(t)) || anyNA(t$row_id) || anyDuplicated(t$row_id) ||
        !setequal(t$row_id, used) || !all(vapply(t[columns[-1L]], is.numeric, logical(1))) ||
        any(!is.finite(t$observed)) || anyNA(t$trials) || any(t$trials != 1) ||
        any(!is.finite(t$weight) | t$weight < 0 | t$weight != floor(t$weight)))
      stop("Distribution targets must retain every used response, unit trial and integer frequency.", call. = FALSE)
    t <- t[match(used, t$row_id), columns, drop = FALSE]; rownames(t) <- NULL
    folds <- lapply(x$folds[sort(names(x$folds))], function(f) lapply(f[c("analysis", "assessment")], sort))
    attributes(folds) <- list(names = names(folds))
    p <- x$predictions
    if (!is.data.frame(p) || !all(c(columns, "fold", "prediction", "eligible") %in% names(p)) ||
        anyNA(p$fold) || any(!p$fold %in% names(folds))) stop("Invalid recorded distribution assessment rows.", call. = FALSE)
    measures <- character()
    for (name in names(folds)) {
      s <- .msr_cv_fold_settings(x, name)
      if (!s$family %in% c("nbinom2", "Gamma", "poisson", "gaussian") ||
          !identical(s$link, if (s$family == "gaussian") "identity" else "log"))
        stop("This distribution comparison requires qualified Gaussian, Poisson, NB2 or Gamma targets.", call. = FALSE)
      rows <- p[p$fold == name, , drop = FALSE]
      if (anyNA(rows$row_id) || anyDuplicated(rows$row_id) || !setequal(rows$row_id, folds[[name]]$assessment))
        stop("Recorded distribution assessment rows disagree with the split plan.", call. = FALSE)
      expected <- t[match(rows$row_id, t$row_id), columns, drop = FALSE]
      actual <- rows[columns]; rownames(actual) <- rownames(expected) <- NULL
      if (!identical(actual, expected)) stop("Recorded distribution outcomes, trials or weights disagree with targets.", call. = FALSE)
      analysis <- t[match(x$folds[[name]]$analysis, t$row_id), , drop = FALSE]
      values <- data.frame(outcome = c(analysis$observed, rows$observed))
      .msr_cv_outcome(outcome ~ 1, values, s$family)
      if (s$family == "gaussian" && any(c(analysis$weight, rows$weight) != 1))
        stop("Gaussian references retain their unit-weight contract.", call. = FALSE)
      rows <- .msr_cv_distribution_rows(x, name)
      score <- if (is.null(metric)) x$settings$metrics else metric
      check <- .msr_cv_distribution_score(rows, score)
      if (!identical(rows$eligible, check$eligible)) stop("Recorded distribution eligibility is inconsistent.", call. = FALSE)
      for (field in score) if (!isTRUE(all.equal(rows[[field]], check[[field]], tolerance = 0, check.attributes = FALSE)))
        stop("Recorded distribution losses disagree with fitted means and training nuisance parameters.", call. = FALSE)
      measures[name] <- if (s$family %in% c("poisson", "nbinom2")) "discrete" else "continuous"
    }
    list(folds = folds, targets = t, measure = measures)
  })
  first <- signatures[[1L]]
  for (s in signatures) {
    if (!identical(first$folds, s$folds) || !identical(first$targets, s$targets))
      stop("Distribution comparison requires identical splits, observations and case frequencies.", call. = FALSE)
    if ((is.null(metric) || !metric %in% c("rmse", "mae")) && !identical(first$measure, s$measure))
      stop("Information criteria and log scores require a common continuous/discrete observation measure.", call. = FALSE)
  }
  invisible(first)
}
