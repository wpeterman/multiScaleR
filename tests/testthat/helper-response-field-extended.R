# One fixed outcome construction per family, reused across covariance layouts.
# Known range/shape isolate native prediction compatibility, not field recovery.
extended_response_field_fixture <- local({
  outcomes <- NULL; cache <- list()
  function(family = c("gaussian", "poisson", "binomial", "nbinom2", "Gamma"),
      covariance = c("exp", "gau", "mat")) {
    family <- match.arg(family); covariance <- match.arg(covariance)
    key <- paste(family, covariance, sep = "-")
    if (!is.null(cache[[key]])) return(cache[[key]])
    base <- native_response_field_fixture("gaussian")
    if (is.null(outcomes)) {
      d <- base$data; locations <- unique(cbind(d$x, d$y_coordinate))
      index <- match(d$field_pos, unique(d$field_pos))
      present <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
      old <- if (present) get(".Random.seed", .GlobalEnv) else NULL
      on.exit(if (present) assign(".Random.seed", old, .GlobalEnv) else
        if (exists(".Random.seed", .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
      set.seed(60261109)
      L <- t(chol(.6^2 * exp(-as.matrix(dist(locations)) / 55)))
      field <- drop(L %*% rnorm(nrow(locations)))
      eta <- .3 * d$mean - .2 * d$rough + .2 * d$z + log(d$exposure) + field[index]
      outcomes <<- list(gaussian = d$gaussian_response, poisson = d$poisson_response,
        binomial = rbinom(nrow(d), 30L, plogis(-.5 + eta)),
        nbinom2 = rnbinom(nrow(d), mu = exp(2 + eta), size = 12),
        Gamma = rgamma(nrow(d), shape = 20, scale = exp(1 + eta) / 20))
    }
    d <- base$data; d$y <- outcomes[[family]]
    if (family == "binomial") d$failure <- 30L - d$y
    response <- if (family == "binomial") "cbind(y, failure)" else "y"
    form <- as.formula(paste0(response,
      " ~ mean + rough + z + offset(log(exposure)) + ", covariance, "(field_pos + 0 | field_group)"))
    provider_family <- switch(family, gaussian = gaussian(), poisson = poisson(),
      binomial = binomial(), nbinom2 = glmmTMB::nbinom2(), Gamma = Gamma(link = "log"))
    theta <- c(log(.6), log(55), if (covariance == "mat") log(1.5))
    mapping <- factor(c(1L, rep(NA_integer_, length(theta) - 1L)))
    warnings <- character(); audit <- gdr_test_audit_path(paste0("response-field-extended-", key, "-"))
    fit <- tryCatch(withCallingHandlers(glmmTMB::glmmTMB(form, data = d,
      family = provider_family, start = list(theta = theta), map = list(theta = mapping),
      control = glmmTMB::glmmTMBControl(parallel = 1L, rank_check = "stop",
        optCtrl = list(iter.max = 1000, eval.max = 1500))),
      warning = function(w) warnings <<- c(warnings, conditionMessage(w))),
      error = function(e) {
        .msr_gdr_atomic_audit(list(status = "native_constructor_failed", family = family,
          covariance = covariance, seed = 60261109L, warnings = warnings, error = conditionMessage(e)), paste0(audit, ".rds"))
        stop(e)
      })
    .msr_gdr_atomic_audit(list(status = "native_constructor_returned", family = family,
      covariance = covariance, seed = 60261109L, warnings = warnings,
      convergence = fit$fit$convergence, pdHess = fit$sdr$pdHess,
      logLik = as.numeric(logLik(fit)), K = attr(logLik(fit), "df"),
      optimizer = fit$fit[c("par", "objective", "convergence", "message", "evaluations")],
      data_sha256 = .msr_gdr_hash(d), range_fixed = 55,
      matern_shape_fixed = if (covariance == "mat") 1.5 else NULL), paste0(audit, ".rds"))
    object <- base$object; object$opt_mod <- fit
    cache[[key]] <<- list(model = fit, object = object, data = d,
      inputs = base$inputs, par = base$par, points = base$points,
      rasters = base$rasters, range_fixed = 55)
    cache[[key]]
  }
})
