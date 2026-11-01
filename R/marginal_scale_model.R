# The private likelihood interface shares only neutral design preparation.
# Fixed numerical starts and GLM initialization do not attach population priors.
# The existing response, scale and zero-SD contracts remain unchanged.
.msr_marginal_prepare <- function(formula, inputs, data, scales = NULL,
                                  family = c("nbinom2", "poisson", "binomial", "gaussian"),
                                  row_id = NULL, reference_scale = NULL) {
  family <- match.arg(family)
  if (is.null(reference_scale)) {
    validate_scalar_numeric(inputs$min_D, "min_D", positive = TRUE)
    validate_scalar_numeric(inputs$max_D, "max_D", positive = TRUE)
    reference_scale <- exp(mean(log(c(inputs$min_D, inputs$max_D))))
  }
  validate_scalar_numeric(reference_scale, "reference_scale", positive = TRUE)
  spec <- .msr_scale_prepare(formula, inputs, data, scales, family, row_id,
    reference_scale = unname(reference_scale))
  .msr_scale_check_capability(spec, "marginal_likelihood", "laplace")
  spec$priors <- NULL
  spec$reference_scale <- unname(reference_scale)
  for (j in seq_along(spec$scales)) {
    sm <- spec$scales[[j]]
    spec$parameters$scale_beta[sm$beta_index] <- ifelse(sm$distance_location,
      log(reference_scale / spec$unit_conv), 0)
    sm$prior_mean <- sm$prior_sd <- NULL
    spec$scales[[j]] <- sm
  }
  # Fixed numerical starting values, not penalties on the fitted SDs.
  spec$parameters$log_scale_sd[] <- log(.2)
  spec$parameters$log_response_sd[] <- log(.2)
  spec$parameters$log_dispersion[] <- 0
  spec$zero_sd <- c(scale = FALSE, response = FALSE)
  spec$group_names <- c(
    scale = if (length(spec$scales[[1L]]$random)) spec$scales[[1L]]$random[[1L]]$group else "",
    response = if (length(spec$response$random)) spec$response$random[[1L]]$group else "")

  p <- spec$parameters
  kernel <- .msr_bayes_components(spec, p)$kernel
  X <- .msr_scale_response_matrix(spec, kernel)
  glm_family <- switch(family, gaussian = stats::gaussian(), binomial = stats::binomial(),
    stats::poisson())
  warnings <- character()
  initial <- tryCatch(withCallingHandlers(stats::glm.fit(X,
    if (family == "binomial") spec$y / spec$trials else spec$y,
    weights = if (family == "binomial") spec$trials else rep(1, length(spec$y)),
    offset = spec$offset, family = glm_family), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
    }), error = function(e) NULL)
  if (!is.null(initial) && all(is.finite(initial$coefficients)))
    spec$parameters$beta <- unname(initial$coefficients)
  spec$initialization <- list(method = if (is.null(initial)) "zero coefficients" else
    paste("reference-scale", glm_family$family, "GLM"), warnings = unique(warnings))
  spec$identification <- character()
  intercept <- attr(spec$response$terms, "intercept") == 1L
  if (intercept && family %in% c("poisson", "nbinom2") && all(spec$y == 0))
    spec$identification <- "No finite count-intercept MLE for an all-zero response."
  if (intercept && family == "binomial" && (all(spec$y == 0) || all(spec$y == spec$trials)))
    spec$identification <- "No finite binomial-intercept MLE for a constant outcome."
  if (family == "gaussian" && !is.null(initial)) {
    residual <- spec$y - drop(X %*% spec$parameters$beta) - spec$offset
    residual_sd <- sqrt(mean(residual^2))
    if (residual_sd <= 1e-10 * max(1, abs(spec$y)))
      spec$identification <- "Gaussian likelihood is unbounded when the mean fits every response exactly."
    spec$parameters$log_dispersion[] <- log(max(residual_sd, .01))
  }
  spec
}

# Exact zero-variance faces omit the corresponding latent block and its density.
# Keep the original representative rows so annual predictions still align.
.msr_marginal_reduce <- function(spec, zero = character()) {
  if (!is.character(zero) || anyDuplicated(zero) || any(!zero %in% c("scale", "response")))
    stop("zero must select distinct scale and/or response SD blocks.", call. = FALSE)
  for (block in zero) {
    field <- if (block == "scale") "log_scale_sd" else "log_response_sd"
    if (!length(spec$parameters[[field]]) && !isTRUE(spec$zero_sd[[block]]))
      stop("Cannot fix an absent SD block to zero.", call. = FALSE)
    if (block == "scale") {
      spec$scales[[1L]]$random <- list()
      spec$parameters$scale_z <- spec$parameters$log_scale_sd <- numeric()
    } else {
      spec$response$random <- list()
      spec$parameters$response_z <- spec$parameters$log_response_sd <- numeric()
    }
    spec$zero_sd[[block]] <- TRUE
  }
  spec
}

.msr_marginal_scale_predictions <- function(spec, parameters) {
  sm <- spec$scales[[1L]]; rows <- sm$representatives
  value <- drop(sm$X[rows, , drop = FALSE] %*% parameters$scale_beta[sm$beta_index])
  for (term in sm$random)
    value <- value + exp(parameters$log_scale_sd[term$sd_index]) *
      parameters$scale_z[term$z_index][term$index[rows]]
  value <- value + log(spec$unit_conv)
  data.frame(covariate = spec$covariates, group = sm$labels, row = rows,
    log_estimate = value, estimate = exp(value), interval_available = FALSE)
}
