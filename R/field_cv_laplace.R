# Likelihood contributions and eta derivatives for one observation at one node.
.msr_field_response <- function(family, y, eta, dispersion) {
  if (family == "poisson") {
    mean <- exp(eta)
    return(list(loglik = y * eta - mean - lgamma(y + 1),
      score = y - mean, information = mean))
  }
  if (family == "binomial") {
    p <- stats::plogis(eta)
    return(list(loglik = .msr_scale_response_loglik("binomial", y, eta, NULL,
      trials = 1), score = y - p, information = p * stats::plogis(-eta)))
  }
  if (family == "nbinom2") {
    eta_size <- eta - log(dispersion)
    p <- stats::plogis(eta_size)
    q <- stats::plogis(-eta_size)
    return(list(loglik = .msr_scale_response_loglik("nbinom2", y, eta,
      log(dispersion)), score = y * q - dispersion * p,
      information = (y + dispersion) * p * q))
  }
  relative <- exp(log(y) - eta)
  list(loglik = dispersion * (log(dispersion) - eta - relative) +
      (dispersion - 1) * log(y) - lgamma(dispersion),
    score = dispersion * (relative - 1), information = dispersion * relative)
}

# Whiten the field so the normalized prior is N(0, I). Its (2*pi) factor
# cancels the Laplace normalizer, including when dimensions differ by block.
.msr_field_laplace <- function(y, eta, group, lower, family, dispersion,
    start = numeric(ncol(lower))) {
  p <- ncol(lower)
  if (!p) {
    value <- sum(.msr_field_response(family, y, eta, dispersion)$loglik)
    if (!is.finite(value)) stop("The zero-field response log likelihood is not finite.", call. = FALSE)
    return(list(log_density = value, mode = numeric(), chol = lower,
      iterations = 0L, max_score = 0))
  }
  at <- function(z) {
    response <- .msr_field_response(family, y, eta + drop(lower %*% z)[group], dispersion)
    list(value = sum(z^2)/2 - sum(response$loglik), response = response)
  }
  z <- start
  point <- at(z)
  tolerance <- 1e-7 * max(1, sqrt(length(y)))
  for (iteration in seq_len(200L)) {
    if (!is.finite(point$value))
      stop("The latent field objective is not finite at its starting point.", call. = FALSE)
    grouped <- rowsum(cbind(point$response$score, point$response$information), group)
    gradient <- z - drop(crossprod(lower, grouped[, 1L]))
    H <- crossprod(lower, lower * grouped[, 2L])
    diag(H) <- diag(H) + 1
    curvature <- chol(H)
    score <- max(abs(gradient))
    if (score <= tolerance) {
      return(list(log_density = -point$value - sum(log(diag(curvature))),
        mode = z, chol = curvature, iterations = iteration - 1L, max_score = score))
    }
    step <- backsolve(curvature, forwardsolve(t(curvature), gradient))
    decrease <- sum(gradient * step)
    accepted <- FALSE
    for (j in 0:39) {
      fraction <- 2^-j
      candidate <- z - fraction * step
      next_point <- at(candidate)
      if (is.finite(next_point$value) && next_point$value <= point$value -
          1e-4 * fraction * decrease + 8 * .Machine$double.eps * max(1, abs(point$value))) {
        accepted <- TRUE
        break
      }
    }
    if (!accepted)
      stop("The latent field Newton search could not improve its objective.", call. = FALSE)
    z <- candidate
    point <- next_point
  }
  stop("The latent field mode did not converge.", call. = FALSE)
}

.msr_field_response_moments <- function(family, eta, variance, dispersion) {
  if (family == "binomial") {
    integrate <- function(link, location, spread) {
      if (spread == 0) return(link(location))
      stats::integrate(function(z) link(location + spread * z) * stats::dnorm(z),
        lower = -Inf, upper = Inf, rel.tol = 1e-9, abs.tol = 1e-11)$value
    }
    spread <- sqrt(variance)
    # Integrate the smaller tail before reconstructing probabilities near one.
    tail <- mapply(integrate, location = -abs(eta), spread = spread,
      MoreArgs = list(link = stats::plogis))
    mean <- ifelse(eta > 0, 1 - tail, tail)
    noise <- mapply(integrate, location = eta, spread = spread,
      MoreArgs = list(link = function(x) stats::plogis(x) * stats::plogis(-x)))
    return(list(mean = mean, noise = noise, variance = tail * (1 - tail)))
  }
  .msr_field_log_response_moments(family, eta + variance/2,
    2 * eta + 2 * variance, dispersion, log_ratio = -variance)
}

.msr_field_log_response_moments <- function(family, log_mean, log_second, dispersion,
    log_ratio = 2 * log_mean - log_second) {
  roundoff <- 128 * .Machine$double.eps * pmax(1, abs(log_second))
  if (any(!is.finite(log_ratio)) || any(log_ratio > roundoff))
    stop("The integrated second response-mean moment is smaller than the square of the mean.", call. = FALSE)
  mean <- exp(log_mean)
  # Apply dispersion and the mixing fraction before exponentiating so that
  # an unrepresentable raw second moment need not overflow a finite variance.
  noise <- switch(family, poisson = mean,
    nbinom2 = mean + exp(log_second - log(dispersion)),
    Gamma = exp(log_second - log(dispersion)))
  mixing <- exp(log_second + log(-expm1(pmin(0, log_ratio))))
  list(mean = mean, noise = noise, variance = noise + mixing)
}

# Completing the normal square tilts the analysis prior by k times the
# assessment projection for moment order k. Move the preceding mode into
# those shifted coordinates before solving with the same covariance factor.
.msr_field_tilted_log_moments <- function(y, eta, group, lower, training_nodes,
    assessment_nodes, assessment_eta, family, dispersion, analysis) {
  training_lower <- lower[training_nodes, training_nodes, drop = FALSE]
  nodes <- unique(assessment_nodes)
  ratios <- vapply(nodes, function(node) {
    tilt <- lower[node, training_nodes]
    offset <- drop(training_lower %*% tilt)[group]
    mode <- analysis$mode
    ratio <- numeric(2L)
    for (order in 1:2) {
      shifted <- .msr_field_laplace(y, eta + order * offset, group, training_lower,
        family, dispersion, start = mode - tilt)
      ratio[order] <- shifted$log_density - analysis$log_density
      mode <- shifted$mode
    }
    ratio
  }, numeric(2L))
  shifted <- t(ratios[, match(assessment_nodes, nodes), drop = FALSE])
  field_variance <- rowSums(lower[assessment_nodes, , drop = FALSE]^2)
  cbind(assessment_eta + field_variance/2,
    2 * assessment_eta + 2 * field_variance) + shifted
}

# Coordinate order puts analysis nodes first. The leading Cholesky block is
# their prior, and the off-diagonal block projects their whitened coordinates.
.msr_field_laplace_block <- function(train_y, test_y, train_eta, test_eta,
    train_xy, test_xy, covariance, field_sd, range, family, dispersion,
    tilted_log_moments = FALSE) {
  all_xy <- rbind(train_xy, test_xy)
  position <- as.character(glmmTMB::numFactor(all_xy))
  first <- !duplicated(position)
  nodes <- all_xy[first, , drop = FALSE]
  group <- match(position, position[first])
  training <- seq_along(train_y)
  test <- length(train_y) + seq_along(test_y)
  if (nrow(unique(unname(all_xy))) != nrow(nodes) || any(group[test] %in% group[training]))
    stop("Use disjoint analysis and assessment locations with exact native coordinate encoding.", call. = FALSE)
  field <- field_sd > 0
  if (field) {
    power <- if (covariance == "exponential") 1 else 2
    lower <- field_sd * t(chol(exp(-(as.matrix(stats::dist(nodes))/range)^power)))
    training_nodes <- seq_len(max(group[training]))
    new_nodes <- setdiff(seq_len(nrow(nodes)), training_nodes)
    training_lower <- lower[training_nodes, training_nodes, drop = FALSE]
  } else {
    lower <- training_lower <- matrix(numeric(), 0L, 0L)
  }
  analysis <- .msr_field_laplace(train_y, train_eta, group[training], training_lower,
    family, dispersion)
  training_variance <- new_variance <- deviation <- numeric(length(test))
  if (field) {
    projection <- lower[group[test], training_nodes, drop = FALSE]
    deviation <- drop(projection %*% analysis$mode)
    loading <- forwardsolve(t(analysis$chol), t(projection))
    training_variance <- colSums(loading^2)
    new_variance <- rowSums(lower[group[test], new_nodes, drop = FALSE]^2)
  }
  combined <- .msr_field_laplace(c(train_y, test_y), c(train_eta, test_eta), group, lower,
    family, dispersion, start = c(analysis$mode, numeric(ncol(lower) - length(analysis$mode))))
  loss <- analysis$log_density - combined$log_density
  if (family %in% c("poisson", "nbinom2", "binomial") && loss < 0) {
    rounding <- 128 * .Machine$double.eps *
      max(1, abs(analysis$log_density), abs(combined$log_density))
    if (loss < -rounding)
      stop("The Laplace ratio produced an impossible discrete block probability; use a smaller assessment block or inspect the fitted field.", call. = FALSE)
    loss <- 0
  }
  moments <- if (field && tilted_log_moments && family != "binomial") {
    log_moments <- .msr_field_tilted_log_moments(train_y, train_eta, group[training], lower,
      training_nodes, group[test], test_eta, family, dispersion, analysis)
    .msr_field_log_response_moments(family, log_moments[, 1L], log_moments[, 2L], dispersion)
  } else {
    .msr_field_response_moments(family, test_eta + deviation,
      training_variance + new_variance, dispersion)
  }
  if (any(!is.finite(c(loss, moments$mean, moments$noise, moments$variance,
      training_variance, new_variance))) || any(moments$variance < 0))
    stop("Non-Gaussian field block moments or score are not finite and valid.", call. = FALSE)
  list(prediction = as.numeric(moments$mean), fitted_training_field_variance = training_variance,
    new_location_field_variance = new_variance, observation_noise_variance = as.numeric(moments$noise),
    observation_predictive_SE = as.numeric(sqrt(moments$variance)),
    joint = list(n = length(test), method = "Laplace_marginal_likelihood_ratio",
      analysis_log_density = analysis$log_density, combined_log_density = combined$log_density,
      analysis_latent_score = analysis$max_score, combined_latent_score = combined$max_score,
      analysis_iterations = analysis$iterations, combined_iterations = combined$iterations,
      neg_log_score = loss))
}

.msr_field_cv_predict_laplace <- function(model, formula, training, assessment,
    train_xy, test_xy, covariance, max_memory_mb, family) {
  basis <- .msr_field_cv_basis(formula, training, model)
  tr <- .msr_unmarked_fixed_design(basis, training)
  te <- .msr_unmarked_fixed_design(basis, assessment)
  if (any(!is.finite(te$X)))
    stop("Assessment fixed designs must be complete and finite.", call. = FALSE)
  n_train <- nrow(unique(unname(train_xy))); n_test <- nrow(unique(unname(test_xy)))
  .msr_field_cv_budget(n_train + n_test, nrow(assessment), ncol(tr$X) + 3L, max_memory_mb)
  beta <- glmmTMB::fixef(model)$cond
  train_eta <- drop(tr$X %*% beta) + if (is.null(tr$offset)) 0 else tr$offset
  test_eta <- drop(te$X %*% beta) + if (is.null(te$offset)) 0 else te$offset
  parameters <- model$obj$env$parList(par = model$fit$parfull)
  theta <- parameters$theta
  if (length(theta) != 2L || any(!is.finite(theta)))
    stop("Field CV requires the original two-parameter covariance layout.", call. = FALSE)
  field_sd <- exp(theta[1L]); range <- exp(theta[2L])
  dispersion <- if (family == "nbinom2") stats::sigma(model) else
    if (family == "Gamma") stats::sigma(model)^-2 else 1
  if (any(!is.finite(c(train_eta, test_eta, field_sd, range, dispersion))) ||
      range <= 0 || dispersion <= 0)
    stop("Field CV requires finite predictors and representable covariance and dispersion parameters.", call. = FALSE)
  response <- as.character(formula[[2L]])
  result <- .msr_field_laplace_block(training[[response]], assessment[[response]], train_eta,
    test_eta, train_xy, test_xy, covariance, field_sd, range, family, dispersion,
    tilted_log_moments = family != "binomial")
  gap <- abs(result$joint$analysis_log_density - .msr_model_loglik(model))
  if (gap > max(1e-6, 1e-8 * abs(result$joint$analysis_log_density)))
    stop("The analysis marginal density disagrees with the normalized native training likelihood.", call. = FALSE)
  result$joint$native_analysis_log_density_gap <- gap
  result$parameters <- list(field_SD = field_sd, range = range,
    dispersion = if (family %in% c("nbinom2", "Gamma")) dispersion else NULL,
    dispersion_parameter = switch(family, nbinom2 = "size", Gamma = "shape", "none"),
    native_K = .msr_parameter_count(model), range_estimated = TRUE,
    estimated_range_parameters = 1L, integrated_state_penalty = 0L)
  result$recipe_sha256 <- .msr_gdr_hash(list(basis$terms, basis$xlevels, basis$contrasts, basis$columns))
  contract <- .msr_prediction_contract("native", "glmmTMB_Laplace_field_CV",
    "expected_response", "response", population = "fixed_at_training_ML_estimates",
    kernel = "fixed_at_training_fold_estimates",
    effects = list(response_field = "Laplace_approximate_integration_given_training_responses"),
    conditioning = c("analysis_responses_only", "fitted_training_parameters", "frozen_training_recipes"),
    integration = list(field = if (family == "binomial") "training_posterior_Laplace_Gaussian_moments" else
        "training_posterior_exponentially_tilted_Laplace_moments",
      assessment = "joint_Laplace_marginal_likelihood_ratio"),
    omitted = c("coefficient_parameter_uncertainty", "range_SD_parameter_uncertainty",
      "ecological_scale_parameter_uncertainty", "selection_uncertainty"),
    coordinates = "native_fitted_parameters_and_projected_locations", diagnostics = list(passed = TRUE),
    fixed_inputs = TRUE, component = paste0(family, "_joint_assessment_response"))
  contract$uncertainty$observation_noise <- TRUE
  contract$uncertainty$standard_errors <- "approximate_observation_predictive_SD_not_native_mean_SE"
  contract$response_spatial_field <- list(covariance = covariance, range_estimated = TRUE,
    family = family, parameter_uncertainty_integrated = FALSE, ecological_scale_uncertainty_integrated = FALSE,
    field_variance_scale = "linear_predictor", observation_variance_scale = "response",
    variance_components = c("remaining_training_field", "new_location_field", "observation_noise"),
    full_grid_covariance_returned = FALSE)
  result$prediction_contract <- contract
  result
}
