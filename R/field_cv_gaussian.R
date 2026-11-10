# Gaussian field scores condition on fitted parameters. Native mean-prediction
# SEs include a different uncertainty target and must not replace this density.
.msr_field_cv_budget <- function(nodes, assessment, columns, max_memory_mb) {
  bytes <- 8 * (64 * as.double(nodes)^2 + 16 * as.double(assessment)^2 +
    32 * as.double(nodes) * assessment + 32 * as.double(columns)^2)
  if (!is.finite(bytes) || bytes > max_memory_mb * 2^20)
    stop("Spatial-field CV exceeds max_memory_mb; use smaller assessment blocks or fewer field nodes.", call. = FALSE)
  bytes
}

.msr_field_cv_provider <- function() {
  if (!requireNamespace("glmmTMB", quietly = TRUE) ||
      !requireNamespace("TMB", quietly = TRUE) ||
      as.character(getNamespaceVersion("glmmTMB")) != "1.1.14" ||
      as.character(getNamespaceVersion("TMB")) != "1.9.25")
    stop("Spatial-field CV requires qualified glmmTMB 1.1.14 and TMB 1.9.25.", call. = FALSE)
  dll <- getLoadedDLLs()[["glmmTMB"]][["path"]]
  list(glmmTMB = "1.1.14", TMB = "1.9.25", DLL_sha256 = digest::digest(file = dll, algo = "sha256"),
    R_sha256 = digest::digest(file = file.path(find.package("glmmTMB"), "R", "glmmTMB.rdb"), algo = "sha256"))
}

.msr_field_cv_formula <- function(formula, covariance) {
  term <- call(if (covariance == "exponential") "exp" else "gau",
    call("|", call("+", as.name(".msr_field_position"), 0), as.name(".msr_field_group")))
  stats::as.formula(call("~", formula[[2L]], call("+", formula[[3L]], term)),
    env = environment(formula))
}

.msr_field_cv_positive <- function(x) is.matrix(x) && all(is.finite(x)) &&
  isTRUE(tryCatch({chol(x); TRUE}, error = function(e) FALSE))

.msr_field_cv_fit <- function(formula, data, xy, covariance, control, audit_dir,
    max_memory_mb, provider) {
  if (!identical(provider, .msr_field_cv_provider()))
    stop("The native spatial-field provider changed during CV.", call. = FALSE)
  X <- stats::model.matrix(formula, data)
  if (any(!is.finite(X))) stop("Training fixed designs must be complete and finite.", call. = FALSE)
  nodes <- unique(unname(as.matrix(xy)))
  .msr_field_cv_budget(nrow(nodes), 0, ncol(X) + 3L, max_memory_mb)
  distance <- as.matrix(stats::dist(nodes))
  if (nrow(nodes) < 3L) stop("Range estimation requires at least three distinct training locations.", call. = FALSE)
  reference_range <- stats::median(distance[distance > 0])
  response <- stats::model.response(stats::model.frame(formula, data, na.action = stats::na.fail))
  reference_sd <- stats::sd(response)
  if (!is.finite(reference_sd) || reference_sd <= 0)
    stop("Training responses need positive finite variation for field/observation SD estimation.", call. = FALSE)
  native_data <- data
  native_data$.msr_field_position <- glmmTMB::numFactor(xy)
  native_data$.msr_field_group <- factor(rep("shared", nrow(data)))
  if (nlevels(native_data$.msr_field_position) != nrow(nodes))
    stop("Native coordinate encoding merged distinct locations; retain representable coordinate precision.", call. = FALSE)
  native_formula <- .msr_field_cv_formula(formula, covariance)
  audit <- tempfile("field-cv-fit-", tmpdir = audit_dir, fileext = ".rds")
  record <- list(status = "started", provider = provider, covariance = covariance,
    training_rows = nrow(data), training_ids_sha256 = .msr_gdr_hash(rownames(data)),
    training_sha256 = .msr_gdr_hash(data),
    coordinates_sha256 = .msr_gdr_hash(xy), control = control, attempts = list())
  on.exit(.msr_gdr_atomic_audit(record, audit))
  .msr_gdr_atomic_audit(record, sub("[.]rds$", "-request.rds", audit))
  selected <- NULL; values <- numeric(2L); points <- list(); eligible <- logical(2L)
  for (i in 1:2) {
    multiplier <- c(.5, 2)[i]
    attempt <- list(start = c("reference", "cold")[i], status = "started", warnings = character(),
      theta_start = log(c(reference_sd, reference_range) * multiplier))
    started <- proc.time()[["elapsed"]]
    fit <- tryCatch(withCallingHandlers(glmmTMB::glmmTMB(native_formula,
      data = native_data, family = stats::gaussian(), REML = FALSE,
      start = list(theta = attempt$theta_start), na.action = stats::na.fail, se = FALSE,
      control = glmmTMB::glmmTMBControl(parallel = 1L, rank_check = "stop", optCtrl = control)),
      warning = function(w) {attempt$warnings <<- c(attempt$warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
      error = function(e) {attempt$error <<- conditionMessage(e); NULL})
    if (!is.null(fit)) {
      tryCatch(withCallingHandlers({
        attempt$construction_optimization_seconds <- proc.time()[["elapsed"]] - started
        # A completed numeric point survives an interruption during sdreport.
        # It is explicitly inference-ineligible and contains no model or tape.
        checkpoint <- structure(list(type = "Gaussian_field_optimizer_point",
          inference_eligible = FALSE, provider = provider, fixed = fit$fit$par,
          value = fit$fit$objective, convergence = fit$fit$convergence,
          coordinate_order = names(fit$fit$par), training_sha256 = record$training_sha256,
          coordinates_sha256 = record$coordinates_sha256, start = attempt$start,
          construction_optimization_seconds = attempt$construction_optimization_seconds),
          class = "msr_field_cv_point")
        .msr_gdr_atomic_audit(checkpoint, sub("[.]rds$", paste0("-point", i, ".rds"), audit))
        covariance_started <- proc.time()[["elapsed"]]
        report <- TMB::sdreport(fit$obj, par.fixed = fit$fit$par, getJointPrecision = TRUE)
        fit$sdr <- report
        # Keep a reproducible native call rather than references to expired
        # local native_formula/native_data/control objects. A replay requests
        # the complete native inference represented by the returned fit.
        fit$call[[1L]] <- quote(glmmTMB::glmmTMB)
        fit$call$formula <- native_formula
        fit$call$data <- native_data
        fit$call$family <- quote(stats::gaussian())
        fit$call$start <- list(theta = attempt$theta_start)
        fit$call$control <- glmmTMB::glmmTMBControl(parallel = 1L, rank_check = "stop", optCtrl = control)
        fit$call$se <- TRUE
        attempt$covariance_seconds <- proc.time()[["elapsed"]] - covariance_started
        returned_value <- fit$obj$fn(fit$fit$par)
        score <- as.vector(fit$obj$gr(fit$fit$par))
        current <- fit$obj$env$last.par
        fixed <- setdiff(seq_along(current), fit$obj$env$random)
        attempt$current_fixed_point_aligned <- isTRUE(all.equal(as.numeric(current[fixed]),
          as.numeric(fit$fit$par), tolerance = 1e-7))
        attempt$returned_value_gap <- abs(returned_value - fit$fit$objective)
        fit$fit$parfull <- current
        precision <- as.matrix(fit$obj$env$spHess(current, random = TRUE))
        attempt$convergence <- fit$fit$convergence
        attempt$value <- fit$fit$objective
        attempt$maximum_score <- max(abs(score))
        attempt$pdHess <- report$pdHess
        attempt$fixed_covariance_positive <- .msr_field_cv_positive(as.matrix(report$cov.fixed))
        attempt$joint_precision_positive <- .msr_field_cv_positive(as.matrix(report$jointPrecision))
        attempt$current_latent_precision_positive <- .msr_field_cv_positive(precision)
        attempt$fixed_coordinates <- length(fit$fit$par)
        attempt$native_K <- .msr_parameter_count(fit)
        latent <- length(fit$obj$env$random)
        attempt$coordinate_layout_verified <- latent == nrow(nodes) &&
          attempt$native_K == attempt$fixed_coordinates &&
          attempt$native_K == length(glmmTMB::fixef(fit)$cond) + 3L &&
          identical(dim(report$cov.fixed), rep(attempt$fixed_coordinates, 2L)) &&
          identical(dim(report$jointPrecision), rep(attempt$fixed_coordinates + latent, 2L)) &&
          identical(dim(precision), rep(latent, 2L))
        attempt$eligible <- identical(as.integer(fit$fit$convergence), 0L) &&
          !length(attempt$warnings) && isTRUE(attempt$pdHess) &&
          attempt$fixed_covariance_positive && attempt$joint_precision_positive &&
          attempt$current_latent_precision_positive && attempt$coordinate_layout_verified &&
          attempt$current_fixed_point_aligned && attempt$returned_value_gap < 1e-7 &&
          attempt$maximum_score < 1e-3
        values[i] <- fit$fit$objective; points[[i]] <- fit$fit$par
        eligible[i] <- attempt$eligible
        if (is.null(selected) || values[i] < selected$fit$objective) selected <- fit
        attempt$status <- "returned"
      }, warning = function(w) {attempt$warnings <<- c(attempt$warnings, conditionMessage(w)); invokeRestart("muffleWarning")}),
        error = function(e) {attempt$error <<- conditionMessage(e)})
    }
    if (!is.null(attempt$error)) attempt$status <- "failed"
    record$attempts[[i]] <- attempt
    .msr_gdr_atomic_audit(attempt, sub("[.]rds$", paste0("-start", i, ".rds"), audit))
    fit <- NULL
    gc(verbose = FALSE)
  }
  if (!all(eligible)) {
    record$status <- "failed"
    stop("Every declared spatial-field training start must pass convergence, score and full curvature gates; inspect audit_dir.", call. = FALSE)
  }
  record$NLL_gap <- abs(diff(values))
  record$coordinate_gap <- max(abs(points[[1L]] - points[[2L]]))
  if (record$NLL_gap > 1e-5 || record$coordinate_gap > 1e-3) {
    record$status <- "start_disagreement"
    stop("Declared spatial-field training starts disagree; retain the failed fold.", call. = FALSE)
  }
  record$status <- "returned"; record$audit_file <- audit
  attr(selected, "msr_field_cv_diagnostics") <- record
  selected
}

.msr_field_cv_basis <- function(formula, training, model) {
  tt <- stats::delete.response(stats::terms(formula, data = training))
  tt <- .msr_cv_prediction_terms(tt, training, qualify = TRUE)
  frame <- stats::model.frame(tt, training, na.action = stats::na.fail, drop.unused.levels = TRUE)
  X <- stats::model.matrix(tt, frame)
  native <- as.matrix(glmmTMB::getME(model, "X"))
  if (!identical(dim(X), dim(native)) || !identical(colnames(X), colnames(native)) ||
      !isTRUE(all.equal(unname(X), unname(native), tolerance = 1e-12, check.attributes = FALSE)))
    stop("Training-only field prediction recipes do not reproduce the native fixed design.", call. = FALSE)
  list(terms = tt, xlevels = lapply(frame[vapply(frame, is.factor, logical(1))], levels),
    contrasts = attr(X, "contrasts"), columns = colnames(X))
}

.msr_field_cv_predict <- function(model, formula, training, assessment, train_xy,
    test_xy, covariance, max_memory_mb) {
  basis <- .msr_field_cv_basis(formula, training, model)
  tr <- .msr_unmarked_fixed_design(basis, training)
  te <- .msr_unmarked_fixed_design(basis, assessment)
  if (any(!is.finite(te$X))) stop("Assessment fixed designs must be complete and finite.", call. = FALSE)
  nodes <- unique(unname(as.matrix(train_xy)))
  .msr_field_cv_budget(nrow(nodes), nrow(assessment), ncol(tr$X) + 3L, max_memory_mb)
  # Native numFactor labels bind rows to nodes, without a dense observation-
  # by-observation training covariance when sites have repeated surveys.
  position <- as.character(glmmTMB::numFactor(train_xy))
  group <- match(position, unique(position))
  model$obj$gr(model$fit$par)
  parameters <- model$obj$env$parList(par = model$obj$env$last.par)
  theta <- parameters$theta
  if (length(theta) != 2L || any(!is.finite(theta)))
    stop("Field CV requires the original two-parameter covariance layout.", call. = FALSE)
  field_sd <- exp(theta[1L]); range <- exp(theta[2L]); noise_sd <- stats::sigma(model)
  if (any(!is.finite(c(field_sd, range, noise_sd))) || any(c(field_sd, range, noise_sd) <= 0))
    stop("Training field range and both SDs must be positive and representable.", call. = FALSE)
  beta <- glmmTMB::fixef(model)$cond
  offset <- function(x, n) if (is.null(x)) numeric(n) else x
  train_mean <- drop(tr$X %*% beta) + offset(tr$offset, nrow(training))
  fixed_mean <- drop(te$X %*% beta) + offset(te$offset, nrow(assessment))
  y <- training[[as.character(formula[[2L]])]]
  counts <- tabulate(group, nrow(nodes))
  residual_mean <- as.numeric(rowsum(y - train_mean, group, reorder = FALSE)) / counts
  distance <- function(a, b) sqrt(outer(a[, 1L], b[, 1L], "-")^2 + outer(a[, 2L], b[, 2L], "-")^2)
  power <- if (covariance == "exponential") 1 else 2
  kernel <- function(d) field_sd^2 * exp(-(d/range)^power)
  K <- kernel(distance(nodes, nodes)); cross <- kernel(distance(test_xy, nodes))
  test <- kernel(distance(test_xy, test_xy))
  noisy <- K; diag(noisy) <- diag(noisy) + noise_sd^2/counts
  L <- t(chol(noisy))
  alpha <- backsolve(t(L), forwardsolve(L, residual_mean))
  prediction <- fixed_mean + drop(cross %*% alpha)
  loading <- forwardsolve(L, t(cross))
  field_covariance <- test - crossprod(loading)
  observation_covariance <- field_covariance
  diag(observation_covariance) <- diag(observation_covariance) + noise_sd^2
  C <- t(chol(observation_covariance))
  observed <- assessment[[as.character(formula[[2L]])]]
  whitened <- forwardsolve(C, observed - prediction)
  quadratic <- sum(whitened^2); logdet <- 2 * sum(log(diag(C)))
  n <- nrow(assessment)
  unit_prior <- t(chol(K))
  prior_loading <- forwardsolve(unit_prior, t(cross))
  new_variance <- diag(test) - colSums(prior_loading^2)
  training_variance <- diag(field_covariance) - new_variance
  allowance <- 128 * .Machine$double.eps * nrow(nodes) * pmax(1, diag(test))
  if (any(!is.finite(c(prediction, quadratic, logdet, new_variance, training_variance))) ||
      any(new_variance < -allowance) || any(training_variance < -allowance))
    stop("Conditional Gaussian field covariance is not numerically valid.", call. = FALSE)
  contract <- .msr_prediction_contract("native", "glmmTMB_Gaussian_field_CV",
    "native_component", "response", population = "fixed_at_training_ML_estimates",
    kernel = "fixed_at_training_fold_estimates",
    effects = list(response_field = "Gaussian_conditional_integration_given_training_responses"),
    conditioning = c("analysis_responses_only", "fitted_training_parameters", "frozen_training_recipes"),
    integration = list(field = "exact_Gaussian_conditioning_at_estimated_parameters", assessment = "joint_block"),
    omitted = c("coefficient_parameter_uncertainty", "range_SD_parameter_uncertainty",
      "ecological_scale_parameter_uncertainty", "selection_uncertainty"),
    coordinates = "native_fitted_parameters_and_projected_locations", diagnostics = list(passed = TRUE),
    fixed_inputs = TRUE, component = "Gaussian_joint_assessment_response")
  contract$uncertainty$observation_noise <- TRUE
  contract$uncertainty$standard_errors <- "conditional_observation_predictive_SD_not_native_mean_SE"
  contract$response_spatial_field <- list(covariance = covariance, range_estimated = TRUE,
    parameter_uncertainty_integrated = FALSE, ecological_scale_uncertainty_integrated = FALSE,
    variance_components = c("remaining_training_field", "new_location_field", "observation_noise"),
    full_grid_covariance_returned = FALSE)
  list(prediction = prediction, fitted_training_field_variance = pmax(0, training_variance),
    new_location_field_variance = pmax(0, new_variance), observation_noise_variance = noise_sd^2,
    observation_predictive_SE = sqrt(diag(observation_covariance)),
    joint = list(n = n, quadratic = quadratic, log_determinant = logdet,
      neg_log_score = .5 * (n * log(2 * pi) + logdet + quadratic),
      covariance_sha256 = .msr_gdr_hash(observation_covariance)),
    parameters = list(field_SD = field_sd, range = range, observation_SD = noise_sd,
      native_K = .msr_parameter_count(model), range_estimated = TRUE,
      estimated_range_parameters = 1L, integrated_state_penalty = 0L),
    recipe_sha256 = .msr_gdr_hash(list(basis$terms, basis$xlevels, basis$contrasts, basis$columns)),
    conditioning = "training_responses_and_ML_parameters",
    parameter_uncertainty = "not_integrated", ecological_scale_uncertainty = "not_integrated",
    observation_noise_included = TRUE, held_out_field_dependence = "joint_Gaussian_block",
    prediction_contract = contract)
}
