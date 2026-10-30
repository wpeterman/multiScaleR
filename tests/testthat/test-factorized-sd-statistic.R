response_sd_statistic_fixture <- function() {
  f <- bayesian_scale_fixture(n_site = 6L, n_year = 5L, seed = 713403)
  prepare <- function(data) multiScaleR:::.msr_marginal_prepare(
    count ~ habitat + temp + (1 | year) + offset(log(effort)), f$inputs, data,
    list(habitat = ~1), "gaussian", "key", 200)
  spec <- prepare(f$data)
  p <- spec$parameters; p$beta <- c(.5, .25, 2); p$scale_beta[] <- log(.2)
  p$response_z[] <- c(-1, -.3, 0, .4, .9); p$log_response_sd[] <- log(.35)
  p$log_dispersion[] <- log(.2)
  eta <- multiScaleR:::.msr_bayes_components(spec, p)$eta
  n <- nrow(f$data)
  residual <- qnorm((((seq_len(n) * 7L - 1L) %% n) + .5) / n)
  f$data$count <- eta + .2 * residual
  list(spec = prepare(f$data), parameters = p)
}

response_statistic_reference <- function(spec, p, sd) {
  eta <- multiScaleR:::.msr_bayes_components(spec, p, include_response = FALSE)$eta
  Z <- model.matrix(~0 + factor(spec$response$random[[1L]]$index))
  covariance <- diag(exp(2 * p$log_dispersion), length(spec$y)) + sd^2 * tcrossprod(Z)
  residual <- spec$y - eta
  length(residual) / 2 * log(2 * pi) + sum(log(diag(chol(covariance)))) +
    drop(crossprod(residual, solve(covariance, residual))) / 2
}

test_that("response SD statistics use a qualified same-target unrestricted candidate", {
  skip_if_not_installed("RTMB")
  f <- response_sd_statistic_fixture()
  before <- serialize(f$spec, NULL); rng <- .Random.seed
  result <- multiScaleR:::.msr_factorized_response_sd_statistic(f$spec, c(0, .2, .5),
    initial_parameters = list(f$parameters))
  observer <- getOption("multiScaleR.statistic_test_observer")
  if (is.function(observer)) observer(result)
  expect_s3_class(result, "msr_response_sd_statistic")
  expect_identical(serialize(f$spec, NULL), before)
  expect_identical(.Random.seed, rng)
  expect_identical(result$curve$sd, c(0, .2, .5))
  expect_identical(result$target, "unpenalized full marginal likelihood")
  expect_false(result$global_fit_eligible)
  expect_false(result$uncertainty$available)
  expect_false(any(c("cutoff", "confidence_level", "interval", "p_value") %in% names(result)))
  expect_true(result$unrestricted$diagnostics$passed)
  expect_true(result$reference$qualified)
  expect_true(all(result$curve$available))
  expect_true(result$unrestricted$diagnostics$boundary_search_complete)
  expect_true(result$unrestricted$diagnostics$limiting_audit_complete)
  expect_identical(result$unrestricted$inference_metadata,
    multiScaleR:::.msr_scale_inference_metadata("marginal_likelihood", "factorized"))
  expect_length(result$profile$points, 3L)
  expect_length(result$warm_starts$from_profile, 3L)
  expect_identical(result$warm_starts$supplied, list(f$parameters))
  full <- result$unrestricted
  M <- response_statistic_reference(f$spec, full$parameters,
    if (full$zero_sd[["response"]]) 0 else exp(full$parameters$log_response_sd))
  expect_equal(result$reference$ordinary_objective, M, tolerance = 1e-8)
  expect_true(result$audits$unrestricted$passed)
  expect_false(result$reference$lower_feasible_profile_witness)
  expect_gte(nrow(result$reference$witnesses), 12L)
  expect_true(all(c("profile", "broad", "nearest") %in% result$reference$witnesses$kind))
  for (k in seq_along(result$profile$points)) {
    point <- result$profile$points[[k]]
    P <- response_statistic_reference(f$spec, point$parameters, point$sd)
    expect_length(point$attempts, 4L)
    expect_equal(result$curve$profile_objective[k], P, tolerance = 1e-8)
    expect_equal(result$curve$raw_statistic[k], 2 * (P - M), tolerance = 1e-7)
    expect_equal(result$curve$statistic[k], result$curve$raw_statistic[k])
    expect_true(result$audits$points[[k]]$passed)
    expect_gte(P, M)
    expect_false(point$global_fit_eligible)
    expect_identical(point$limits$fixed_sd, point$sd)
    for (limit in point$limits$references) expect_length(limit$attempts, 3L)
    audit <- result$audits$points[[k]]
    if (point$sd > 0) {
      expect_identical(audit$parameter_blocks[audit$fixed_indices], "log_response_sd")
      expect_gt(abs(audit$full_gradient[audit$fixed_indices]), 1e-4)
      expect_true(result$warm_starts$mapping[[k]][[1L]]$present_after_original_deduplication)
    } else {
      expect_length(audit$fixed_indices, 0L)
      expect_false("log_response_sd" %in% audit$parameter_blocks)
      expect_false(result$warm_starts$mapping[[k]][[1L]]$present_after_original_deduplication)
    }
    expect_true(result$warm_starts$mapping[[k]][[2L]]$present_after_original_deduplication)
  }
  expect_true(all(c(FALSE, TRUE) %in% vapply(full$attempts,
    function(a) isTRUE(a$zero_sd[["response"]]), TRUE)))
  faces <- vapply(full$attempts, function(a) isTRUE(a$zero_sd[["response"]]), TRUE)
  expect_gte(sum(faces), 3L)
  expect_gte(sum(!faces), 3L)
  expect_identical(full$diagnostics$attempted_starts, length(full$attempts))
  zero <- result$profile$points[[1L]]
  zero_full <- multiScaleR:::.msr_response_sd_statistic_audit(f$spec, zero$parameters,
    TRUE, zero$objective, fixed = FALSE)
  expect_true(zero_full$passed)
  expect_equal(zero_full$ordinary_objective, zero$objective, tolerance = 1e-8)
  expect_length(zero_full$fixed_indices, 0L)
  expect_false("log_response_sd" %in% zero_full$parameter_blocks)
  expect_identical(zero_full$stationarity, "all retained population coordinates")
  # Cold failures and retained-vector discrepancies remain explicit, independently
  # of whether the fitted candidate originally passed its own numerical gates.
  bad <- multiScaleR:::.msr_response_sd_statistic_audit(f$spec, full$parameters,
    isTRUE(full$zero_sd[["response"]]), M + .01)
  expect_false(bad$passed)
  expect_gt(bad$value_error, 1e-6)
  local_mocked_bindings(.msr_factorized_adfun = function(...) stop("retained cold construction failure"),
    .package = "multiScaleR")
  failed <- multiScaleR:::.msr_response_sd_statistic_audit(f$spec, full$parameters,
    isTRUE(full$zero_sd[["response"]]), M)
  expect_false(failed$passed)
  expect_match(failed$error, "retained cold")
})

test_that("negative discrepancies, lower witnesses and original rejections are retained", {
  f <- response_sd_statistic_fixture()
  p <- f$parameters; p$scale_z <- p$response_z <- numeric()
  profile <- list(target = "unpenalized full marginal likelihood", parameter = "response_group_sd",
    units = "response linear predictor", curve = data.frame(sd = c(0, .2)),
    points = lapply(c(0, .2), function(sd) list(sd = sd, parameters = p, objective = 10,
      available = TRUE, reason = "original conditional reason", attempts = list(
        list(audited_value = TRUE, objective = 10, parameters = p)),
      limits = list(references = list()))))
  fit <- list(parameters = p, zero_sd = c(scale = FALSE, response = FALSE),
    log_likelihood = -9, attempts = list(), diagnostics = list(passed = TRUE, reason = "original fit reason"))
  full_value <- 9; profile_value <- 10; full_pass <- TRUE; point_pass <- TRUE
  local_mocked_bindings(.msr_factorized_response_sd_profile = function(...) profile,
    .msr_factorized_fit = function(...) fit,
    .msr_response_sd_statistic_audit = function(spec, parameters, zero, objective, fixed = FALSE)
      list(passed = if (fixed) point_pass else full_pass,
        ordinary_objective = if (fixed) profile_value else full_value), .package = "multiScaleR")
  call <- function() multiScaleR:::.msr_factorized_response_sd_statistic(f$spec, c(0, .2))
  good <- call()
  expect_true(all(good$curve$available))
  expect_equal(good$curve$statistic, c(2, 2))
  # A selected exact-zero face uses equality without an epsilon SD or clamping.
  full_value <- 10; fit$zero_sd[["response"]] <- TRUE
  fit$parameters$log_response_sd <- numeric()
  equal <- call()
  expect_true(equal$reference$qualified)
  expect_true(all(equal$curve$available))
  expect_identical(equal$curve$raw_statistic, c(0, 0))
  expect_identical(equal$curve$statistic, c(0, 0))
  expect_true(equal$unrestricted$zero_sd[["response"]])
  full_value <- 9; fit$zero_sd[["response"]] <- FALSE; fit$parameters <- p
  # Even a negative discrepancy smaller than the value-audit tolerance is retained.
  profile_value <- 9 - 1e-10
  negative <- call()
  expect_false(negative$reference$qualified)
  expect_true(negative$reference$lower_feasible_profile_witness)
  expect_true(all(negative$curve$raw_statistic < 0))
  expect_true(all(is.na(negative$curve$statistic)))
  expect_false(any(negative$curve$available))
  expect_match(negative$curve$reason[1L], "Negative")
  profile_value <- 10
  profile$points[[2L]]$attempts[[2L]] <- list(audited_value = TRUE,
    eligible_nuisance = FALSE, objective = 8.5, error = "retained rejected witness")
  lower <- call()
  expect_false(lower$reference$qualified)
  expect_true(lower$reference$lower_feasible_profile_witness)
  expect_equal(sum(lower$reference$witnesses$kind == "profile"), 3L)
  expect_false(any(lower$curve$available))
  expect_match(lower$profile$points[[2L]]$attempts[[2L]]$error, "retained rejected")
  profile$points[[2L]]$attempts <- profile$points[[2L]]$attempts[1L]
  fit$diagnostics$passed <- FALSE
  rejected <- call()
  expect_identical(rejected$unrestricted$diagnostics, fit$diagnostics)
  expect_false(any(rejected$curve$available))
  expect_match(rejected$curve$reason[1L], "original fit reason")
  fit$diagnostics$passed <- TRUE; full_pass <- FALSE
  expect_match(call()$curve$reason[1L], "unrestricted physical-vector")
  full_pass <- TRUE; point_pass <- FALSE
  expect_match(call()$curve$reason[1L], "profile physical-vector")
  point_pass <- TRUE; profile$points[[2L]]$available <- FALSE
  partial <- call()
  expect_identical(partial$curve$available, c(TRUE, FALSE))
  expect_match(partial$curve$reason[2L], "original conditional reason")
  local_mocked_bindings(.msr_factorized_fit = function(spec, control, initial_parameters, checkpoint) {
    checkpoint(list(list(error = "earlier completed attempt")))
    stop("terminal fitting failure")
  }, .package = "multiScaleR")
  failed <- call()
  expect_false(failed$reference$qualified)
  expect_match(failed$unrestricted$error, "terminal fitting")
  expect_match(failed$unrestricted$attempts[[1L]]$error, "earlier completed")
  expect_false(any(failed$curve$available))
})

test_that("statistic wrapper delegates unchanged scope and control validation", {
  f <- response_sd_statistic_fixture()
  statistic <- multiScaleR:::.msr_factorized_response_sd_statistic
  expect_error(statistic(f$spec, .2), "sd")
  expect_error(statistic(f$spec, c(0, .2), list(cutoff = 3.84)), "control")
  posterior <- f$spec; posterior$priors <- list(scale_median = 200)
  expect_error(statistic(posterior, c(0, .2)), "unpenalized")
  constant <- f$spec; constant$offset[] <- 0
  constant$y <- as.numeric(factor(constant$response$random[[1L]]$index))
  rejected <- statistic(constant, c(0, .2))
  expect_false(rejected$reference$qualified)
  expect_false(any(rejected$curve$available))
  expect_length(rejected$unrestricted$attempts, 0L)
  expect_match(paste(rejected$unrestricted$diagnostics$reason, collapse = " "), "unbounded")
})
