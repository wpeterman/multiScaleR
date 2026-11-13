# Run from the package root with Rscript tools/benchmarks/spatial_field_cv_quadrature/run.R.
# Arguments optionally set a fresh output directory and a smoke-test row limit.

# Standard-normal quadrature from the probabilists' Hermite Jacobi matrix.
normal_rule <- local({
  cache <- new.env(parent = emptyenv())
  function(n) {
    key <- as.character(n)
    if (!exists(key, cache, inherits = FALSE)) {
      j <- matrix(0, n, n)
      ij <- cbind(seq_len(n - 1L), seq.int(2L, n))
      j[ij] <- sqrt(seq_len(n - 1L))
      e <- eigen(j + t(j), symmetric = TRUE)
      cache[[key]] <- list(nodes = e$values, weights = e$vectors[1L, ]^2)
    }
    cache[[key]]
  }
})

tensor_rule <- local({
  cache <- new.env(parent = emptyenv())
  function(n, d) {
    key <- paste(n, d)
    if (!exists(key, cache, inherits = FALSE)) {
      rule <- normal_rule(n)
      ix <- as.matrix(expand.grid(rep(list(seq_len(n)), d),
        KEEP.OUT.ATTRS = FALSE))
      w <- rowSums(matrix(log(rule$weights[ix]), ncol = d))
      cache[[key]] <- list(nodes = matrix(rule$nodes[ix], ncol = d), log_weight = w)
    }
    cache[[key]]
  }
})

log_sum <- function(x) {
  largest <- max(x)
  if (!is.finite(largest)) return(largest)
  largest + log(sum(exp(x - largest)))
}

response_logdens <- function(y, eta, family, dispersion) {
  switch(family, poisson = dpois(y, exp(eta), log = TRUE),
    binomial = dbinom(y, 1, plogis(eta), log = TRUE),
    nbinom2 = dnbinom(y, size = dispersion, mu = exp(eta), log = TRUE),
    Gamma = dgamma(y, shape = dispersion, scale = exp(eta)/dispersion, log = TRUE),
    gaussian = dnorm(y, eta, dispersion, log = TRUE))
}

likelihood <- function(z, y, eta, groups, lower, family, dispersion) {
  field <- z %*% t(lower)
  total <- numeric(nrow(z))
  for (i in seq_along(y))
    total <- total + response_logdens(y[i], eta[i] + field[, groups[i]],
      family, dispersion)
  total
}

prepare_reference <- function(y, eta, groups, lower, family, dispersion) {
  loglik <- function(z) likelihood(z, y, eta, groups, lower, family, dispersion)
  objective <- function(z) sum(z^2)/2 - loglik(matrix(z, nrow = 1L))
  fit <- optim(numeric(ncol(lower)), objective, method = "BFGS",
    hessian = TRUE, control = list(reltol = 1e-12, maxit = 500))
  if (fit$convergence != 0L || !is.finite(fit$value))
    stop("Numerical quadrature proposal did not converge.")
  list(loglik = loglik, center = fit$par, scale = chol(chol2inv(chol(fit$hessian))))
}

reference_integral <- function(reference, n, spread) {
  rule <- tensor_rule(n, length(reference$center))
  scale <- spread * reference$scale
  draws <- sweep(rule$nodes %*% scale, 2L, reference$center, "+")
  # The determinant and both Gaussian normalizing constants are included
  # through the normalized prior/proposal ratio in whitened coordinates.
  log_weight <- rule$log_weight + reference$loglik(draws) +
    (rowSums(rule$nodes^2) - rowSums(draws^2))/2 + sum(log(diag(scale)))
  value <- log_sum(log_weight)
  if (!is.finite(value)) stop("Quadrature returned a non-finite marginal likelihood.")
  list(value = value, draws = draws, weight = exp(log_weight - value))
}

response_moments <- function(analysis, fixed, projection, variance, family, n) {
  shifts <- sweep(analysis$draws %*% t(projection), 2L, fixed, "+")
  if (family == "gaussian") return(drop(crossprod(analysis$weight, shifts)))
  if (family != "binomial") {
    conditional <- exp(sweep(shifts, 2L, variance/2, "+"))
    return(drop(crossprod(analysis$weight, conditional)))
  }
  rule <- normal_rule(n)
  vapply(seq_along(fixed), function(i) {
    probabilities <- plogis(outer(shifts[, i], sqrt(variance[i]) * rule$nodes, "+"))
    sum(analysis$weight * drop(probabilities %*% rule$weights))
  }, 0)
}

reference_score <- function(d) {
  training <- seq_along(d$train_y)
  assessment <- length(training) + seq_along(d$test_y)
  nodes <- seq_len(d$n_train)
  new_nodes <- setdiff(seq_len(ncol(d$lower)), nodes)
  ref_train <- prepare_reference(d$train_y, d$train_eta, d$groups[training],
    d$lower[nodes, nodes, drop = FALSE], d$family, d$dispersion)
  ref_joint <- prepare_reference(c(d$train_y, d$test_y), c(d$train_eta, d$test_eta),
    d$groups, d$lower, d$family, d$dispersion)
  projection <- d$lower[d$groups[assessment], nodes, drop = FALSE]
  innovation <- rowSums(d$lower[d$groups[assessment], new_nodes, drop = FALSE]^2)
  at_order <- function(order, spread) {
    analysis <- reference_integral(ref_train, order, spread)
    joint <- reference_integral(ref_joint, order, spread)
    list(analysis = analysis$value, combined = joint$value,
      loss = analysis$value - joint$value,
      mean = response_moments(analysis, d$test_eta, projection, innovation,
        d$family, max(48L, order)))
  }
  gap <- function(a, b) c(
    log = max(abs(c(a$analysis - b$analysis, a$combined - b$combined))),
    mean = max(abs(a$mean - b$mean)/pmax(1, abs(b$mean))))
  previous <- at_order(14L, 1.15)
  alternate <- c(log = NA_real_, mean = NA_real_)
  stable <- FALSE
  tolerance <- 5e-6
  for (order in c(22L, 32L, 48L, 64L, 80L)) {
    current <- at_order(order, 1.15)
    refinement <- gap(previous, current)
    if (all(refinement < tolerance) || order == 80L) {
      alternate <- gap(current, at_order(order, 1.4))
      stable <- all(refinement < tolerance) && all(alternate < tolerance)
      if (stable) break
    }
    previous <- current
  }
  c(current, list(converged = stable, order = order,
    log_refinement_gap = unname(refinement["log"]),
    mean_refinement_gap = unname(refinement["mean"]),
    log_proposal_gap = unname(alternate["log"]),
    mean_proposal_gap = unname(alternate["mean"])))
}

simulate_case <- function(case) {
  set.seed(case$seed)
  layouts <- list("1+1" = c(0, 1), "2+1" = c(0, 2, 1), "1+2" = c(0, 1, 2))
  position <- layouts[[case$layout]]
  n_train <- if (case$layout == "2+1") 2L else 1L
  coordinates <- cbind(position, 0)
  power <- if (case$covariance == "exponential") 1 else 2
  correlation <- exp(log(case$rho) * as.matrix(dist(coordinates))^power)
  lower <- case$field_sd * t(chol(correlation))
  groups <- rep(seq_along(position), each = case$visits)
  n <- length(groups)
  eta <- if (case$family == "binomial") -.4 + rnorm(n, sd = .65) else
    .45 + rnorm(n, sd = .45)
  field <- drop(lower %*% rnorm(length(position)))
  link <- eta + field[groups]
  dispersion <- switch(case$family, nbinom2 = 3, Gamma = 4, gaussian = .65, 1)
  y <- switch(case$family, poisson = rpois(n, exp(link)),
    binomial = rbinom(n, 1L, plogis(link)), nbinom2 = rnbinom(n, size = dispersion, mu = exp(link)),
    Gamma = rgamma(n, shape = dispersion, scale = exp(link)/dispersion),
    gaussian = rnorm(n, link, dispersion))
  training <- which(groups <= n_train)
  test <- which(groups > n_train)
  list(train_y = y[training], test_y = y[test], train_eta = eta[training],
    test_eta = eta[test], train_xy = coordinates[groups[training], , drop = FALSE],
    test_xy = coordinates[groups[test], , drop = FALSE], groups = groups,
    lower = lower, n_train = n_train, family = case$family, dispersion = dispersion,
    range = (-1/log(case$rho))^(1/power))
}

direct_zero <- function(d) {
  analysis <- sum(response_logdens(d$train_y, d$train_eta, d$family, d$dispersion))
  loss <- -sum(response_logdens(d$test_y, d$test_eta, d$family, d$dispersion))
  list(analysis = analysis, combined = analysis - loss, loss = loss,
    mean = if (d$family == "binomial") plogis(d$test_eta) else exp(d$test_eta),
    converged = TRUE)
}

dense_gaussian <- function(d) {
  field <- tcrossprod(d$lower)[d$groups, d$groups, drop = FALSE]
  diag(field) <- diag(field) + d$dispersion^2
  a <- seq_along(d$train_y)
  b <- length(a) + seq_along(d$test_y)
  training <- field[a, a, drop = FALSE]
  cross <- field[b, a, drop = FALSE]
  mu <- d$test_eta + drop(cross %*% solve(training, d$train_y - d$train_eta))
  variance <- field[b, b, drop = FALSE] - cross %*% solve(training, t(cross))
  L <- t(chol(variance))
  loss <- (length(b) * log(2 * pi) + sum(forwardsolve(L, d$test_y - mu)^2))/2 +
    sum(log(diag(L)))
  list(loss = loss, mean = mu)
}

capture <- function(fun) {
  messages <- character()
  value <- tryCatch(withCallingHandlers(fun(), warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = if (inherits(value, "error")) NULL else value,
    message = if (inherits(value, "error")) conditionMessage(value) else "",
    warnings = paste(unique(messages), collapse = " | "))
}

run_study <- function(out, limit) {
  score_file <- file.path(out, "scenarios.csv")
  prediction_file <- file.path(out, "predictions.csv")
  if (any(file.exists(c(score_file, prediction_file))))
    stop("Choose a fresh output directory for the run to preserve previous scenario records.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  pkgload::load_all(compile = FALSE, quiet = TRUE)
  stopifnot(requireNamespace("glmmTMB", quietly = TRUE))
  writeLines(c(capture.output(sessionInfo()), "",
    "MD5 of package sources used for the calculation",
    capture.output(print(tools::md5sum(file.path("R", c("field_cv_laplace.R",
      "field_cv_gaussian.R", "scale_response_likelihood.R")))))),
    file.path(out, "session.txt"))

  cases <- expand.grid(family = c("poisson", "binomial", "nbinom2", "Gamma"),
    field_sd = c(.25, 1, 2), rho = c(.2, .85), visits = c(1L, 8L),
    layout = c("1+1", "2+1", "1+2"), replicate = seq_len(8L),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  cases$covariance <- ifelse(cases$replicate <= 4L, "exponential", "gaussian")
  zero <- subset(cases, field_sd == .25 & rho == .2)
  zero$field_sd <- 0
  gaussian <- subset(cases, family == "poisson")
  gaussian$family <- "gaussian"
  cases <- rbind(transform(cases, cohort = "Laplace"),
    transform(zero, cohort = "zero_field"), transform(gaussian, cohort = "exact_gaussian"))
  cases$id <- sprintf("field_%04d", seq_len(nrow(cases)))
  cases$seed <- 742607L + seq_len(nrow(cases))
  cases <- head(cases, min(limit, nrow(cases)))
  write.csv(cases, file.path(out, "design.csv"), row.names = FALSE)
  cat("Executing", nrow(cases), "scenarios\n")
  pred_header <- TRUE
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, , drop = FALSE]
    d <- simulate_case(case)
    reference <- capture(function() {
      if (case$cohort == "zero_field") direct_zero(d) else reference_score(d)
    })
    evaluated <- capture(function() {
      if (case$cohort == "exact_gaussian") return(dense_gaussian(d))
      fit <- multiScaleR:::.msr_field_laplace_block(d$train_y, d$test_y, d$train_eta, d$test_eta,
        d$train_xy, d$test_xy, case$covariance, case$field_sd, d$range,
        d$family, d$dispersion)
      list(loss = fit$joint$neg_log_score, mean = fit$prediction,
        analysis = fit$joint$analysis_log_density, combined = fit$joint$combined_log_density,
        analysis_score = fit$joint$analysis_latent_score, combined_score = fit$joint$combined_latent_score)
    })
    ref <- reference$value
    fit <- evaluated$value
    get <- function(x, nm) if (is.null(x[[nm]])) NA_real_ else x[[nm]]
    ok <- isTRUE(ref$converged) && !is.null(fit)
    mean_error <- if (ok) fit$mean - ref$mean else rep(NA_real_, length(d$test_y))
    scaled_error <- if (d$family %in% c("binomial", "gaussian")) abs(mean_error) else
      abs(mean_error)/ref$mean
    if (!ok) scaled_error <- mean_error
    row <- cbind(case, data.frame(n_analysis = length(d$train_y), n_assessment = length(d$test_y),
      reference_status = if (is.null(ref)) "failed" else if (isTRUE(ref$converged)) "converged" else "unresolved",
      evaluation_status = if (is.null(fit)) "failed" else "returned",
      reference_message = reference$message, reference_warnings = reference$warnings,
      evaluation_message = evaluated$message, evaluation_warnings = evaluated$warnings,
      quadrature_order = get(ref, "order"), log_refinement_gap = get(ref, "log_refinement_gap"),
      mean_refinement_gap = get(ref, "mean_refinement_gap"), log_proposal_gap = get(ref, "log_proposal_gap"),
      mean_proposal_gap = get(ref, "mean_proposal_gap"),
      reference_loss = get(ref, "loss"), evaluated_loss = get(fit, "loss"),
      reference_analysis = get(ref, "analysis"), evaluated_analysis = get(fit, "analysis"),
      reference_combined = get(ref, "combined"), evaluated_combined = get(fit, "combined"),
      analysis_mode_score = get(fit, "analysis_score"), combined_mode_score = get(fit, "combined_score"),
      loss_error = if (ok) fit$loss - ref$loss else NA_real_,
      max_mean_error = max(abs(mean_error)), max_scaled_mean_error = max(scaled_error)))
    write.table(row, score_file, row.names = FALSE, sep = ",", col.names = i == 1L,
      append = i > 1L)
    means <- function(x) if (is.null(x)) rep(NA_real_, length(d$test_y)) else x$mean
    pred <- data.frame(id = case$id, assessment_row = seq_along(d$test_y), observed = d$test_y,
      reference_mean = means(ref), evaluated_mean = means(fit),
      mean_error = mean_error, scaled_abs_mean_error = scaled_error)
    write.table(pred, prediction_file, row.names = FALSE, sep = ",", col.names = pred_header,
      append = !pred_header)
    pred_header <- FALSE
    if (i %% 48L == 0L || i == nrow(cases)) cat(i, "of", nrow(cases), "recorded\n")
  }
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- if (length(args)) args[1L] else
    "tools/benchmarks/spatial_field_cv_quadrature/results"
  limit <- if (length(args) > 1L) as.integer(args[2L]) else Inf
  run_study(out, limit)
}
