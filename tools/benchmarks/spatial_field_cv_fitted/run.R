source("tools/benchmarks/spatial_field_cv_quadrature/run.R")

fitted_design <- function() {
  d <- expand.grid(family = c("gaussian", "poisson", "nbinom2"),
    design = c("interpolation", "sparse_strong", "buffered_transfer"), replicate = seq_len(40L),
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  d$id <- sprintf("fit_%03d", seq_len(nrow(d)))
  d$seed <- 671200L + seq_len(nrow(d))
  d
}

fitted_simulate <- function(case) {
  set.seed(case$seed)
  if (case$design == "interpolation") {
    coords <- as.matrix(expand.grid(x = 0:5, y = 0:5))
    test <- c(2L, 5L, 10L, 13L, 18L, 23L, 28L, 33L)
    analysis <- setdiff(seq_len(nrow(coords)), test)
    visits <- 3L; field_sd <- .65; range <- 1.6; covariance <- "exponential"
  } else if (case$design == "sparse_strong") {
    coords <- as.matrix(expand.grid(x = 0:5, y = 0:3))
    test <- c(2L, 5L, 9L, 14L, 17L, 21L, 24L)
    analysis <- setdiff(seq_len(nrow(coords)), test)
    visits <- 2L; field_sd <- 1.25; range <- 3; covariance <- "exponential"
  } else {
    coords <- as.matrix(expand.grid(x = 0:7, y = 0:5))
    test <- which(coords[, 1L] >= 6)
    analysis <- which(coords[, 1L] < 4)
    visits <- 3L; field_sd <- .85; range <- 2.5; covariance <- "gaussian"
  }
  power <- if (covariance == "exponential") 1 else 2
  K <- field_sd^2 * exp(-(as.matrix(dist(coords))/range)^power)
  u <- drop(t(chol(K)) %*% rnorm(nrow(coords)))
  node <- rep(seq_len(nrow(coords)), each = visits)
  xy <- coords[node, , drop = FALSE]
  data <- data.frame(habitat = rnorm(nrow(xy)), exposure = exp(runif(nrow(xy), -.45, .45)))
  ids <- rownames(data) <- rownames(xy) <- sprintf("location_%02d_visit_%d", node, rep(seq_len(visits), length(u)))
  beta <- if (case$family == "gaussian") c(.8, .55) else c(.55, .5)
  eta <- beta[1L] + beta[2L] * data$habitat + log(data$exposure)
  link <- eta + u[node]
  dispersion <- switch(case$family, gaussian = .65, poisson = 1, nbinom2 = 3.5)
  data$y <- switch(case$family, gaussian = rnorm(length(link), link, dispersion),
    poisson = rpois(length(link), exp(link)),
    nbinom2 = rnbinom(length(link), mu = exp(link), size = dispersion))
  train_rows <- which(node %in% analysis)
  test_rows <- which(node %in% test)
  list(data = data, xy = xy, node = node, coords = coords, eta = eta, family = case$family,
    train_rows = train_rows, test_rows = test_rows, beta = beta, field_sd = field_sd,
    range = range, dispersion = dispersion, covariance = covariance, K = K,
    folds = list(split = list(analysis = ids[train_rows], assessment = ids[test_rows])))
}

generating_prediction <- function(d) {
  ai <- d$train_rows; ti <- d$test_rows
  data <- d$data
  if (d$family != "gaussian") {
    ans <- multiScaleR:::.msr_field_laplace_block(data$y[ai], data$y[ti], d$eta[ai], d$eta[ti],
      d$xy[ai, , drop = FALSE], d$xy[ti, , drop = FALSE], d$covariance, d$field_sd, d$range,
      d$family, d$dispersion, tilted_log_moments = TRUE)
    return(list(mean = ans$prediction, variance = ans$observation_predictive_SE^2,
      loss = ans$joint$neg_log_score))
  }
  take <- d$node[c(ai, ti)]
  V <- d$K[take, take, drop = FALSE]
  diag(V) <- diag(V) + d$dispersion^2
  a <- seq_along(ai); t <- length(a) + seq_along(ti)
  cross <- V[t, a, drop = FALSE]
  train <- V[a, a, drop = FALSE]
  mean <- d$eta[ti] + drop(cross %*% solve(train, data$y[ai] - d$eta[ai]))
  covariance <- V[t, t, drop = FALSE] - cross %*% solve(train, t(cross))
  lower <- t(chol(covariance))
  loss <- (length(ti) * log(2*pi) + sum(forwardsolve(lower, data$y[ti] - mean)^2))/2 +
    sum(log(diag(lower)))
  list(mean = mean, variance = diag(covariance), loss = loss)
}

append_frame <- function(x, file) {
  append <- file.exists(file)
  write.table(x, file, append = append, col.names = !append, sep = ",", row.names = FALSE)
}

audit_starts <- function(path, id) {
  paths <- list.files(path, pattern = "-start[0-9]+[.]rds$", recursive = TRUE, full.names = TRUE)
  if (!length(paths)) return(NULL)
  extract <- function(x, field) if (is.null(x[[field]])) NA else x[[field]]
  do.call(rbind, lapply(paths, function(path) {
    a <- readRDS(path)
    data.frame(id = id, attempt = extract(a, "start"), status = extract(a, "status"),
      eligible = extract(a, "eligible"), value = extract(a, "value"),
      convergence = extract(a, "convergence"), pdHess = extract(a, "pdHess"),
      maximum_score = extract(a, "maximum_score"), fixed_positive = extract(a, "fixed_covariance_positive"),
      joint_positive = extract(a, "joint_precision_positive"),
      latent_positive = extract(a, "current_latent_precision_positive"),
      current_point_aligned = extract(a, "current_fixed_point_aligned"),
      returned_value_gap = extract(a, "returned_value_gap"), layout = extract(a, "coordinate_layout_verified"),
      warnings = paste(a$warnings, collapse = " | "), error = extract(a, "error"))
  }))
}

public_prediction <- function(d, audit) {
  cv <- cross_validate_field(y ~ habitat + offset(log(exposure)), d$data, d$xy, d$folds,
    crs = 32617, audit_dir = audit, covariance = d$covariance, family = d$family)
  if (!cv$complete) {
    return(list(complete = FALSE, message = cv$errors$split$message,
      stage = cv$errors$split$stage, warnings = paste(cv$warnings$split, collapse = " | ")))
  }
  train <- cv$training$split
  pars <- train$parameters
  diagnostic <- cv$diagnostics$split
  list(complete = TRUE, message = "", stage = "", warnings = paste(cv$warnings$split, collapse = " | "),
    mean = cv$predictions$prediction, variance = cv$predictions$observation_predictive_SE^2,
    loss = cv$joint_scores$split$neg_log_score, intercept = unname(train$coefficients[1L]),
    slope = unname(train$coefficients[2L]), field_sd = pars$field_SD, range = pars$range,
    dispersion = if (d$family == "gaussian") pars$observation_SD else
      if (d$family == "poisson") 1 else pars$dispersion,
    start_value_gap = diagnostic$NLL_gap, start_parameter_gap = diagnostic$coordinate_gap)
}

fitted_run <- function(out, part, total, limit) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  scenario_file <- file.path(out, "scenarios.csv")
  prediction_file <- file.path(out, "predictions.csv")
  attempts_file <- file.path(out, "starts.csv")
  if (any(file.exists(c(scenario_file, prediction_file, attempts_file))))
    stop("Use a fresh directory to preserve original records.")
  pkgload::load_all(compile = FALSE, quiet = TRUE)
  stopifnot(as.character(packageVersion("glmmTMB")) == "1.1.14",
    as.character(packageVersion("TMB")) == "1.9.21")
  design <- fitted_design()
  cases <- design[((seq_len(nrow(design)) - 1L) %% total) + 1L == part, , drop = FALSE]
  cases <- head(cases, min(limit, nrow(cases)))
  write.csv(cases, file.path(out, "design.csv"), row.names = FALSE)
  writeLines(c(capture.output(sessionInfo()), "", "Source MD5",
    capture.output(print(tools::md5sum(c("R/cross_validate_field.R", "R/field_cv_laplace.R",
      "R/field_cv_gaussian.R", "tools/benchmarks/spatial_field_cv_fitted/run.R"))))),
    file.path(out, "session.txt"))
  for (i in seq_len(nrow(cases))) {
    started <- proc.time()[["elapsed"]]
    case <- cases[i, , drop = FALSE]
    d <- fitted_simulate(case)
    audit <- file.path(out, "audit", case$id)
    fitted_call <- capture(function() public_prediction(d, audit))
    fit <- fitted_call$value
    oracle_call <- capture(function() generating_prediction(d))
    oracle <- oracle_call$value
    get <- function(x, nm, default = NA_real_) if (is.null(x[[nm]])) default else x[[nm]]
    y <- d$data$y[d$test_rows]
    pred <- data.frame(id = case$id, row = seq_along(y), site = d$node[d$test_rows], observed = y,
      fitted_mean = get(fit, "mean"), fitted_variance = get(fit, "variance"),
      generating_mean = get(oracle, "mean"), generating_variance = get(oracle, "variance"))
    if (case$family == "gaussian") {
      pred$fitted_pit <- pnorm(y, pred$fitted_mean, sqrt(pred$fitted_variance))
      pred$generating_pit <- pnorm(y, pred$generating_mean, sqrt(pred$generating_variance))
    } else {
      pred$fitted_pit <- pred$generating_pit <- NA_real_
    }
    starts <- audit_starts(audit, case$id)
    if (!is.null(starts)) append_frame(starts, attempts_file)
    append_frame(pred, prediction_file)
    append_frame(cbind(case, data.frame(covariance = d$covariance,
      analysis_locations = length(unique(d$node[d$train_rows])),
      assessment_locations = length(unique(d$node[d$test_rows])),
      n_analysis = length(d$train_rows), n_assessment = length(y),
      fit_status = if (isTRUE(fit$complete)) "scored" else "failed",
      fit_stage = get(fit, "stage", "call"), fit_message = get(fit, "message", fitted_call$message),
      fit_warnings = paste(fitted_call$warnings, get(fit, "warnings", "")),
      generating_status = if (is.null(oracle)) "failed" else "scored",
      generating_message = oracle_call$message, generating_warnings = oracle_call$warnings,
      fitted_loss = get(fit, "loss"), generating_loss = get(oracle, "loss"),
      fitted_intercept = get(fit, "intercept"), fitted_slope = get(fit, "slope"),
      fitted_field_sd = get(fit, "field_sd"), fitted_range = get(fit, "range"),
      fitted_dispersion = get(fit, "dispersion"), true_intercept = d$beta[1L],
      true_slope = d$beta[2L], true_field_sd = d$field_sd, true_range = d$range,
      true_dispersion = d$dispersion, start_value_gap = get(fit, "start_value_gap"),
      start_parameter_gap = get(fit, "start_parameter_gap"),
      elapsed_seconds = proc.time()[["elapsed"]] - started)), scenario_file)
    if (i %% 10L == 0L || i == nrow(cases))
      cat("part", part, i, "of", nrow(cases), "scenarios recorded\n")
  }
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("Supply a fresh output directory.")
  fitted_run(args[1L], if (length(args) > 1L) as.integer(args[2L]) else 1L,
    if (length(args) > 2L) as.integer(args[3L]) else 1L,
    if (length(args) > 3L) as.integer(args[4L]) else Inf)
}
