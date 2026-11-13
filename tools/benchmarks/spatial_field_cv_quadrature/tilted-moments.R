# Follow-up comparison of tilted first and second moments at frozen parameters.
source("tools/benchmarks/spatial_field_cv_quadrature/run.R")

# Exercise the package's tilted moment solver at the original study parameters.
tilted_log_moments <- function(d, analysis) {
  ti <- seq_along(d$train_y)
  ai <- length(ti) + seq_along(d$test_y)
  multiScaleR:::.msr_field_tilted_log_moments(d$train_y, d$train_eta,
    d$groups[ti], d$lower, seq_len(d$n_train), d$groups[ai], d$test_eta,
    d$family, d$dispersion, analysis)
}

independent_log_moments <- function(d, initial_order) {
  ti <- seq_along(d$train_y)
  ai <- length(ti) + seq_along(d$test_y)
  training_nodes <- seq_len(d$n_train)
  assessment_nodes <- setdiff(seq_len(ncol(d$lower)), training_nodes)
  projection <- d$lower[d$groups[ai], training_nodes, drop = FALSE]
  innovation <- rowSums(d$lower[d$groups[ai], assessment_nodes, drop = FALSE]^2)
  prior <- prepare_reference(d$train_y, d$train_eta, d$groups[ti],
    d$lower[training_nodes, training_nodes, drop = FALSE], d$family, d$dispersion)
  at_order <- function(order, spread) {
    ref <- reference_integral(prior, order, spread)
    shifts <- sweep(ref$draws %*% t(projection), 2L, d$test_eta, "+")
    log_weight <- log(ref$weight)
    matrix(vapply(1:2, function(k) {
      values <- sweep(k * shifts, 2L, k^2 * innovation/2, "+")
      apply(values, 2L, function(value) log_sum(log_weight + value))
    }, numeric(length(ai))), ncol = 2L)
  }
  tolerance <- 5e-6
  previous <- at_order(max(22L, min(48L, initial_order)), 1.15)
  orders <- c(32L, 48L, 64L, 80L, 112L)
  orders <- orders[orders > max(22L, min(48L, initial_order))]
  converged <- FALSE
  alternate <- NA_real_
  for (order in orders) {
    current <- at_order(order, 1.15)
    gap <- max(abs(current - previous))
    if (gap < tolerance || order == max(orders)) {
      alternate <- max(abs(current - at_order(order, 1.4)))
      converged <- gap < tolerance && alternate < tolerance
      if (converged) break
    }
    previous <- current
  }
  list(log_moment = current, converged = converged, order = order, gap = gap, alternate = alternate)
}

predictive_moments <- function(log_moment, family, dispersion) {
  mean <- exp(log_moment[, 1L])
  second <- exp(log_moment[, 2L])
  log_gap <- 2 * log_moment[, 1L] - log_moment[, 2L]
  rounding <- 128 * .Machine$double.eps * pmax(1, abs(log_moment[, 2L]))
  if (any(log_gap > rounding))
    stop("The integrated second moment is smaller than the square of the mean.")
  noise <- switch(family, poisson = mean, nbinom2 = mean + second/dispersion,
    Gamma = second/dispersion)
  list(mean = mean, second = second,
    variance = noise + second * (-expm1(pmin(0, log_gap))))
}

tilted_approximation <- function(d, analysis) {
  log_moment <- tilted_log_moments(d, analysis)
  response <- multiScaleR:::.msr_field_log_response_moments(d$family,
    log_moment[, 1L], log_moment[, 2L], d$dispersion)
  list(mean = response$mean, variance = response$variance, second = exp(log_moment[, 2L]))
}

gaussian_approximation <- function(d, field_sd, covariance) {
  fit <- multiScaleR:::.msr_field_laplace_block(d$train_y, d$test_y, d$train_eta,
    d$test_eta, d$train_xy, d$test_xy, covariance, field_sd, d$range, d$family, d$dispersion)
  list(mean = fit$prediction, variance = fit$observation_predictive_SE^2,
    second = fit$prediction^2 * exp(fit$fitted_training_field_variance + fit$new_location_field_variance))
}

analysis_laplace <- function(d) {
  nodes <- seq_len(d$n_train)
  multiScaleR:::.msr_field_laplace(d$train_y, d$train_eta, d$groups[seq_along(d$train_y)],
    d$lower[nodes, nodes, drop = FALSE], d$family, d$dispersion)
}

fit_vignette_fixture <- function() {
  set.seed(412)
  locations <- as.matrix(expand.grid(x = seq(0, 70, by = 10), y = seq(0, 60, by = 10)))
  site <- rep(seq_len(nrow(locations)), each = 3)
  xy <- locations[site, , drop = FALSE]
  K <- 1.3^2 * exp(-as.matrix(dist(locations)) / 25)
  field <- drop(t(chol(K)) %*% rnorm(nrow(locations)))
  data <- data.frame(habitat = rnorm(nrow(xy)), exposure = exp(runif(nrow(xy), -.3, .3)))
  eta <- .4 + .5 * data$habitat + log(data$exposure) + field[site]
  data$count <- rpois(nrow(data), exp(eta))
  ids <- rownames(data) <- rownames(xy) <- sprintf("visit_%03d", seq_len(nrow(xy)))
  ai <- which(site %in% c(1:3, 33:35))
  ti <- setdiff(seq_along(site), ai)
  folds <- list(block = list(analysis = ids[ti], assessment = ids[ai]))
  fit <- cross_validate_field(count ~ habitat + offset(log(exposure)), data, xy, folds,
    crs = 32617, family = "poisson", audit_dir = tempfile("count-field-cv-cost-"))
  stopifnot(fit$complete)
  parameters <- fit$training$block$parameters
  beta <- fit$training$block$coefficients
  fixed <- drop(model.matrix(~ habitat, data) %*% beta) + log(data$exposure)
  order <- unique(site[c(ti, ai)])
  lower <- parameters$field_SD *
    t(chol(exp(-as.matrix(dist(locations[order, , drop = FALSE]))/parameters$range)))
  list(d = list(train_y = data$count[ti], test_y = data$count[ai], train_eta = fixed[ti],
    test_eta = fixed[ai], train_xy = xy[ti, , drop = FALSE], test_xy = xy[ai, , drop = FALSE],
    groups = match(site[c(ti, ai)], order), lower = lower,
    n_train = length(unique(site[ti])), family = "poisson", dispersion = 1,
    range = parameters$range), field_sd = parameters$field_SD, covariance = "exponential")
}

runtime_check <- function(path) {
  fixture <- fit_vignette_fixture()
  d <- fixture$d
  analysis <- analysis_laplace(d)
  invisible(gaussian_approximation(d, fixture$field_sd, fixture$covariance))
  invisible(tilted_log_moments(d, analysis))
  order <- rep(c("current_first", "tilted_first"), 15)
  results <- lapply(seq_along(order), function(i) {
    current <- function() invisible(gaussian_approximation(d, fixture$field_sd, fixture$covariance))
    tilted <- function() invisible(tilted_approximation(d, analysis))
    time <- function(fun) unname(system.time(fun(), gcFirst = FALSE)[["elapsed"]])
    if (order[i] == "current_first") {
      a <- time(current)
      b <- time(tilted)
    } else {
      b <- time(tilted)
      a <- time(current)
    }
    data.frame(replicate = i, order = order[i], training_nodes = d$n_train,
      assessment_nodes = ncol(d$lower) - d$n_train, current_seconds = a, added_tilt_seconds = b)
  })
  out <- do.call(rbind, results)
  write.csv(out, file.path(path, "tilted-runtime.csv"), row.names = FALSE)
  print(vapply(out[c("current_seconds", "added_tilt_seconds")], median, 0))
}

run_moment_study <- function(path) {
  pkgload::load_all(compile = FALSE, quiet = TRUE)
  cases <- read.csv(file.path(path, "scenarios.csv"), stringsAsFactors = FALSE)
  cases <- subset(cases, cohort == "Laplace" & family %in% c("poisson", "nbinom2", "Gamma"))
  results <- vector("list", nrow(cases))
  for (i in seq_len(nrow(cases))) {
    case <- cases[i, , drop = FALSE]
    d <- simulate_case(case)
    reference <- capture(function() {
      m <- independent_log_moments(d, case$quadrature_order)
      c(m, list(predictive = predictive_moments(m$log_moment, d$family, d$dispersion)))
    })
    current <- capture(function() gaussian_approximation(d, case$field_sd, case$covariance))
    tilted <- capture(function() tilted_approximation(d, analysis_laplace(d)))
    ref <- reference$value
    n <- length(d$test_y)
    value <- function(x, name) if (is.null(x[[name]])) rep(NA_real_, n) else x[[name]]
    accuracy <- function(x, name) {
      if (!isTRUE(ref$converged) || is.null(x)) return(rep(NA_real_, n))
      abs(x[[name]]/ref$predictive[[name]] - 1)
    }
    results[[i]] <- cbind(case[rep(1L, n), c("id", "family", "field_sd", "rho", "visits",
      "layout", "covariance")], data.frame(assessment_row = seq_len(n),
      reference_status = if (is.null(ref)) "failed" else if (ref$converged) "converged" else "unresolved",
      reference_message = reference$message, reference_warnings = reference$warnings,
      quadrature_order = value(ref, "order"), log_refinement_gap = value(ref, "gap"),
      log_proposal_gap = value(ref, "alternate"),
      reference_mean = value(ref$predictive, "mean"), reference_second = value(ref$predictive, "second"),
      reference_variance = value(ref$predictive, "variance"),
      current_mean = value(current$value, "mean"), current_second = value(current$value, "second"),
      current_variance = value(current$value, "variance"),
      current_mean_error = accuracy(current$value, "mean"), current_second_error = accuracy(current$value, "second"),
      current_variance_error = accuracy(current$value, "variance"), current_message = current$message,
      current_warnings = current$warnings, tilted_mean = value(tilted$value, "mean"),
      tilted_second = value(tilted$value, "second"), tilted_variance = value(tilted$value, "variance"),
      tilted_mean_error = accuracy(tilted$value, "mean"), tilted_second_error = accuracy(tilted$value, "second"),
      tilted_variance_error = accuracy(tilted$value, "variance"), tilted_message = tilted$message,
      tilted_warnings = tilted$warnings))
    if (i %% 96L == 0L) cat(i, "of", nrow(cases), "moment cases processed\n")
  }
  all <- do.call(rbind, results)
  rownames(all) <- NULL
  write.csv(all, file.path(path, "tilted-moments.csv"), row.names = FALSE)
  groups <- split(all, list(all$family, all$field_sd), drop = TRUE)
  summary <- do.call(rbind, lapply(groups, function(x) {
    do.call(rbind, lapply(c("mean", "second", "variance"), function(moment) {
      old <- x[[paste0("current_", moment, "_error")]]
      new <- x[[paste0("tilted_", moment, "_error")]]
      keep <- is.finite(old) & is.finite(new)
      q <- function(error) if (any(keep)) unname(quantile(error[keep], c(.5, .95, 1))) else rep(NA_real_, 3L)
      a <- q(old); b <- q(new)
      data.frame(family = x$family[1L], field_sd = x$field_sd[1L], moment = moment,
        rows = nrow(x), paired = sum(keep), reference_unresolved = length(unique(x$id[x$reference_status != "converged"])),
        tilted_failed = length(unique(x$id[nzchar(x$tilted_message)])), current_p50 = a[1L],
        current_p95 = a[2L], current_max = a[3L], tilted_p50 = b[1L], tilted_p95 = b[2L], tilted_max = b[3L])
    }))
  }))
  rownames(summary) <- NULL
  write.csv(summary, file.path(path, "tilted-moments-summary.csv"), row.names = FALSE)
  print(summary, row.names = FALSE)
  runtime_check(path)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  path <- if (length(args)) args[1L] else
    "tools/benchmarks/spatial_field_cv_quadrature/results"
  run_moment_study(path)
}
