# Run from the package root. The optional argument names a fresh results folder.
source("tools/benchmarks/spatial_field_cv_quadrature/run.R")

large_design <- function() {
  settings <- data.frame(setting = c("weak_sparse", "strong_sparse", "strong_repeated"),
    n_analysis_nodes = c(8L, 16L, 16L), n_assessment_nodes = c(4L, 6L, 6L),
    visits = c(1L, 1L, 4L), rho = c(.25, .75, .75))
  count <- expand.grid(family = c("poisson", "nbinom2"), field_sd = c(.5, 1.8),
    setting = settings$setting, replicate = 1:2, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  other <- expand.grid(family = c("binomial", "Gamma"), field_sd = 1.8,
    setting = settings$setting[2:3], replicate = 1L, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  cases <- rbind(count, other)
  cases <- cbind(cases, settings[match(cases$setting, settings$setting), -1L])
  cases$covariance <- ifelse(cases$replicate == 2L |
      cases$family %in% c("binomial", "Gamma") & cases$setting == "strong_repeated",
    "gaussian", "exponential")
  cases$id <- sprintf("larger_%03d", seq_len(nrow(cases)))
  cases$seed <- 913040L + seq_len(nrow(cases))
  cases$cohort <- "larger"
  cases$small_id <- NA_character_
  small <- read.csv("tools/benchmarks/spatial_field_cv_quadrature/results/design.csv",
    stringsAsFactors = FALSE)
  small <- small[match(c("field_0357", "field_0167", "field_0142", "field_0792"), small$id), ]
  n_analysis <- ifelse(small$layout == "2+1", 2L, 1L)
  n_assessment <- ifelse(small$layout == "1+2", 2L, 1L)
  controls <- data.frame(family = small$family, field_sd = small$field_sd,
    setting = small$layout, replicate = small$replicate, n_analysis_nodes = n_analysis,
    n_assessment_nodes = n_assessment, visits = small$visits, rho = small$rho,
    covariance = small$covariance, id = sprintf("control_%03d", seq_len(nrow(small))),
    seed = small$seed, cohort = "quadrature_control", small_id = small$id)
  list(cases = rbind(cases, controls), small = small)
}

simulate_larger <- function(case) {
  set.seed(case$seed)
  n_train <- case$n_analysis_nodes
  count <- n_train + case$n_assessment_nodes
  xy <- if (n_train == 8L) as.matrix(expand.grid(x = 0:3, y = 0:2)) else
    as.matrix(expand.grid(x = 0:5, y = 0:3))[seq_len(count), ]
  assess <- if (n_train == 8L) c(3, 6, 9, 12) else c(5, 8, 11, 14, 17, 20)
  xy <- xy[c(setdiff(seq_len(count), assess), assess), , drop = FALSE]
  power <- if (case$covariance == "exponential") 1 else 2
  range <- (-1/log(case$rho))^(1/power)
  covariance <- case$field_sd^2 * exp(-(as.matrix(dist(xy))/range)^power)
  lower <- t(chol(covariance))
  groups <- rep(seq_len(count), each = case$visits)
  eta <- if (case$family == "binomial") -.35 + rnorm(length(groups), sd = .65) else
    .35 + rnorm(length(groups), sd = .45)
  link <- eta + drop(lower %*% rnorm(count))[groups]
  dispersion <- switch(case$family, nbinom2 = 3, Gamma = 4, 1)
  y <- switch(case$family, poisson = rpois(length(link), exp(link)),
    binomial = rbinom(length(link), 1, plogis(link)),
    nbinom2 = rnbinom(length(link), size = dispersion, mu = exp(link)),
    Gamma = rgamma(length(link), shape = dispersion, scale = exp(link)/dispersion))
  train <- which(groups <= n_train)
  test <- which(groups > n_train)
  list(train_y = y[train], test_y = y[test], train_eta = eta[train], test_eta = eta[test],
    train_xy = xy[groups[train], , drop = FALSE], test_xy = xy[groups[test], , drop = FALSE],
    groups = groups, lower = lower, n_train = n_train, family = case$family,
    dispersion = dispersion, range = range)
}

score_functions <- function(d) {
  training <- seq_along(d$train_y)
  assessment <- length(training) + seq_along(d$test_y)
  list(train = function(field)
    sum(response_logdens(d$train_y, d$train_eta + field[d$groups[training]], d$family, d$dispersion)),
    test = function(field)
      sum(response_logdens(d$test_y, d$test_eta + field[d$groups[assessment]], d$family, d$dispersion)))
}

build_proposal <- function(lower, scores, temperature) {
  dimension <- ncol(lower)
  target <- function(z) {
    field <- drop(lower %*% z)
    nll <- sum(z^2)/2 - scores$train(field)
    if (temperature > 0) nll <- nll - temperature * scores$test(field)
    nll
  }
  fitting <- capture(function() {
    fit <- optim(numeric(dimension), target, method = "BFGS",
      hessian = TRUE, control = list(reltol = 1e-11, maxit = 400))
    if (fit$convergence != 0L || !is.finite(fit$value)) stop("Numerical proposal optimization failed.")
    list(center = fit$par, factor = 1.15 * t(chol(chol2inv(chol(fit$hessian)))))
  })
  if (is.null(fitting$value)) {
    list(center = numeric(dimension), factor = diag(dimension),
      method = "prior", message = fitting$message, warnings = fitting$warnings)
  } else {
    c(fitting$value, list(method = "local_Gaussian", message = "", warnings = fitting$warnings))
  }
}

# The ellipse rotates draws under the proposal N(0, I). The slice threshold
# uses the exact Gaussian prior and normalized conditional response densities
# divided by that proposal density. The numerical proposal does not set the
# posterior target.
ess_chain <- function(d, lower, scores, temperature, proposal, seed, warmup, draws) {
  set.seed(seed)
  dimension <- ncol(lower)
  z_factor <- proposal$factor
  field_center <- drop(lower %*% proposal$center)
  field_factor <- lower %*% z_factor
  residual <- function(v, z, field) {
    value <- (sum(v^2) - sum(z^2))/2 + scores$train(field)
    if (temperature > 0) value <- value + temperature * scores$test(field)
    value
  }
  v <- rnorm(dimension)
  z <- proposal$center + drop(z_factor %*% v)
  field <- field_center + drop(field_factor %*% v)
  current <- residual(v, z, field)
  if (!is.finite(current)) stop("The posterior chain starting target is not finite.")
  prior_v <- matrix(rnorm(dimension * (warmup + draws)), nrow = dimension)
  prior_z <- z_factor %*% prior_v
  prior_field <- field_factor %*% prior_v
  states <- matrix(NA_real_, draws, dimension)
  n_proposals <- integer(draws)
  for (i in seq_len(warmup + draws)) {
    threshold <- current - rexp(1L)
    angle <- runif(1L, 0, 2 * pi)
    left <- angle - 2 * pi
    right <- angle
    for (j in seq_len(1000L)) {
      ca <- cos(angle)
      sa <- sin(angle)
      next_v <- ca * v + sa * prior_v[, i]
      next_z <- proposal$center + ca * (z - proposal$center) + sa * prior_z[, i]
      next_field <- field_center + ca * (field - field_center) + sa * prior_field[, i]
      value <- residual(next_v, next_z, next_field)
      if (is.finite(value) && value >= threshold) break
      if (angle < 0) left <- angle else right <- angle
      angle <- runif(1L, left, right)
    }
    if (!is.finite(value) || value < threshold) stop("Elliptical slice search did not find a valid point.")
    v <- next_v
    z <- next_z
    field <- next_field
    current <- value
    if (i > warmup) {
      states[i - warmup, ] <- field
      n_proposals[i - warmup] <- j
    }
  }
  if (temperature == 0) {
    train <- seq_len(d$n_train)
    test <- setdiff(seq_len(ncol(d$lower)), train)
    cond_map <- forwardsolve(d$lower[train, train, drop = FALSE], t(states))
    conditional <- t(d$lower[test, train, drop = FALSE] %*% cond_map)
    innovation <- matrix(rnorm(length(test) * draws), ncol = length(test))
    states <- cbind(states, conditional + innovation %*% t(d$lower[test, test, drop = FALSE]))
  }
  list(field = states, log_assessment = apply(states, 1L, scores$test),
    average_proposals = mean(n_proposals))
}

split_rank_rhat <- function(values) {
  half <- nrow(values) %/% 2L
  split <- cbind(values[seq_len(half), , drop = FALSE],
    values[seq.int(nrow(values) - half + 1L, nrow(values)), , drop = FALSE])
  rank_rhat <- function(x) {
    ranks <- rank(x, ties.method = "average")
    normal <- matrix(qnorm((ranks - 3/8)/(length(ranks) + 1/4)), nrow = half)
    within <- mean(apply(normal, 2L, var))
    if (within == 0) return(if (var(colMeans(normal)) == 0) 1 else Inf)
    sqrt(((half - 1)/half * within + var(colMeans(normal)))/within)
  }
  max(rank_rhat(split), rank_rhat(abs(split - median(split))))
}

estimate_bridge <- function(left, right, delta) {
  x <- delta * left
  y <- delta * right
  bounds <- range(c(x, y)) + c(-20, 20)
  equation <- function(log_ratio) mean(plogis(x - log_ratio)) - mean(plogis(log_ratio - y))
  log_ratio <- uniroot(equation, bounds, tol = 1e-10)$root
  list(log_ratio = log_ratio,
    overlap = min(mean(plogis(x - log_ratio)), mean(plogis(log_ratio - y))))
}

sample_temperatures <- function(d, seed, temperatures, warmup, draws, chains) {
  scores <- score_functions(d)
  n <- seq_len(d$n_train)
  states <- vector("list", length(temperatures))
  for (i in seq_along(temperatures)) {
    temperature <- temperatures[i]
    lower <- if (temperature == 0) d$lower[n, n, drop = FALSE] else d$lower
    proposal <- build_proposal(lower, scores, temperature)
    states[[i]] <- lapply(seq_len(chains), function(j)
      ess_chain(d, lower, scores, temperature, proposal, seed + 17L * i + j, warmup, draws))
    chain_columns <- lapply(seq_len(ncol(states[[i]][[1L]]$field)), function(k)
      do.call(cbind, lapply(states[[i]], function(x) x$field[, k])))
    la <- do.call(cbind, lapply(states[[i]], function(x) x$log_assessment))
    rhat <- vapply(c(chain_columns, list(la)), split_rank_rhat, 0)
    attr(states[[i]], "diagnostics") <- data.frame(temperature = temperature,
      max_rhat = max(rhat), mean_proposals = mean(vapply(states[[i]], function(x) x$average_proposals, 0)),
      proposal = proposal$method, proposal_message = proposal$message, proposal_warnings = proposal$warnings)
  }
  states
}

monte_carlo_replicates <- function(draws, chains, batch = 100L) {
  list(full = list(rows = seq_len(draws), chains = seq_len(chains)),
    individual = lapply(seq_len(chains), function(i) list(rows = seq_len(draws), chains = i)),
    batches = lapply(seq_len(chains * (draws %/% batch)), function(i) {
      group <- (i - 1L) %% (draws %/% batch)
      list(rows = group * batch + seq_len(batch), chains = (i - 1L) %/% (draws %/% batch) + 1L)
    }))
}

reference_score_mc <- function(states, temperatures, subsets) {
  likelihoods <- lapply(states, function(x) do.call(cbind, lapply(x, function(z) z$log_assessment)))
  estimate <- function(s) {
    loglik <- lapply(likelihoods, function(x) as.vector(x[s$rows, s$chains, drop = FALSE]))
    parts <- lapply(seq_along(diff(temperatures)), function(i)
      estimate_bridge(loglik[[i]], loglik[[i + 1L]], temperatures[i + 1L] - temperatures[i]))
    list(loss = -sum(vapply(parts, function(x) x$log_ratio, 0)),
      overlap = vapply(parts, function(x) x$overlap, 0))
  }
  point <- estimate(subsets$full)
  se <- function(groups) sd(vapply(groups, function(x) estimate(x)$loss, 0))/sqrt(length(groups))
  c(point, list(mcse = max(se(subsets$batches), se(subsets$individual))))
}

reference_moments_mc <- function(d, states, subsets) {
  train <- seq_len(d$n_train)
  test <- setdiff(seq_len(ncol(d$lower)), train)
  assessment <- length(d$train_y) + seq_along(d$test_y)
  nodes <- d$groups[assessment] - d$n_train
  covariance <- tcrossprod(d$lower)
  projection <- solve(covariance[train, train, drop = FALSE], covariance[train, test, drop = FALSE])
  innovation <- diag(covariance[test, test, drop = FALSE] - t(projection) %*% covariance[train, test, drop = FALSE])
  moments <- lapply(states[[1L]], function(chain) {
    fixed <- sweep((chain$field[, train, drop = FALSE] %*% projection)[, nodes, drop = FALSE],
      2L, d$test_eta, "+")
    if (d$family == "binomial") {
      rule <- normal_rule(48L)
      mean <- vapply(seq_along(nodes), function(i) drop(
        plogis(outer(fixed[, i], sqrt(innovation[nodes[i]]) * rule$nodes, "+")) %*% rule$weights),
        numeric(nrow(fixed)))
      return(list(mean = mean))
    }
    list(mean = exp(sweep(fixed, 2L, innovation[nodes]/2, "+")),
      second = exp(sweep(2 * fixed, 2L, 2 * innovation[nodes], "+")))
  })
  estimate <- function(s) {
    moments <- moments[s$chains]
    average <- function(name) colMeans(do.call(rbind,
      lapply(moments, function(x) x[[name]][s$rows, , drop = FALSE])))
    mean <- average("mean")
    if (d$family == "binomial") return(cbind(mean = mean, variance = mean * (1 - mean)))
    second <- average("second")
    noise <- switch(d$family, poisson = mean, nbinom2 = mean + second/d$dispersion,
      Gamma = second/d$dispersion)
    cbind(mean = mean, variance = noise + second - mean^2)
  }
  point <- estimate(subsets$full)
  se <- function(groups) apply(simplify2array(lapply(groups, estimate)), c(1L, 2L), sd) / sqrt(length(groups))
  list(point = point, mcse = pmax(se(subsets$batches), se(subsets$individual)))
}

run_larger <- function(out, limit) {
  scenario_file <- file.path(out, "scenarios.csv")
  prediction_file <- file.path(out, "predictions.csv")
  temperature_file <- file.path(out, "temperatures.csv")
  if (any(file.exists(c(scenario_file, prediction_file, temperature_file))))
    stop("Use a fresh output directory to retain previous records.")
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  pkgload::load_all(compile = FALSE, quiet = TRUE)
  stopifnot(requireNamespace("glmmTMB", quietly = TRUE))
  design <- large_design()
  cases <- head(design$cases, min(limit, nrow(design$cases)))
  quadrature_path <- "tools/benchmarks/spatial_field_cv_quadrature/results"
  small_scores <- read.csv(file.path(quadrature_path, "scenarios.csv"), stringsAsFactors = FALSE)
  small_moments <- read.csv(file.path(quadrature_path, "tilted-moments.csv"), stringsAsFactors = FALSE)
  small_means <- read.csv(file.path(quadrature_path, "predictions.csv"), stringsAsFactors = FALSE)
  temperatures <- seq(0, 1, by = .25)
  chains <- 3L
  warmup <- 350L
  draws <- 1600L
  subsets <- monte_carlo_replicates(draws, chains)
  write.csv(cases, file.path(out, "design.csv"), row.names = FALSE)
  writeLines(c(capture.output(sessionInfo()), "",
    sprintf("Temperatures %s, %s chains, %s warmup, %s draws",
      paste(temperatures, collapse = ", "), chains, warmup, draws),
    "", capture.output(print(tools::md5sum(c("R/field_cv_laplace.R",
      "tools/benchmarks/spatial_field_cv_larger/run.R",
      "tools/benchmarks/spatial_field_cv_quadrature/run.R"))))),
    file.path(out, "session.txt"))
  append_rows <- function(x, file) {
    append <- file.exists(file)
    write.table(x, file, sep = ",", append = append, col.names = !append, row.names = FALSE)
  }
  cat("Executing", nrow(cases), "larger-field and reference-control scenarios\n")
  for (i in seq_len(nrow(cases))) {
    started <- proc.time()[["elapsed"]]
    case <- cases[i, , drop = FALSE]
    d <- if (case$cohort == "larger") simulate_larger(case) else
      simulate_case(design$small[match(case$small_id, design$small$id), , drop = FALSE])
    evaluated <- capture(function() .msr_field_laplace_block(d$train_y, d$test_y,
      d$train_eta, d$test_eta, d$train_xy, d$test_xy, case$covariance, case$field_sd,
      d$range, d$family, d$dispersion, tilted_log_moments = TRUE))
    sampled <- capture(function() sample_temperatures(d, case$seed + 651741L,
      temperatures, warmup, draws, chains))
    ref_score <- ref_moment <- list(value = NULL, message = sampled$message, warnings = sampled$warnings)
    max_rhat <- analysis_rhat <- overlap <- NA_real_
    if (!is.null(sampled$value)) {
      states <- sampled$value
      metadata <- do.call(rbind, lapply(states, attr, "diagnostics"))
      append_rows(cbind(id = case$id, metadata), temperature_file)
      max_rhat <- max(metadata$max_rhat)
      analysis_rhat <- metadata$max_rhat[1L]
      ref_score <- capture(function() reference_score_mc(states, temperatures, subsets))
      ref_moment <- capture(function() reference_moments_mc(d, states, subsets))
      if (!is.null(ref_score$value)) overlap <- min(ref_score$value$overlap)
    }
    n <- length(d$test_y)
    estimate <- function(x, field) if (is.null(x[[field]])) NA_real_ else x[[field]]
    ref_loss <- estimate(ref_score$value, "loss")
    loss_se <- estimate(ref_score$value, "mcse")
    fit_loss <- estimate(evaluated$value$joint, "neg_log_score")
    score_resolved <- is.finite(max_rhat) && max_rhat <= 1.05 && is.finite(overlap) &&
      overlap >= .02 && is.finite(loss_se) && loss_se/n <= .015
    expected <- if (case$cohort == "quadrature_control") small_scores$reference_loss[
      match(case$small_id, small_scores$id)] else NA_real_
    moment <- ref_moment$value
    rvalue <- function(field, column) if (is.null(moment[[field]]))
      rep(NA_real_, n) else moment[[field]][, column]
    quad_mean <- quad_var <- rep(NA_real_, n)
    if (case$cohort == "quadrature_control") {
      quad_mean <- small_means$reference_mean[small_means$id == case$small_id]
      quad_var <- if (case$family == "binomial") quad_mean * (1 - quad_mean) else
        small_moments$reference_variance[small_moments$id == case$small_id]
    }
    ref_mean <- rvalue("point", "mean")
    ref_variance <- rvalue("point", "variance")
    mean_se <- rvalue("mcse", "mean")
    variance_se <- rvalue("mcse", "variance")
    mean_threshold <- if (case$family == "binomial") .01 else .03 * abs(ref_mean)
    ok <- is.finite(analysis_rhat) && analysis_rhat <= 1.05
    mean_resolved <- ok & is.finite(mean_se) & mean_se <= mean_threshold
    variance_resolved <- ok & is.finite(variance_se) & variance_se <= .1 * ref_variance
    fitted <- evaluated$value
    preds <- data.frame(id = case$id, assessment_row = seq_len(n), observed = d$test_y,
      reference_mean = ref_mean, reference_mean_mcse = mean_se, mean_resolved = mean_resolved,
      reference_variance = ref_variance, reference_variance_mcse = variance_se,
      variance_resolved = variance_resolved, fitted_mean = estimate(fitted, "prediction"),
      fitted_variance = estimate(fitted, "observation_predictive_SE")^2,
      quadrature_mean = quad_mean, quadrature_variance = quad_var)
    append_rows(preds, prediction_file)
    append_rows(cbind(case, data.frame(n_analysis = length(d$train_y), n_assessment = n,
      max_analysis_response = max(d$train_y), max_assessment_response = max(d$test_y),
      reference_loss = ref_loss, reference_loss_mcse = loss_se, joint_resolved = score_resolved,
      quadrature_loss = expected, fitted_loss = fit_loss, max_rhat = max_rhat,
      analysis_rhat = analysis_rhat, minimum_overlap = overlap,
      means_resolved = sum(mean_resolved), variances_resolved = sum(variance_resolved),
      sampling_message = sampled$message, sampling_warnings = sampled$warnings,
      reference_score_message = ref_score$message, reference_moment_message = ref_moment$message,
      reference_warnings = paste(ref_score$warnings, ref_moment$warnings),
      evaluation_status = if (is.null(fitted)) "failed" else "returned",
      evaluation_message = evaluated$message, evaluation_warnings = evaluated$warnings,
      elapsed_seconds = proc.time()[["elapsed"]] - started)), scenario_file)
    cat(i, "of", nrow(cases), case$id, "score resolved", score_resolved,
      "mean", sum(mean_resolved), "variance", sum(variance_resolved), "\n")
  }
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  out <- if (length(args)) args[1L] else "tools/benchmarks/spatial_field_cv_larger/results"
  limit <- if (length(args) > 1L) as.integer(args[2L]) else Inf
  run_larger(out, limit)
}
