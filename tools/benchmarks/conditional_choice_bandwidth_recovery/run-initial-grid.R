source("tools/benchmarks/spatial_field_cv_quadrature/run.R")

choice_design <- function() {
  d <- expand.grid(kernel = c("gaussian", "exp"), n_analysis = c(100L, 450L),
    signal = c("strong_interior", "weak_interior", "near_lower", "null"),
    replicate = seq_len(24L), KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  d$id <- sprintf("choice_%03d", seq_len(nrow(d)))
  d$seed <- 191850L + seq_len(nrow(d))
  d
}

choice_kernel <- function(radial, distances, bandwidth, kernel) {
  q <- distances / bandwidth
  a <- if (kernel == "gaussian") q^2/2 else q
  w <- exp(-a)
  w <- w/sum(w)
  slope <- if (kernel == "gaussian") 2*a else a
  list(value = drop(radial %*% w), derivative = drop(radial %*% (w * (slope - sum(w * slope)))))
}

choice_softmax <- function(eta, used) {
  peak <- pmax(eta[, 1L], eta[, 2L], eta[, 3L], eta[, 4L])
  weight <- exp(eta - peak)
  log_p <- eta - peak - log(rowSums(weight))
  list(nll = -sum(log_p[used]), log_p = log_p, p = exp(log_p))
}

choice_simulate <- function(case) {
  set.seed(case$seed)
  n_set <- case$n_analysis + 160L
  n <- 4L * n_set
  distances <- seq(.0125, .7875, length.out = 32)
  bands <- cbind(exp(-distances/.08), sin(11*distances), cos(19*distances),
    exp(-((distances - .27)/.10)^2), distances)
  bands <- sweep(bands, 2L, c(1, .7, .55, .6, .35), "*")
  random <- matrix(rnorm(n * ncol(bands)), n)
  radial <- random %*% t(bands) + matrix(rnorm(n * length(distances), sd = .15), n)
  length_cov <- .3 * random[, 2L] + rnorm(n, sd = .95)
  turn <- -.2 * random[, 3L] + rnorm(n, sd = .98)
  offset <- runif(n, -.22, .22)
  sigma <- if (case$signal == "near_lower") .075 else .25
  beta <- switch(case$signal, strong_interior = 1.2, weak_interior = .25, near_lower = 1.2, null = 0)
  kernel <- choice_kernel(radial, distances, sigma, case$kernel)
  q <- distances/sigma
  w <- if (case$kernel == "gaussian") exp(-q^2/2) else exp(-q)
  w <- w/sum(w)
  known_sd <- sqrt(sum(drop(crossprod(w, bands))^2) + .15^2 * sum(w^2))
  raw_beta <- c(beta/known_sd, -.43, .34)
  eta <- raw_beta[1L] * kernel$value + raw_beta[2L] * length_cov + raw_beta[3L] * turn + offset
  eta <- matrix(eta, ncol = 4L, byrow = TRUE)
  p <- choice_softmax(eta, cbind(seq_len(n_set), rep(1L, n_set)))$p
  chosen <- apply(p, 1L, function(prob) sample.int(4L, 1L, prob = prob))
  y <- integer(n); y[seq.int(0L, n - 1L, by = 4L) + chosen] <- 1L
  ids <- sprintf("alternative_%05d", seq_len(n))
  data <- data.frame(used = y, choice = rep(seq_len(n_set), each = 4L),
    hab = 0, move_length = length_cov, turn = turn, exposure = exp(offset), row.names = ids)
  raw <- lapply(seq_len(n), function(i) matrix(radial[i, ], ncol = 1L, dimnames = list(NULL, "hab")))
  names(raw) <- ids
  initial <- choice_kernel(radial, distances, .22, case$kernel)$value
  inputs <- structure(list(kernel_dat = data.frame(hab = initial, row.names = ids),
    raw_cov = raw, d_list = setNames(rep(list(distances), n), ids),
    binned = list(nbins = length(distances)), scale_vars = msr_vars(hab = kernel_var("hab")),
    kernel = case$kernel, sigma = .22, shape = NULL, min_D = 60, max_D = 800,
    unit_conv = 1000, resolution = 10, cell_size = c(10, 10), n_cols = 100L,
    n_covs = 1L, cell_data_stored = TRUE, scl_params = list(mean = mean(initial), sd = sd(initial))),
    class = "multiScaleR_data")
  ai <- seq_len(4L * case$n_analysis); ti <- setdiff(seq_len(n), ai)
  list(data = data, inputs = inputs, folds = list(fold = list(analysis = ids[ai], assessment = ids[ti])),
    ai = ai, ti = ti, radial = radial, distances = distances, chosen = chosen, offset = offset,
    known_eta = eta, sigma = sigma, beta = beta, raw_beta = raw_beta, kernel = case$kernel)
}

choice_matrix <- function(x) matrix(x, ncol = 4L, byrow = TRUE)

choice_nuisance <- function(x, offset, chosen) {
  used <- cbind(seq_along(chosen), chosen)
  indicator <- matrix(0, nrow(offset), 4L)
  indicator[used] <- 1
  at <- function(beta) {
    eta <- offset
    for (j in seq_along(x)) eta <- eta + beta[j] * x[[j]]
    value <- choice_softmax(eta, used)
    residual <- value$p - indicator
    gradient <- vapply(x, function(X) sum(residual * X), 0)
    list(value = value$nll, p = value$p, gradient = gradient)
  }
  beta <- numeric(length(x))
  state <- at(beta)
  for (iteration in seq_len(60L)) {
    if (max(abs(state$gradient)) < 1e-7)
      return(list(beta = beta, value = state$value, score = state$gradient,
        residual = state$p - indicator))
    mean <- lapply(x, function(X) rowSums(X * state$p))
    H <- matrix(0, length(x), length(x))
    for (j in seq_along(x)) for (k in seq_len(j)) {
      H[j, k] <- sum(x[[j]] * x[[k]] * state$p) - sum(mean[[j]] * mean[[k]])
      H[k, j] <- H[j, k]
    }
    direction <- drop(solve(H, state$gradient))
    for (j in 0:28) {
      next_beta <- beta - 2^-j * direction
      next_state <- at(next_beta)
      if (next_state$value <= state$value + 16 * .Machine$double.eps * max(1, state$value)) break
    }
    if (j == 28L) stop("Independent nuisance search did not improve.")
    beta <- next_beta
    state <- next_state
  }
  stop("Independent nuisance score did not converge.")
}

choice_reference <- function(d) {
  training <- seq_len(length(d$ai)/4L)
  assess <- setdiff(seq_along(d$chosen), training)
  movement <- c("move_length", "turn")
  x <- lapply(d$data[movement], choice_matrix)
  x_train <- lapply(x, function(X) X[training, , drop = FALSE])
  x_test <- lapply(x, function(X) X[assess, , drop = FALSE])
  offset <- choice_matrix(d$offset)
  cache <- new.env(parent = emptyenv())
  at <- function(log_s) {
    key <- sprintf("%.17g", log_s)
    if (!is.null(cache[[key]])) return(cache[[key]])
    kernel <- choice_kernel(d$radial[d$ai, , drop = FALSE], d$distances, exp(log_s), d$kernel)
    hab <- choice_matrix(kernel$value)
    fit <- choice_nuisance(c(list(hab), x_train), offset[training, , drop = FALSE], d$chosen[training])
    fit$scale_score <- fit$beta[1L] * sum(fit$residual * choice_matrix(kernel$derivative))
    cache[[key]] <- fit[c("beta", "value", "score", "scale_score")]
    cache[[key]]
  }
  grid <- seq(log(.06), log(.8), length.out = 121L)
  fits <- lapply(grid, at)
  nll <- vapply(fits, function(x) x$value, 0)
  s_score <- vapply(fits, function(x) x$scale_score, 0)
  minima <- which(nll[-c(1L, length(nll))] <= nll[-c(length(nll) - 1L, length(nll))] &
    nll[-c(1L, length(nll))] <= nll[-c(1L, 2L)]) + 1L
  refined <- vapply(minima, function(i) optimize(function(s) at(s)$value,
    grid[c(i - 1L, i + 1L)], tol = 1e-8)$minimum, 0)
  candidates <- unique(c(grid[1L], grid[length(grid)], grid[which.min(nll)], refined))
  results <- lapply(candidates, at)
  values <- vapply(results, function(x) x$value, 0)
  best_i <- which.min(values)
  best <- results[[best_i]]
  best_sigma <- exp(candidates[best_i])
  best$sigma <- best_sigma
  best$projected_scale_score <- if (best_sigma <= .06000001) min(0, best$scale_score) else
    if (best_sigma >= .79999999) max(0, best$scale_score) else best$scale_score
  keep <- nll <= best$value + 1.92
  support_regions <- sum(diff(c(FALSE, keep)) == 1L)
  support_span <- if (any(keep)) range(exp(grid[keep])) else c(NA_real_, NA_real_)
  used <- cbind(seq_along(assess), d$chosen[assess])
  evaluate <- function(beta, sigma = NULL) {
    eta <- offset[assess, , drop = FALSE]
    des <- x_test
    if (!is.null(sigma)) {
      hab <- choice_kernel(d$radial[d$ti, , drop = FALSE], d$distances, sigma, d$kernel)$value
      des <- c(list(choice_matrix(hab)), des)
    }
    for (i in seq_along(des)) eta <- eta + beta[i] * des[[i]]
    choice_softmax(eta, used)
  }
  no_hab <- choice_nuisance(x_train, offset[training, , drop = FALSE], d$chosen[training])
  known_scale <- at(log(d$sigma))
  no_hab_loss <- evaluate(no_hab$beta)$nll
  list(best = best, profile = data.frame(bandwidth = 1000 * exp(grid), nll = nll,
      scale_score = s_score, nuisance_score = vapply(fits, function(x) max(abs(x$score)), 0), support = keep),
    refined = data.frame(bandwidth = 1000 * exp(candidates), nll = values,
      scale_score = vapply(results, function(x) x$scale_score, 0)),
    support_regions = support_regions, support_lower = support_span[1L], support_upper = support_span[2L],
    best_loss = evaluate(best$beta, best_sigma)$nll,
    known_scale_loss = if (d$beta == 0) no_hab_loss else evaluate(known_scale$beta, d$sigma)$nll,
    no_hab_loss = no_hab_loss,
    oracle_loss = choice_softmax(d$known_eta[assess, , drop = FALSE], used)$nll,
    evaluate = evaluate, at = at)
}

choice_package <- function(d, starts = NULL) {
  control <- if (is.null(starts)) list() else list(starts = matrix(starts, ncol = 1L))
  fit <- cross_validate_clogit(used ~ hab + move_length + turn + offset(log(exposure)) + strata(choice),
    d$data, d$folds, d$inputs, optim_args = list(profile_control = control),
    keep_fits = TRUE, keep_selection = TRUE)
  if (!fit$complete) return(list(scored = FALSE, reason = fit$errors$fold$message,
    stage = fit$errors$fold$stage, starts = fit$errors$fold$start_results,
    warnings = paste(fit$warnings$fold, collapse = " | ")))
  x <- fit$fits$fold
  b <- stats::coef(x$opt_mod)
  c <- x$scl_params$sd[["hab"]]
  list(scored = TRUE, stage = "", reason = "", warnings = paste(fit$warnings$fold, collapse = " | "),
    starts = x$profile_fit$start_results, winner = x$profile_fit$winner,
    sigma = x$scale_est["hab", "Mean"]/1000, scale_se = x$scale_est["hab", "SE"],
    nll = x$optim_results$value, beta = unname(b[c("hab", "move_length", "turn")]) / c(c, 1, 1),
    loss = sum(fit$predictions$neg_log_score[fit$predictions$eligible]),
    p = choice_matrix(fit$predictions$prediction), scale_boundary = x$diagnostics$conditional_logit$boundary,
    K = fit$selection_records$fold$data$K, nobs = fit$selection_records$fold$data$nobs)
}

append_choice <- function(x, file) {
  if (is.null(x)) return(invisible(NULL))
  a <- file.exists(file)
  write.table(x, file, append = a, col.names = !a, sep = ",", row.names = FALSE)
}

choice_case <- function(case, out, starts = NULL) {
  start <- proc.time()[["elapsed"]]
  d <- choice_simulate(case)
  target <- capture(function() choice_reference(d))
  package <- capture(function() choice_package(d, starts))
  ref <- target$value; fit <- package$value
  get <- function(x, nm) if (is.null(x[[nm]])) NA else x[[nm]]
  if (!is.null(ref)) {
    append_choice(cbind(id = case$id, ref$profile), file.path(out, "profiles.csv"))
    append_choice(cbind(id = case$id, ref$refined), file.path(out, "refined.csv"))
  }
  if (!is.null(fit$starts)) append_choice(cbind(id = case$id, fit$starts), file.path(out, "starts.csv"))
  agreement <- NULL
  if (isTRUE(fit$scored) && !is.null(ref)) {
    nuisance <- ref$at(log(fit$sigma))
    ref_pred <- ref$evaluate(fit$beta, fit$sigma)
    agreement <- list(nll = fit$nll - nuisance$value, p = max(abs(fit$p - ref_pred$p)),
      beta = max(abs(fit$beta - nuisance$beta)),
      scale_score = nuisance$scale_score)
  }
  prediction <- if (isTRUE(fit$scored)) fit$p else matrix(NA_real_, nrow = 160L, ncol = 4L)
  actual <- d$chosen[setdiff(seq_along(d$chosen), seq_len(case$n_analysis))]
  append_choice(data.frame(id = case$id, choice = seq_len(160L), used = actual,
    p1 = prediction[, 1L], p2 = prediction[, 2L], p3 = prediction[, 3L], p4 = prediction[, 4L]),
    file.path(out, "assessment.csv"))
  row <- cbind(case, data.frame(true_bandwidth = if (d$beta == 0) NA_real_ else d$sigma * 1000,
    true_standardized_beta = d$beta, fitted_status = if (isTRUE(fit$scored)) "scored" else "failed",
    fitted_stage = get(fit, "stage"), fitted_reason = if (is.null(fit)) package$message else fit$reason,
    fitted_warnings = paste(package$warnings, get(fit, "warnings")),
    reference_status = if (is.null(ref)) "failed" else "returned", reference_reason = target$message,
    reference_warnings = target$warnings, fitted_bandwidth = get(fit, "sigma") * 1000,
    independent_bandwidth = get(ref$best, "sigma") * 1000,
    fitted_nll = get(fit, "nll"), independent_nll = get(ref$best, "value"),
    reference_constrained_score = get(ref$best, "projected_scale_score"),
    reference_nuisance_score = if (is.null(ref$best$score)) NA_real_ else max(abs(ref$best$score)),
    reference_support_regions = get(ref, "support_regions"), reference_support_lower = get(ref, "support_lower") * 1000,
    reference_support_upper = get(ref, "support_upper") * 1000,
    public_vs_own_profile_nll = get(agreement, "nll"), max_probability_gap = get(agreement, "p"),
    max_raw_coefficient_gap = get(agreement, "beta"), profile_scale_score_at_public = get(agreement, "scale_score"),
    fitted_boundary = get(fit, "scale_boundary"), fitted_winner = get(fit, "winner"),
    eligible_starts = if (is.null(fit$starts)) NA_integer_ else sum(fit$starts$eligible),
    fitted_scale_se_available = is.finite(get(fit, "scale_se")), training_K = get(fit, "K"), training_nobs = get(fit, "nobs"),
    fitted_assessment_nll = get(fit, "loss"), independent_assessment_nll = get(ref, "best_loss"),
    fixed_truth_assessment_nll = get(ref, "known_scale_loss"), oracle_assessment_nll = get(ref, "oracle_loss"),
    no_habitat_assessment_nll = get(ref, "no_hab_loss"), elapsed_seconds = proc.time()[["elapsed"]] - start))
  append_choice(row, file.path(out, "scenarios.csv"))
  row
}

run_choice <- function(out, part, total, limit, primary = NULL) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(out, "scenarios.csv"))) stop("Choose a fresh output directory.")
  pkgload::load_all(compile = FALSE, quiet = TRUE)
  cases <- choice_design()
  starts <- NULL
  if (!is.null(primary)) {
    previous <- read.csv(primary, stringsAsFactors = FALSE)
    stopifnot(nrow(previous) == 384L, setequal(previous$id, cases$id))
    selected <- with(previous, reference_status == "returned" & (fitted_status == "failed" |
      fitted_nll - independent_nll > .0001))
    cases <- cases[cases$id %in% previous$id[selected], ]
    starts <- c(.06, .09, .18, .35, .62, .8)
  }
  cases <- cases[((seq_len(nrow(cases)) - 1L) %% total) + 1L == part, , drop = FALSE]
  cases <- head(cases, min(limit, nrow(cases)))
  write.csv(cases, file.path(out, "design.csv"), row.names = FALSE)
  writeLines(c(capture.output(sessionInfo()), "", "Source MD5",
    capture.output(print(tools::md5sum(c("R/conditional_logit.R", "R/cv_choice.R",
      "tools/benchmarks/conditional_choice_bandwidth_recovery/run.R"))))),
    file.path(out, "session.txt"))
  for (i in seq_len(nrow(cases))) {
    row <- choice_case(cases[i, , drop = FALSE], out, starts)
    if (i %% 8L == 0L || i == nrow(cases))
      cat("part", part, i, "of", nrow(cases), row$id, row$fitted_status, row$reference_status, "\n")
  }
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args)) stop("Supply a fresh output directory.")
  run_choice(args[1L], if (length(args) > 1L) as.integer(args[2L]) else 1L,
    if (length(args) > 2L) as.integer(args[3L]) else 1L,
    if (length(args) > 3L) as.integer(args[4L]) else Inf,
    if (length(args) > 4L) args[5L] else NULL)
}
