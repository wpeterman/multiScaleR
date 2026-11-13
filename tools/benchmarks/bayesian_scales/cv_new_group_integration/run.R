# Run from the package root with Rscript. An optional output directory permits
# independent smoke runs without replacing the archived assessment.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 5L)
  stop("Usage: Rscript run.R [output_dir] [replicates] [posterior_vectors] [draw_budgets] [cases]")
directory <- "tools/benchmarks/bayesian_scales/cv_new_group_integration"
output <- if (length(args)) args[[1L]] else file.path(directory, "results")
n_repeat <- if (length(args) >= 2L) as.integer(args[[2L]]) else 32L
n_vectors <- if (length(args) >= 3L) as.integer(args[[3L]]) else 96L
budgets <- if (length(args) >= 4L) as.integer(strsplit(args[[4L]], ",", fixed = TRUE)[[1L]]) else
  c(1L, 4L, 16L, 64L)
stopifnot(is.finite(n_repeat), n_repeat > 0L, is.finite(n_vectors), n_vectors > 1L,
  length(budgets) > 0L, all(is.finite(budgets)), all(budgets > 0L), !anyDuplicated(budgets))
dir.create(output, recursive = TRUE, showWarnings = FALSE)
pkgload::load_all(compile = FALSE, quiet = TRUE)
ns <- asNamespace("multiScaleR")
cases <- list(
  gaussian_small = list(family = "gaussian", tau = .15, sd = .8, heavy = 0,
    y = c(.2, -.4, .3, .1, -.2, .5, 0, -.1), block = rep("fold", 8)),
  gaussian_large = list(family = "gaussian", tau = 2, sd = .65, heavy = 0,
    y = c(3.1, 2.7, 3.4, .1, -2.2, -2.6, -1.9, -.1), block = rep(c("east", "west"), each = 4)),
  gaussian_extreme_weighted = list(family = "gaussian", tau = 2.5, sd = .4, heavy = .25,
    y = c(6, 5.8, 6.1, .1, -5.5, -5.4, -5.8, -.1), block = rep("fold", 8)),
  poisson_sparse_weighted = list(family = "poisson", tau = 1.25, heavy = .15,
    y = c(0, 0, 1, 0, 0, 1, 0, 0), block = rep(c("east", "west"), each = 4)),
  poisson_extreme_weighted = list(family = "poisson", tau = 1.8, heavy = .25,
    y = c(36, 44, 39, 0, 0, 0, 0, 0), block = rep("fold", 8)),
  gaussian_large_diffuse_weighted = list(family = "gaussian", tau = 2, sd = .65,
    heavy = .14, mass_count = 32L, y = c(3.1, 2.7, 3.4, .1, -2.2, -2.6, -1.9, -.1),
    block = rep(c("east", "west"), each = 4)),
  gaussian_extreme_diffuse_weighted = list(family = "gaussian", tau = 2.5, sd = .4,
    heavy = .14, mass_count = 32L, y = c(6, 5.8, 6.1, .1, -5.5, -5.4, -5.8, -.1),
    block = rep("fold", 8)),
  poisson_extreme_diffuse_weighted = list(family = "poisson", tau = 1.8,
    heavy = .14, mass_count = 32L, y = c(36, 44, 39, 0, 0, 0, 0, 0), block = rep("fold", 8))
)

selected_cases <- if (length(args) == 5L) {
  match(strsplit(args[[5L]], ",", fixed = TRUE)[[1L]], names(cases))
} else seq_len(5L)
stopifnot(length(selected_cases) > 0L, !anyNA(selected_cases), !anyDuplicated(selected_cases))

# Distances, finite-cell sums, held-out groups and posterior vectors are all
# defined independently of the scoring routine. The posterior is synthetic.
distance <- c(.06, .18, .35, .55, .8)
points <- paste0("p", seq_len(8L))
counts <- 1 + outer(seq_along(points), seq_along(distance), function(i, j) (i + j) %% 5)
cover <- outer(seq_along(points), distance, function(i, d) sin(.7 * i + 5 * d) + i * d / 4)
inputs <- list(kernel = "gaussian", unit_conv = 1000, min_D = 35, max_D = 1000,
  binned = list(bin_rep = distance, vsum = list(habitat = cover * counts),
    csum = list(habitat = counts), point_ids = points))
train <- data.frame(key = points[c(1, 2, 3, 4, 1, 3)],
  site = rep(c("known_a", "known_b"), each = 3), count = c(0, 3, 1, 0, 2, 1))
assessment <- data.frame(key = points[c(5, 6, 5, 1, 7, 8, 7, 3)],
  site = rep(c("new_a", "known_a", "new_b", "known_b"), c(3, 1, 3, 1)))
new <- startsWith(assessment$site, "new")
kernel_mean <- function(key, scale) {
  ix <- match(key, points)
  w <- exp(-.5 * (distance / scale)^2)
  drop((cover[ix, , drop = FALSE] * counts[ix, , drop = FALSE]) %*% w /
    (counts[ix, , drop = FALSE] %*% w))
}
prior_scale <- .25
training_reference <- kernel_mean(train$key, prior_scale)
log_mix <- function(log_values, weight) {
  top <- max(log_values)
  top + log(sum(weight * exp(log_values - top)))
}

# With one independent Gaussian response intercept per group, the conditional
# covariance after integrating unseen groups is D + tau^2 ZZ'. For Poisson,
# integrate each scalar intercept in the standard-normal measure. Center and
# scale around the log-integrand mode to resolve narrow rare-event densities.
gaussian_log_density <- function(y, eta, tau, sd, site, unseen) {
  same <- outer(site, site, "==") & outer(unseen, unseen, "&")
  covariance <- diag(sd^2, length(y)) + tau^2 * same
  R <- chol(covariance)
  residual <- backsolve(R, y - eta, transpose = TRUE)
  -.5 * length(y) * log(2 * pi) - sum(log(diag(R))) - .5 * sum(residual^2)
}
poisson_log_density <- function(y, eta, tau, site, unseen, tolerance) {
  known <- sum(dpois(y[!unseen], exp(eta[!unseen]), log = TRUE))
  groups <- split(which(unseen), site[unseen])
  one_group <- function(i) {
    n <- sum(y[i]); rate <- sum(exp(eta[i]))
    constant <- sum(y[i] * eta[i] - lgamma(y[i] + 1))
    mode <- uniroot(function(z) tau * (n - rate * exp(tau * z)) - z,
      c(-30, 30), tol = 1e-13)$root
    width <- 1 / sqrt(1 + tau^2 * rate * exp(tau * mode))
    logf <- function(z) constant + n * tau * z - rate * exp(tau * z) + dnorm(z, log = TRUE)
    height <- logf(mode)
    out <- integrate(function(x) exp(logf(mode + width * x) - height),
      -Inf, Inf, rel.tol = tolerance, abs.tol = tolerance * .01, subdivisions = 1000L)
    stopifnot(out$message == "OK", out$value > 0, is.finite(out$value))
    height + log(width * out$value)
  }
  known + sum(vapply(groups, one_group, numeric(1)))
}

make_case <- function(case) {
  fit <- multiScale_bayes(count ~ habitat + (1 | site), inputs, train,
    row_id = "key", priors = list(scale_median = 250), family = case$family, method = "prepare")
  # Fixed nonuniform grid makes this a deterministic discrete mixture.
  z <- qnorm((seq_len(n_vectors) - .5) / n_vectors)
  permuted <- z[((seq_len(n_vectors) - 1L) * 37L) %% n_vectors + 1L]
  intercept <- (if (case$family == "gaussian") 0 else log(.4)) + .18 * z
  slope <- .3 + .14 * permuted
  scale <- prior_scale * exp(.12 * z)
  tau <- case$tau * exp(.12 * permuted)
  sd <- if (case$family == "gaussian") case$sd * exp(.09 * z) else rep(NA_real_, n_vectors)
  known_z <- cbind(.1 + .12 * z, -.15 + .12 * permuted)
  draws <- list(beta = cbind(intercept, slope), scale_beta = cbind(log(scale)),
    scale_z = matrix(numeric(), n_vectors, 0L),
    log_scale_sd = matrix(numeric(), n_vectors, 0L), response_z = known_z,
    log_response_sd = cbind(log(tau)),
    log_dispersion = if (case$family == "gaussian") cbind(log(sd)) else matrix(numeric(), n_vectors, 0L))
  stopifnot(identical(names(draws), names(fit$parameters)),
    identical(vapply(draws, ncol, integer(1)), lengths(fit$parameters)))
  weight <- rep((1 - case$heavy) / n_vectors, n_vectors)
  n_mass <- if (is.null(case$mass_count)) 1L else case$mass_count
  stopifnot(n_mass <= n_vectors)
  add <- seq_len(n_mass)
  weight[add] <- weight[add] + case$heavy / n_mass
  if (case$heavy > 0) {
    fit$importance <- list(candidates = do.call(cbind, draws), weights = weight)
    # Deliberately unavailable, so the scorer must use weighted candidates.
    fit$draws <- NULL
  } else fit$draws <- draws
  context <- ns$.msr_bayes_prediction_prepare(fit, assessment, inputs, "key",
    "sample", "sample", TRUE, TRUE)
  target <- data.frame(observed = case$y, trials = 1, assessment_group = case$block)
  eta <- vapply(seq_len(n_vectors), function(i) {
    fixed <- intercept[i] + slope[i] * (kernel_mean(assessment$key, scale[i]) -
      mean(training_reference)) / stats::sd(training_reference)
    known <- match(assessment$site, c("known_a", "known_b"))
    fixed[!new] <- fixed[!new] + tau[i] * known_z[i, known[!new]]
    fixed
  }, numeric(nrow(assessment)))
  blocks <- split(seq_len(nrow(target)), factor(case$block, levels = unique(case$block)))
  rows <- as.list(seq_len(nrow(target)))
  normal <- identical(case$family, "gaussian")
  marginalize <- function(indices, tol) {
    vapply(seq_len(n_vectors), function(i) {
      y <- target$observed[indices]; lp <- eta[indices, i]
      site <- assessment$site[indices]; unseen <- new[indices]
      if (normal) gaussian_log_density(y, lp, tau[i], sd[i], site, unseen) else
        poisson_log_density(y, lp, tau[i], site, unseen, tol)
    }, numeric(1))
  }
  reference <- function(groups, tol) -vapply(groups, function(g)
    log_mix(marginalize(g, tol), weight), numeric(1))
  truth <- list(joint = reference(blocks, 1e-11), marginal = reference(rows, 1e-11))
  check <- c(truth$joint - reference(blocks, 1e-9), truth$marginal - reference(rows, 1e-9))
  stopifnot(all(abs(check) < 1e-8))
  if (normal) {
    # Direct univariate identities independently check the marginal covariance.
    for (j in seq_along(rows)) {
      single <- dnorm(case$y[j], eta[j, ], sqrt(sd^2 + new[j] * tau^2), log = TRUE)
      stopifnot(abs(truth$marginal[j] + log_mix(single, weight)) < 1e-11)
    }
  }
  true_response <- if (normal) drop(eta %*% weight) else
    drop(exp(eta + outer(new, tau^2 / 2)) %*% weight)
  list(fit = fit, context = context, target = target, truth = truth, response = true_response,
    weight_ess = 1 / sum(weight^2), max_weight = max(weight), quadrature_delta = max(abs(check)),
    blocks = blocks)
}

run_case <- function(i) {
  name <- names(cases)[i]
  case <- cases[[i]]
  ref <- make_case(case)
  measures <- c("fold", rep("joint_block", length(ref$blocks)), rep("marginal_row", nrow(assessment)),
    rep("expected_response", nrow(assessment)))
  targets <- c("pooled", names(ref$blocks), seq_len(nrow(assessment)), seq_len(nrow(assessment)))
  denom <- c(nrow(assessment), lengths(ref$blocks), rep(1L, nrow(assessment) * 2L))
  truth <- c(sum(ref$truth$joint), ref$truth$joint, ref$truth$marginal, ref$response) / denom
  result <- vector("list", length(budgets) * n_repeat)
  for (b in seq_along(budgets)) {
    budget <- budgets[b]
    for (r in seq_len(n_repeat)) {
      seed <- 321000L + i * 10000L + b * 100L + r +
        if (n_vectors == 96L) 0L else n_vectors * 100L
      start <- proc.time()[["elapsed"]]
      score <- tryCatch(ns$.msr_bayes_cv_scores(ref$fit, ref$context, ref$target, budget, seed),
        error = identity)
      elapsed <- proc.time()[["elapsed"]] - start
      failure <- inherits(score, "error")
      values <- if (failure) rep(NA_real_, length(truth)) else
        c(sum(score$joint), score$joint, score$marginal_neg_log_score, score$prediction) / denom
      result[[(b - 1L) * n_repeat + r]] <- data.frame(case = name, family = case$family,
        weight_ess = ref$weight_ess, max_weight = ref$max_weight, quadrature_delta = ref$quadrature_delta,
        budget = budget, replicate = r, seed = seed, measure = measures, target = targets,
        reference = truth, estimate = values, error = values - truth, elapsed = elapsed,
        status = if (failure) "failed" else ifelse(is.finite(values), "scored", "nonfinite"),
        message = if (failure) conditionMessage(score) else "", row.names = NULL)
    }
    message(name, " at ", budget, " new draws complete")
  }
  do.call(rbind, result)
}

results <- do.call(rbind, lapply(selected_cases, run_case))
connection <- gzfile(file.path(output, "all-results.csv.gz"), "wt")
write.csv(results, connection, row.names = FALSE)
close(connection)
groups <- split(results, list(results$case, results$budget, results$measure), drop = TRUE)
summary <- do.call(rbind, lapply(groups, function(x) {
  good <- x$status == "scored"
  e <- x$error[good]
  finite_stat <- function(f) if (length(e)) f(e) else NA_real_
  data.frame(case = x$case[1L], family = x$family[1L], budget = x$budget[1L],
    measure = x$measure[1L], weight_ess = x$weight_ess[1L], max_weight = x$max_weight[1L],
    quadrature_delta = x$quadrature_delta[1L], n = nrow(x), failure = sum(!good),
    bias = finite_stat(mean), sd = finite_stat(stats::sd),
    rmse = finite_stat(function(e) sqrt(mean(e^2))),
    abs_q95 = finite_stat(function(e) unname(quantile(abs(e), .95))),
    max_abs = finite_stat(function(e) max(abs(e))),
    fraction_abs_over_0.1 = finite_stat(function(e) mean(abs(e) > .1)))
}))
summary <- summary[order(summary$case, summary$measure, summary$budget), ]
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output, "session-info.txt"))
print(summary[summary$measure == "fold", ], row.names = FALSE)
