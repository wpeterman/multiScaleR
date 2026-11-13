# Run from the package source root. Outputs are plain CSV plus run metadata.
output <- "tools/benchmarks/mixed_cv_replicated"
pkgload::load_all(compile = NA, quiet = TRUE)
stopifnot(requireNamespace("glmmTMB", quietly = TRUE))
dir.create(output, showWarnings = FALSE, recursive = TRUE)

models <- c(gaussian_intercept = "gaussian", gaussian_slope = "gaussian",
  poisson = "poisson", bernoulli = "binomial", grouped_binomial = "binomial",
  nbinom2 = "nbinom2", Gamma = "Gamma")
scenarios <- c(weak = .12, strong = 1.15, sparse = .65, extreme = .75)
seeds <- c(21L, 39L, 72L, 105L)
targets <- c("new_group", "existing_group")
plan <- expand.grid(model = names(models), scenario = names(scenarios), seed = seeds,
  prediction = targets, KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
plan$id <- sprintf("fold%03d", seq_len(nrow(plan)))
utils::write.csv(plan, file.path(output, "plan.csv"), row.names = FALSE)
sources <- c("R/cross_validate_mixed.R", "R/mixed_cv_density.R", "R/mixed_cv_comparison.R")
writeLines(c(paste("date UTC", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  paste("R", getRversion()), paste("glmmTMB", packageVersion("glmmTMB")),
  paste("TMB", packageVersion("TMB")), capture.output(tools::md5sum(sources)),
  paste("benchmark MD5", tools::md5sum(file.path(output, "run.R")))),
  file.path(output, "run-metadata.txt"))

simulate <- function(model, scenario, seed) {
  set.seed(seed)
  group <- factor(rep(seq_len(18L), each = 6L))
  visit <- rep(seq_len(6L), nlevels(group))
  n <- length(group)
  tau <- scenarios[[scenario]]
  u <- stats::rnorm(nlevels(group))
  b <- tau * u
  slope <- .75 * tau * (.35 * u + sqrt(1 - .35^2) * stats::rnorm(nlevels(group)))
  x <- stats::rnorm(n)
  exposure <- exp(stats::runif(n, -.3, .3))
  intercept <- if (scenario == "sparse") -2.1 else .6
  eta <- intercept + .55 * x + b[group] + log(exposure)
  trials <- switch(model, grouped_binomial = rep(c(4L, 12L, 25L), length.out = n),
    rep(1L, n))
  y <- switch(model,
    gaussian_intercept = eta + stats::rnorm(n, sd = .65),
    gaussian_slope = eta + slope[group] * x + stats::rnorm(n, sd = .65),
    poisson = stats::rpois(n, exp(eta)),
    bernoulli = stats::rbinom(n, trials, stats::plogis(eta - .8)),
    grouped_binomial = stats::rbinom(n, trials, stats::plogis(eta - .8)),
    nbinom2 = stats::rnbinom(n, size = .9, mu = exp(eta)),
    Gamma = stats::rgamma(n, shape = .65, scale = exp(eta) / .65))
  if (scenario == "extreme") {
    selected <- visit == 6L & as.integer(group) <= 4L
    y[selected] <- switch(model,
      gaussian_intercept = y[selected] + 10,
      gaussian_slope = y[selected] - 10,
      poisson = y[selected] + 450L,
      bernoulli = 1L - y[selected],
      grouped_binomial = ifelse(y[selected] > trials[selected] / 2, 0L, trials[selected]),
      nbinom2 = y[selected] + 450L,
      Gamma = y[selected] + 180)
  }
  data.frame(y = y, successes = y, failures = trials - y, trials = trials,
    group = group, x = x, visit = visit, exposure = exposure,
    row.names = sprintf("g%02dv%02d", as.integer(group), visit))
}

fold_rows <- function(d, scenario, prediction) {
  if (prediction == "new_group") {
    assessment <- as.integer(d$group) <= 4L
    analysis <- !assessment
  } else {
    assessment <- d$visit >= 5L
    analysis <- d$visit <= 4L
    if (scenario == "sparse") analysis <- analysis &
      (as.integer(d$group) > 9L | d$visit == 1L)
  }
  list(analysis = rownames(d)[analysis], assessment = rownames(d)[assessment])
}

form <- function(model) switch(model,
  gaussian_slope = y ~ I(scale(x)) + offset(log(exposure)) + (1 + x | group),
  grouped_binomial = cbind(successes, failures) ~ I(scale(x)) + offset(log(exposure)) + (1 | group),
  y ~ I(scale(x)) + offset(log(exposure)) + (1 | group))

gaussian_reference <- function(fit, train, test, eta, slope) {
  covariance <- glmmTMB::VarCorr(fit)$cond$group
  sigma2 <- stats::sigma(fit)^2
  random_design <- function(d) if (slope) cbind(1, d$x) else matrix(1, nrow(d), 1L)
  Ztest <- random_design(test)
  mean <- eta(test)
  cov <- Ztest %*% covariance %*% t(Ztest) + diag(sigma2, nrow(test))
  if (nrow(train)) {
    Ztrain <- random_design(train)
    train_cov <- Ztrain %*% covariance %*% t(Ztrain) + diag(sigma2, nrow(train))
    between <- Ztest %*% covariance %*% t(Ztrain)
    mean <- mean + drop(between %*% solve(train_cov, train$y - eta(train)))
    cov <- cov - between %*% solve(train_cov, t(between))
  }
  residual <- test$y - mean
  score <- (nrow(test) * log(2 * pi) + as.numeric(determinant(cov, logarithm = TRUE)$modulus) +
    drop(crossprod(residual, solve(cov, residual)))) / 2
  if (!is.finite(score) || any(!is.finite(mean))) stop("Nonfinite Gaussian reference.")
  list(mean = mean, score = score, refinement = NA_real_)
}

# Composite eight-point Gauss-Legendre quadrature on the standardized Normal
# effect. Both finite ranges and resolutions are fixed in advance. The native
# implementation uses a mode-centered adaptive infinite-interval integral.
panel_rule <- function(limit, step) {
  J <- matrix(0, 8L, 8L)
  i <- seq_len(7L)
  J[cbind(i, i + 1L)] <- J[cbind(i + 1L, i)] <- i / sqrt(4 * i^2 - 1)
  rule <- eigen(J, symmetric = TRUE)
  centers <- seq(-limit + step / 2, limit - step / 2, by = step)
  nodes <- as.vector(outer(rule$values * step / 2, centers, "+"))
  weights <- rep(rule$vectors[1L, ]^2 * step, length(centers))
  list(z = nodes, log_weight = log(weights) + stats::dnorm(nodes, log = TRUE))
}
quadrature <- list(base = panel_rule(16, .125), fine = panel_rule(20, .0625))
logsum <- function(x) {
  high <- max(x)
  if (!is.finite(high)) return(high)
  high + log(sum(exp(x - high)))
}

glmm_reference <- function(fit, train, test, eta, family) {
  tau <- sqrt(glmmTMB::VarCorr(fit)$cond$group[1L, 1L])
  phi <- stats::sigma(fit)
  likelihood <- function(data, effect) {
    if (!nrow(data)) return(rep(0, length(effect)))
    link <- outer(eta(data), effect, "+")
    y <- rep(data$y, ncol(link))
    mu <- exp(link)
    ld <- switch(family,
      poisson = stats::dpois(y, mu, log = TRUE),
      nbinom2 = stats::dnbinom(y, size = phi, mu = mu, log = TRUE),
      Gamma = stats::dgamma(y, shape = 1 / phi^2, scale = mu * phi^2, log = TRUE),
      binomial = {
        n <- rep(data$trials, ncol(link))
        # Complement both the outcome and probability where the latter is
        # large, so the native dbinom call retains tiny failure probabilities.
        flipped <- as.vector(link) > 0
        stats::dbinom(ifelse(flipped, n - y, y), n,
          stats::plogis(-abs(link)), log = TRUE)
      })
    colSums(matrix(ld, nrow(data), length(effect)))
  }
  integrate <- function(rule) {
    effect <- tau * rule$z
    weighted_train <- rule$log_weight + likelihood(train, effect)
    base <- logsum(weighted_train)
    score <- base - logsum(weighted_train + likelihood(test, effect))
    links <- outer(eta(test), effect, "+")
    log_mean <- if (family == "binomial") stats::plogis(links, log.p = TRUE) else links
    mean <- vapply(seq_len(nrow(test)), function(i)
      exp(logsum(weighted_train + log_mean[i, ]) - base), numeric(1))
    list(mean = mean, score = score)
  }
  lower <- integrate(quadrature$base)
  higher <- integrate(quadrature$fine)
  refinement <- max(abs(lower$score - higher$score) / (1 + abs(higher$score)),
    abs(lower$mean - higher$mean) / (1 + abs(higher$mean)))
  if (!is.finite(refinement) || refinement > 1e-8)
    stop("Quadrature refinement disagreed. Scaled difference = ", format(refinement, digits = 8), ".")
  c(higher, list(refinement = refinement))
}

one <- function(case) {
  model <- case$model
  family <- models[[model]]
  scenario <- case$scenario
  d <- simulate(model, scenario, case$seed)
  split <- fold_rows(d, scenario, case$prediction)
  groups <- unique(as.character(d[split$assessment, "group"]))
  warnings <- character()
  cv <- tryCatch(withCallingHandlers(
    cross_validate_mixed(form(model), d, list(one = split), family = family,
      prediction = case$prediction, keep_fits = TRUE),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = identity)
  fit_ok <- !inherits(cv, "error") && isTRUE(cv$complete)
  failure <- if (inherits(cv, "error")) conditionMessage(cv) else cv$errors$one$message
  stage <- if (inherits(cv, "error")) "call" else cv$errors$one$stage
  if (inherits(cv, "multiScaleR_mixed_cv")) warnings <- c(warnings, unlist(cv$warnings, use.names = FALSE))
  rows <- list()
  refs <- lapply(groups, function(g) {
    training <- d[split$analysis, , drop = FALSE]
    assessment <- d[split$assessment, , drop = FALSE]
    a <- training[as.character(training$group) == g, , drop = FALSE]
    t <- assessment[as.character(assessment$group) == g, , drop = FALSE]
    result <- data.frame(id = case$id, group = g, n_train = nrow(a), n_test = nrow(t),
      reference_status = "not_run_fit_failed", reference_message = NA_character_, observed_score = NA_real_,
      reference_score = NA_real_, absolute_score_error = NA_real_, refinement = NA_real_)
    if (fit_ok) {
      fit <- cv$fits$one
      beta <- glmmTMB::fixef(fit)$cond
      eta <- function(x) beta[1L] + beta[2L] * (x$x - mean(training$x)) /
        stats::sd(training$x) + log(x$exposure)
      oracle <- tryCatch(withCallingHandlers(if (family == "gaussian")
        gaussian_reference(fit, a, t, eta, model == "gaussian_slope") else
          glmm_reference(fit, a, t, eta, family),
        warning = function(w) stop(conditionMessage(w), call. = FALSE)), error = identity)
      result$observed_score <- cv$joint_scores$one$neg_log_score[match(g, cv$joint_scores$one$group)]
      if (inherits(oracle, "error")) {
        result$reference_status <- "failed"
        result$reference_message <- conditionMessage(oracle)
      } else {
        result$reference_status <- "available"
        result$reference_score <- oracle$score
        result$absolute_score_error <- abs(result$observed_score - oracle$score)
        result$refinement <- oracle$refinement
      }
      estimate <- cv$predictions$prediction[match(rownames(t), cv$predictions$row_id)]
      reference <- if (inherits(oracle, "error")) rep(NA_real_, nrow(t)) else oracle$mean
      rows[[g]] <<- data.frame(id = case$id, group = g, row_id = rownames(t),
        estimate = estimate, reference = reference, absolute_error = abs(estimate - reference),
        scaled_error = abs(estimate - reference) / (1 + abs(reference)))
    }
    result
  })
  random_sd <- if (fit_ok) sqrt(diag(glmmTMB::VarCorr(cv$fits$one)$cond$group)) else NA_real_
  fold <- cbind(case, family = family, fit_status = if (fit_ok) "available" else "failed",
    failure_stage = if (length(stage)) stage else NA_character_,
    failure_message = if (length(failure)) failure else NA_character_,
    warnings = if (length(warnings)) paste(unique(warnings), collapse = " | ") else NA_character_,
    n_analysis = length(split$analysis), n_assessment = length(split$assessment),
    fitted_intercept_sd = random_sd[1L], fitted_slope_sd = random_sd[2L])
  list(fold = fold, group = do.call(rbind, refs), row = do.call(rbind, rows))
}

# Materialize the complete cohort plan before fitting. Retain the rows from
# every attempt, including native, scoring and reference failures.
all <- vector("list", nrow(plan))
for (i in seq_len(nrow(plan))) {
  all[[i]] <- one(plan[i, , drop = FALSE])
  if (i %% 14L == 0L) {
    cat("Completed", i, "of", nrow(plan), "folds\n")
    for (item in c("fold", "group", "row"))
      utils::write.csv(do.call(rbind, lapply(all[seq_len(i)], function(x) x[[item]])),
        file.path(output, paste0(item, "-results.csv")), row.names = FALSE)
  }
}
