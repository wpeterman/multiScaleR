# Reproducible numerical qualification for Bayesian scale approximations.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) %in% 2:3, args[1] %in% c("pilot", "qualification"))
mode <- args[1L]
output <- args[2L]
dir.create(output, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(output))
pkgload::load_all(quiet = TRUE)
for (p in c("RTMB", "tmbstan", "rstan", "posterior"))
  if (!requireNamespace(p, quietly = TRUE)) stop("Missing study dependency ", p)
ns <- asNamespace("multiScaleR")
family <- c("gaussian", "poisson", "binomial", "nbinom2")
design <- c("annual", "repeated")
plan <- expand.grid(family = family, design = design,
  replicate = seq_len(if (mode == "pilot") 1L else 48L), stringsAsFactors = FALSE)
plan$case <- sprintf("%s-%s-%02d", plan$family, plan$design, plan$replicate)
plan$seed <- (if (mode == "pilot") 710000L else 810000L) + seq_len(nrow(plan)) * 20L
selector <- if (length(args) == 3L) args[3L] else "all"
parts <- strsplit(selector, ":", fixed = TRUE)[[1L]]
selected <- if (selector == "assess") integer() else if (selector == "all")
  seq_len(nrow(plan)) else which(startsWith(plan$case, parts[1L]))
if (!length(selected) && selector != "assess") stop("Unknown case or panel ", selector)
if (length(parts) > 1L) {
  if (length(parts) != 2L || !parts[2L] %in% c("odd", "even"))
    stop("Use :odd or :even to split independent replicates.")
  parity <- if (parts[2L] == "odd") 1L else 0L
  selected <- selected[plan$replicate[selected] %% 2L == parity]
}
methods <- c("nuts", "joint_importance", "conditional_beta")
source <- sort(c(list.files("R", full.names = TRUE), list.files("src", full.names = TRUE,
  pattern = "[.](c|cpp|h)$"), "DESCRIPTION", "NAMESPACE",
  "tools/benchmarks/bayesian_scales/qualification_20261112/run.R"))
fingerprint <- digest::digest(lapply(source, digest::digest, file = TRUE, algo = "sha256"))
config <- list(mode = mode, source = fingerprint, R = R.version.string,
  packages = vapply(c("RTMB", "tmbstan", "rstan", "posterior", "TMB"),
    function(p) as.character(packageVersion(p)), character(1)))
atomic <- function(x, file) {
  tmp <- tempfile(tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp, compress = "gzip")
  if (!file.rename(tmp, file)) stop("Could not publish record ", file)
}
metadata <- file.path(output, "configuration.rds")
if (file.exists(metadata)) {
  if (!identical(readRDS(metadata), config)) stop("Study configuration has changed. Choose a new output directory.")
} else atomic(config, metadata)

simulate_case <- function(case) {
  set.seed(case$seed)
  repeated <- case$design == "repeated"
  kernel <- if (repeated) "exp" else "gaussian"
  n_year <- if (repeated) 4L else 3L
  sites <- if (repeated) 9L else 6L
  train <- expand.grid(site = seq_len(sites), year = seq_len(n_year))
  if (repeated) train <- train[(train$site + train$year) %% 5L != 0, , drop = FALSE]
  future <- data.frame(site = c(1L, sites + 1L, sites + 1L),
    year = c(n_year, n_year + 1L, n_year + 1L))
  data <- rbind(train, future)
  data$time <- (data$year - 2) / 2
  data$key <- sprintf("row%02d", seq_len(nrow(data)))
  data$offset <- runif(nrow(data), -.2, .2)
  data$trials <- sample(4:9, nrow(data), replace = TRUE)
  d <- seq(.025, 1, length.out = 12L)
  n <- nrow(data)
  count <- matrix(sample(2:6, n * length(d), replace = TRUE), n)
  habitat <- outer(rnorm(n), sin(8 * d)) + outer(rnorm(n), cos(12 * d)) +
    matrix(rnorm(n * length(d), sd = .25), n)
  powers <- if (kernel == "gaussian") 2 else 1
  means <- function(log_scale) {
    w <- exp(-outer(exp(-powers * log_scale), d^powers - min(d)^powers) / powers)
    rowSums(w * habitat * count) / rowSums(w * count)
  }
  analysis <- seq_len(nrow(train))
  assessment <- setdiff(seq_len(n), analysis)
  p <- list(scale_median = 250, scale_intercept_sd = .45, scale_coef_sd = .25,
    scale_group_sd = .30, beta_sd = .70, response_group_sd = .35,
    dispersion_log_mean = if (case$family == "gaussian") log(.6) else log(2),
    dispersion_log_sd = .35)
  base <- rnorm(1, log(p$scale_median / 1000), p$scale_intercept_sd)
  trend <- if (repeated) rnorm(1, sd = p$scale_coef_sd) else 0
  scale_sd <- abs(rnorm(1, sd = p$scale_group_sd))
  log_scale <- base + trend * data$time + rnorm(n_year + 1L, sd = scale_sd)[data$year]
  beta <- rnorm(2, sd = p$beta_sd)
  response_sd <- if (repeated) abs(rnorm(1, sd = p$response_group_sd)) else 0
  re <- rnorm(sites + 1L, sd = response_sd)[data$site]
  dispersion <- exp(rnorm(1, p$dispersion_log_mean, p$dispersion_log_sd))
  ref <- means(rep(log(p$scale_median / 1000), n))[analysis]
  eta <- beta[1L] + beta[2L] * (means(log_scale) - mean(ref)) / sd(ref) + re + data$offset
  data$y <- switch(case$family, gaussian = rnorm(n, eta, dispersion),
    poisson = rpois(n, exp(eta)), binomial = rbinom(n, data$trials, plogis(eta)),
    nbinom2 = rnbinom(n, mu = exp(eta), size = dispersion))
  data$miss <- data$trials - data$y
  data$year <- factor(data$year)
  data$site <- factor(data$site)
  response <- if (case$family == "binomial") "cbind(y, miss)" else "y"
  formula <- as.formula(paste(response, "~ habitat + offset(offset)",
    if (repeated) "+ (1 | site)" else ""))
  scale <- list(habitat = if (repeated) ~ 1 + time + (1 | year) else ~ 1 + (1 | year))
  # Separate prepared inputs from training responses. Bins have fixed geometry.
  inputs <- list(kernel = kernel, min_D = 25, max_D = 1000, unit_conv = 1000,
    binned = list(bin_rep = d, vsum = list(habitat = habitat * count),
      csum = list(habitat = count), point_ids = data$key, n_points = n))
  list(data = data, analysis = analysis, assessment = assessment,
    formula = formula, scales = scale, inputs = inputs, priors = p, family = case$family,
    reference = c(mean = mean(ref), sd = sd(ref)),
    truth = c(beta = beta[2L], log_distance = log_scale[1L] + log(1000), scale_sd = scale_sd))
}

score_fit <- function(fit, input, case, importance) {
  spec <- fit$specification
  stopifnot(isTRUE(all.equal(unname(c(spec$center, spec$spread)),
    unname(input$reference), tolerance = 1e-10)))
  derived <- if (importance) NULL else diagnose_bayes_scales(fit)
  eligible <- isTRUE(fit$diagnostics$passed) && (importance || isTRUE(derived$passed))
  if (importance) {
    weight <- fit$importance$weights
    fit$draws <- ns$.msr_bayes_parameter_blocks(fit$importance$candidates, fit$parameters)
  } else weight <- rep(1 / nrow(fit$draws$beta), nrow(fit$draws$beta))
  # Use coherent parameter vectors for the three future responses.
  future <- input$data[input$assessment, , drop = FALSE]
  eta <- predict(fit, newdata = future, kernel_inputs = input$inputs, row_id = "key",
    type = "link", summary = FALSE, new_scale_levels = "sample", new_response_levels = "sample",
    seed = case$seed + 6L)
  observed <- rep(future$y, each = length(weight))
  trials <- rep(future$trials, each = length(weight))
  dispersion <- if (ncol(fit$draws$log_dispersion))
    rep(exp(fit$draws$log_dispersion[, 1L]), nrow(future)) else 1
  ll <- switch(input$family, gaussian = dnorm(observed, eta, dispersion, log = TRUE),
    poisson = dpois(observed, exp(eta), log = TRUE),
    binomial = dbinom(observed, trials, plogis(eta), log = TRUE),
    nbinom2 = dnbinom(observed, mu = exp(eta), size = dispersion, log = TRUE))
  joint <- rowSums(matrix(ll, nrow = length(weight))) + log(weight)
  stopifnot(!anyNA(joint), any(is.finite(joint)))
  joint_log_score <- -(max(joint) + log(sum(exp(joint - max(joint)))))
  sm <- spec$scales$habitat
  term <- sm$random[[1L]]
  rows <- sm$representatives
  log_distance <- sweep(fit$draws$scale_beta[, sm$beta_index, drop = FALSE] %*%
    t(sm$X[rows, , drop = FALSE]), 2L, log(spec$unit_conv), "+") +
    fit$draws$scale_z[, term$z_index[term$index[rows]], drop = FALSE] *
    exp(fit$draws$log_scale_sd[, term$sd_index])
  rank <- c(beta = sum(weight * (fit$draws$beta[, 2L] <= input$truth[["beta"]])),
    log_distance = sum(weight * (log_distance[, 1L] <= input$truth[["log_distance"]])),
    scale_sd = sum(weight * (exp(fit$draws$log_scale_sd[, 1L]) <= input$truth[["scale_sd"]])))
  thresholds <- c(150, 350)
  if (importance) {
    probability <- do.call(rbind, lapply(seq_along(rows), function(j)
      do.call(rbind, lapply(thresholds, function(threshold) {
        indicator <- log_distance[, j] <= log(threshold)
        estimate <- sum(weight * indicator)
        mcse <- if (any(indicator) && !all(indicator))
          sqrt(sum(weight^2 * (indicator - estimate)^2)) else NA_real_
        data.frame(row = rows[j], threshold, estimate, mcse)
      }))))
  } else probability <- diagnose_bayes_probabilities(fit,
    thresholds = list(habitat = thresholds))$probabilities[, c("row", "threshold", "estimate", "mcse")]
  list(status = if (eligible) "eligible" else "rejected", diagnostics = fit$diagnostics,
    derived = if (importance) NULL else derived$parameters, rank = rank,
    probability = probability, joint_log_score = joint_log_score)
}

fit_case <- function(input, case, method) {
  start <- proc.time()[["elapsed"]]
  warnings <- character()
  stage <- "fit"
  fit <- NULL
  ans <- tryCatch(withCallingHandlers({
    training <- input$data[input$analysis, , drop = FALSE]
    fit <- if (method == "nuts") multiScale_bayes(input$formula, input$inputs, training,
      scales = input$scales, family = input$family, row_id = "key", priors = input$priors,
      method = "sample", chains = 4L, iter = 3000L, warmup = 1500L, seed = case$seed + 1L,
      control = list(adapt_delta = .99, max_treedepth = 12L))
    else multiScale_bayes(input$formula, input$inputs, training,
      scales = input$scales, family = input$family, row_id = "key", priors = input$priors,
      method = "marginal", seed = case$seed + 2L, laplace_draws = 256L,
      control = list(correction = "importance", conditional_beta = method == "conditional_beta",
        importance_min = 32768L, importance_max = 65536L, importance_df = 5,
        importance_refinements = 1L))
    stage <- "assessment"
    score_fit(fit, input, case, method != "nuts")
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = function(e) list(status = "failed", stage = stage, error = conditionMessage(e),
    diagnostics = if (is.null(fit)) NULL else fit$diagnostics))
  c(list(case = case$case, family = case$family, design = case$design, method = method,
    elapsed = proc.time()[["elapsed"]] - start, warnings = unique(warnings)), ans)
}

for (i in selected) {
  case <- plan[i, ]
  input <- tryCatch(simulate_case(case), error = identity)
  for (method in methods) {
    target <- file.path(output, paste(case$case, method, "rds", sep = "."))
    if (file.exists(target)) next
    cat(case$case, method, "\n")
    result <- if (!inherits(input, "error")) fit_case(input, case, method) else
      list(case = case$case, family = case$family, design = case$design, method = method,
        status = "failed", error = conditionMessage(input))
    atomic(result, target)
    cat("  ", result$status, if (is.null(result$error)) "" else result$error, "\n")
  }
}

read_case <- function(i, method) {
  case <- plan[i, ]
  file <- file.path(output, paste(case$case, method, "rds", sep = "."))
  if (file.exists(file)) readRDS(file) else
    list(case = case$case, family = case$family, design = case$design, method = method,
      status = "unrun")
}
results <- unlist(lapply(seq_len(nrow(plan)), function(i)
  lapply(methods, function(method) read_case(i, method))), recursive = FALSE)
summary <- do.call(rbind, lapply(results, function(x)
  data.frame(case = x$case, family = x$family, design = x$design, method = x$method,
    status = x$status, elapsed = if (is.null(x$elapsed)) NA_real_ else x$elapsed,
    joint_log_score = if (is.null(x$joint_log_score)) NA_real_ else x$joint_log_score,
    error = if (is.null(x$error)) NA_character_ else x$error)))
tag <- if (!selector %in% c("all", "assess")) paste0("-", sub(":", "-", selector)) else ""
write.csv(summary, file.path(output, paste0("fits", tag, ".csv")), row.names = FALSE)
cells <- split(summary, interaction(summary$family, summary$design, summary$method, drop = TRUE))
eligibility <- do.call(rbind, lapply(cells, function(x) {
  n <- nrow(x)
  eligible <- sum(x$status == "eligible")
  unrun <- sum(x$status == "unrun")
  lower <- if (eligible) qbeta(.05, eligible, n - eligible + 1L) else 0
  data.frame(family = x$family[1L], design = x$design[1L], method = x$method[1L],
    planned = n, eligible = eligible, rejected = sum(x$status == "rejected"),
    failed = sum(x$status == "failed"), unrun = unrun, lower_95 = lower,
    availability_qualified = mode == "qualification" && unrun == 0L && lower >= .90)
}))
write.csv(eligibility, file.path(output, paste0("eligibility", tag, ".csv")), row.names = FALSE)
ranks <- do.call(rbind, lapply(results, function(x) if (length(x$rank))
  data.frame(case = x$case, family = x$family, design = x$design, method = x$method,
    status = x$status, quantity = names(x$rank), cdf_at_truth = unname(x$rank))))
if (nrow(as.data.frame(ranks))) {
  write.csv(ranks, file.path(output, paste0("ranks", tag, ".csv")), row.names = FALSE)
  selected <- ranks[ranks$status == "eligible", , drop = FALSE]
  cells <- split(selected, interaction(selected$family, selected$design, selected$method,
    selected$quantity, drop = TRUE))
  if (length(cells)) {
    coverage <- do.call(rbind, lapply(cells, function(x) {
      count <- sum(x$cdf_at_truth >= .05 & x$cdf_at_truth <= .95)
      limits <- binom.test(count, nrow(x))$conf.int
      cbind(x[1L, c("family", "design", "method", "quantity")], n = nrow(x),
        coverage_90 = count / nrow(x), lower_95 = limits[1L], upper_95 = limits[2L])
    }))
    write.csv(coverage, file.path(output, paste0("coverage", tag, ".csv")), row.names = FALSE)
  }
}
comparisons <- list()
for (i in seq_len(nrow(plan))) {
  ref <- read_case(i, "nuts")
  for (method in methods[-1L]) {
    alt <- read_case(i, method)
    if (is.null(ref$probability) || is.null(alt$probability)) next
    cmp <- merge(ref$probability, alt$probability, by = c("row", "threshold"),
      suffixes = c("_nuts", "_importance"))
    cmp$absolute_difference <- abs(cmp$estimate_nuts - cmp$estimate_importance)
    cmp$combined_mcse <- sqrt(cmp$mcse_nuts^2 + cmp$mcse_importance^2)
    cmp$material_difference <- cmp$absolute_difference > pmax(.05, 5 * cmp$combined_mcse)
    comparisons[[length(comparisons) + 1L]] <- cbind(case = ref$case, method = method,
      reference_status = ref$status, importance_status = alt$status, cmp)
  }
}
if (length(comparisons))
  write.csv(do.call(rbind, comparisons), file.path(output, paste0("probabilities", tag, ".csv")),
    row.names = FALSE)
print(with(summary, table(family, method, status)))
