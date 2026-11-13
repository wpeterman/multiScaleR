args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% c("fit", "assess") || length(args) > 2L)
  stop("Usage: Rscript run.R fit [1|2] OR Rscript run.R assess")
root <- "tools/benchmarks/bayesian_scales/uncorrected_laplace_20261112"
archive <- file.path(root, "results")
dir.create(archive, recursive = TRUE, showWarnings = FALSE)
plan <- data.frame(family = rep(c("gaussian", "gaussian", "nbinom2"), each = 4L),
  design = rep(c("annual", "repeated", "repeated"), each = 4L),
  replicate = rep(seq_len(4L), 3L))
plan$case <- with(plan, sprintf("%s-%s-%02d", family, design, replicate))
plan$seed <- 910000L + seq_len(nrow(plan)) * 20L
methods <- c("nuts", "joint_laplace", "marginal_laplace")
save_record <- function(record, path) {
  temporary <- tempfile(tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(record, temporary, compress = "gzip")
  if (!file.rename(temporary, path)) stop("Cannot publish record ", path)
}
record_path <- function(case, method) file.path(archive, paste(case$case, method, "rds", sep = "."))
source_paths <- c(file.path(root, c("run.R", "simulate.R")), "DESCRIPTION",
  "R/multiScale_bayes.R", "R/bayes_scale_laplace.R", "R/bayes_scale_marginal.R",
  "R/bayes_scale_prediction.R", "R/bayes_scale_diagnostics.R",
  "R/scale_model_prepare.R", "R/scale_response_likelihood.R")
hash <- vapply(source_paths, function(path) digest::digest(file = path, algo = "sha256"), character(1))
config <- list(hash = hash, session = sessionInfo(), plan = plan, methods = methods,
  versions = vapply(c("TMB", "RTMB", "tmbstan", "rstan", "posterior"),
    function(p) as.character(utils::packageVersion(p)), character(1)))
check_hash <- function(result) {
  if (!identical(result$hash, hash)) stop("Existing record uses a different frozen source: ", result$case)
}

if (args[1L] == "fit") {
  pkgload::load_all(compile = FALSE, quiet = TRUE)
  for (p in c("RTMB", "tmbstan", "rstan", "posterior"))
    if (!requireNamespace(p, quietly = TRUE)) stop("Missing study dependency ", p)
  source(file.path(root, "simulate.R"))
  shard <- if (length(args) == 2L) as.integer(args[2L]) else 1L
  if (!shard %in% 1:2) stop("Choose shard 1 or 2.")
  config_file <- file.path(archive, paste0("configuration-", shard, ".rds"))
  if (file.exists(config_file)) check_hash(readRDS(config_file)) else save_record(config, config_file)

  score_fit <- function(fit, input, case, method) {
    spec <- fit$specification
    stopifnot(isTRUE(all.equal(unname(c(spec$center, spec$spread)),
      unname(input$reference), tolerance = 1e-10)))
    derived <- if (method == "nuts") diagnose_bayes_scales(fit) else NULL
    eligible <- isTRUE(fit$diagnostics$passed) && (method != "nuts" || isTRUE(derived$passed))
    sm <- spec$scales$habitat
    term <- sm$random[[1L]]
    rows <- sm$representatives
    log_distance <- sweep(fit$draws$scale_beta[, sm$beta_index, drop = FALSE] %*%
      t(sm$X[rows, , drop = FALSE]), 2L, log(spec$unit_conv), "+") +
      fit$draws$scale_z[, term$z_index[term$index[rows]], drop = FALSE] *
      exp(fit$draws$log_scale_sd[, term$sd_index])
    future <- input$data[input$assessment, , drop = FALSE]
    future_log <- predict_scales(fit, future, new_levels = "sample", scale = "log",
      summary = FALSE, seed = case$seed + 4L)$habitat[, 2L]
    values <- cbind(beta = fit$draws$beta[, 2L], log_distance = log_distance[, 1L],
      scale_sd = exp(fit$draws$log_scale_sd[, 1L]), future_log_distance = future_log)
    truth <- c(input$truth, future_log_distance = input$future_log_distance)
    q <- t(apply(values, 2L, stats::quantile, probs = c(.05, .5, .95), names = FALSE))
    intervals <- data.frame(target = colnames(values), truth = unname(truth[colnames(values)]),
      lower = q[, 1L], median = q[, 2L], upper = q[, 3L], row.names = NULL)
    intervals$cdf_at_truth <- colMeans(sweep(values, 2L, intervals$truth, "<="))
    intervals$covered <- with(intervals, truth >= lower & truth <= upper)
    thresholds <- c(150, 350)
    probs <- if (method == "nuts") diagnose_bayes_probabilities(fit,
      thresholds = list(habitat = thresholds))$probabilities[, c("row", "threshold", "estimate", "mcse")] else
        do.call(rbind, lapply(seq_along(rows), function(j)
          do.call(rbind, lapply(thresholds, function(threshold) {
            z <- log_distance[, j] <= log(threshold)
            estimate <- mean(z)
            data.frame(row = rows[j], threshold, estimate,
              mcse = if (any(z) && !all(z)) sqrt(estimate * (1 - estimate) / length(z)) else NA_real_)
          }))))
    eta <- predict(fit, future, kernel_inputs = input$inputs, row_id = "key", type = "link",
      summary = FALSE, new_scale_levels = "sample", new_response_levels = "sample",
      seed = case$seed + 6L)
    size <- exp(fit$draws$log_dispersion[, 1L])
    lp <- vapply(seq_len(nrow(future)), function(j) if (input$family == "gaussian")
      dnorm(future$y[j], eta[, j], size, log = TRUE) else
      dnbinom(future$y[j], mu = exp(eta[, j]), size = size, log = TRUE), numeric(nrow(eta)))
    joint <- rowSums(lp)
    stopifnot(all(is.finite(values)), !anyNA(joint), any(is.finite(joint)))
    max_joint <- max(joint)
    list(status = if (eligible) "eligible" else "rejected",
      derived_passed = if (is.null(derived)) NA else derived$passed, diagnostics = fit$diagnostics,
      intervals = intervals, probability = probs,
      joint_log_score = -(max_joint + log(mean(exp(joint - max_joint)))))
  }

  fit_case <- function(input, case, method) {
    warnings <- character(); stage <- "fit"; fit <- NULL
    start <- proc.time()[["elapsed"]]
    ans <- tryCatch(withCallingHandlers({
      common <- list(formula = input$formula, inputs = input$inputs,
        data = input$data[input$analysis, , drop = FALSE], scales = input$scales,
        family = input$family, row_id = "key", priors = input$priors)
      control <- switch(method,
        nuts = list(method = "sample", chains = 4L, iter = 3000L, warmup = 1500L,
          seed = case$seed + 1L, control = list(adapt_delta = .99, max_treedepth = 12L)),
        joint_laplace = list(method = "laplace", laplace_draws = 6000L, seed = case$seed + 2L),
        marginal_laplace = list(method = "marginal", laplace_draws = 6000L, seed = case$seed + 3L))
      fit <- do.call(multiScale_bayes, c(common, control))
      stage <- "assessment"
      score_fit(fit, input, case, method)
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) list(status = "failed", stage = stage,
      error = conditionMessage(e), diagnostics = if (is.null(fit)) NULL else fit$diagnostics))
    c(as.list(case), list(method = method, hash = hash,
      elapsed = proc.time()[["elapsed"]] - start, warnings = unique(warnings)), ans)
  }

  for (i in which(plan$replicate %% 2L == shard %% 2L)) {
    case <- plan[i, ]
    input <- simulate_case(case)
    for (method in methods) {
      path <- record_path(case, method)
      if (file.exists(path)) {
        check_hash(readRDS(path))
        next
      }
      message(case$case, " ", method)
      result <- fit_case(input, case, method)
      save_record(result, path)
      message("  ", result$status, if (!is.null(result$error)) paste0(" (", result$error, ")"))
    }
  }
} else {
  get_record <- function(case, method) {
    path <- record_path(case, method)
    if (!file.exists(path)) return(c(as.list(case), list(method = method, status = "unrun")))
    result <- readRDS(path)
    check_hash(result)
    result
  }
  records <- unlist(lapply(seq_len(nrow(plan)), function(i)
    lapply(methods, function(m) get_record(plan[i, ], m))), recursive = FALSE)
  value_or <- function(x, missing) if (is.null(x)) missing else x
  results <- do.call(rbind, lapply(records, function(x) data.frame(
    case = x$case, family = x$family, design = x$design, method = x$method, status = x$status,
    stage = value_or(x$stage, NA_character_), error = value_or(x$error, NA_character_),
    elapsed = value_or(x$elapsed, NA_real_), joint_log_score = value_or(x$joint_log_score, NA_real_),
    warnings = paste(x$warnings, collapse = " | "))))
  write.csv(results, file.path(archive, "fits.csv"), row.names = FALSE)
  compact <- function(item) do.call(rbind, lapply(records, function(x)
    if (!is.null(x[[item]])) cbind(results[match(paste(x$case, x$method),
      paste(results$case, results$method)), c("case", "family", "design", "method", "status")],
      x[[item]])))
  intervals <- compact("intervals"); probabilities <- compact("probability")
  if (nrow(as.data.frame(intervals))) write.csv(intervals, file.path(archive, "intervals.csv"), row.names = FALSE)
  if (nrow(as.data.frame(probabilities))) write.csv(probabilities, file.path(archive, "probabilities.csv"), row.names = FALSE)
  paired <- function(tab, keys) {
    if (is.null(tab)) return(NULL)
    a <- tab[tab$method == "nuts", , drop = FALSE]
    b <- tab[tab$method != "nuts", , drop = FALSE]
    merge(b, a, by = c("case", "family", "design", keys), suffixes = c("", "_nuts"))
  }
  p <- paired(probabilities, c("row", "threshold"))
  if (nrow(as.data.frame(p))) {
    p$absolute_difference <- abs(p$estimate - p$estimate_nuts)
    p$combined_mcse <- sqrt(p$mcse^2 + p$mcse_nuts^2)
    p$material <- p$absolute_difference > pmax(.05, 5 * p$combined_mcse)
    write.csv(p, file.path(archive, "paired_probabilities.csv"), row.names = FALSE)
  }
  q <- paired(intervals, "target")
  if (nrow(as.data.frame(q))) {
    q$width_ratio <- (q$upper - q$lower) / (q$upper_nuts - q$lower_nuts)
    q$median_error <- q$median - q$median_nuts
    write.csv(q, file.path(archive, "paired_intervals.csv"), row.names = FALSE)
  }
  print(with(results, table(family, design, method, status)))
}
