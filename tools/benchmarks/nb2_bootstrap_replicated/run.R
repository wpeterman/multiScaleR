# Run from the radius worktree. Supply the corrected package source and output.
# Modes are separate studies and require separate output directories.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 3L, args[3L] %in% c("bounded", "precision"))
source_root <- normalizePath(args[1L], mustWork = TRUE)
output <- normalizePath(args[2L], mustWork = FALSE)
mode <- args[3L]
n_each <- if (mode == "bounded") 8L else 48L
B <- if (mode == "bounded") 49L else 99L
level <- .90
scenarios <- c("interior", "weak", "near_bound")
plan <- expand.grid(scenario = scenarios, replicate = seq_len(n_each),
  KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
plan$seed <- (if (mode == "bounded") 1500000L else 1600000L) +
  seq_len(nrow(plan)) * 100L
plan$bootstrap_seed <- plan$seed + 1L
plan$id <- sprintf("%s-%02d", plan$scenario, plan$replicate)

capture <- function(expr) {
  messages <- character()
  value <- tryCatch(withCallingHandlers(expr, warning = function(w) {
    messages <<- c(messages, conditionMessage(w))
    invokeRestart("muffleWarning")
  }), error = identity)
  list(value = value, warnings = unique(messages))
}

atomic <- function(x, file) {
  tmp <- tempfile(tmpdir = dirname(file))
  on.exit(unlink(tmp), add = TRUE)
  saveRDS(x, tmp)
  if (!file.rename(tmp, file)) stop("Cannot publish ", file)
}

script <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1L])
stopifnot(file.exists(script), requireNamespace("digest", quietly = TRUE),
  requireNamespace("MASS", quietly = TRUE))
package_files <- c("DESCRIPTION", "NAMESPACE",
  file.path("R", list.files(file.path(source_root, "R"), pattern = "[.]R$")),
  file.path("src", list.files(file.path(source_root, "src"),
    pattern = "[.](c|cpp|h|dll|so)$")))
hash_files <- function() {
  files <- file.path(source_root, package_files)
  stats::setNames(vapply(files, digest::digest, "", file = TRUE,
    algo = "sha256"), package_files)
}

# This is the pilot's independent, fixed generator with new dataset seeds.
simulate_case <- function(case) {
  set.seed(case$seed)
  n <- 110L
  sigma <- if (case$scenario == "near_bound") .73 else .24
  effect <- if (case$scenario == "weak") .12 else .95
  size <- if (case$scenario == "near_bound") Inf else 2
  d <- seq(.015, 1, length.out = 21L)
  count <- matrix(sample(1:5, n * length(d), replace = TRUE), n)
  radial <- outer(rnorm(n), sin(7 * d)) + outer(rnorm(n), cos(12 * d)) +
    outer(rnorm(n), d) + matrix(rnorm(n * length(d), sd = .25), n)
  V <- radial * count
  habitat <- function(scale) {
    w <- exp(-(d^2 - min(d)^2) / (2 * scale^2))
    as.numeric(base::scale(drop(V %*% w) / drop(count %*% w)))
  }
  data <- data.frame(hab = habitat(.23), z = rnorm(n), effort = runif(n, .8, 1.2),
    site = factor(rep(seq_len(20L), each = 5L, length.out = n)))
  mu <- exp(.4 + effect * habitat(sigma) + .3 * data$z + log(data$effort))
  data$y <- if (is.infinite(size)) rpois(n, mu) else rnbinom(n, size = size, mu = mu)
  inputs <- list(raw_cov = NULL, d_list = NULL, min_D = 50, max_D = 800,
    unit_conv = 1000, kernel = "gaussian",
    binned = list(bin_rep = d, vsum = list(hab = V), csum = list(hab = count),
      point_ids = rownames(data), n_points = n))
  list(data = data, inputs = inputs, sigma = sigma * 1000, effect = effect, size = size)
}

ci95 <- function(k, n) if (n == 0L) rep(NA_real_, 2L) else
  as.numeric(stats::binom.test(k, n)$conf.int)

fit_case <- function(case) {
  clock <- proc.time()[["elapsed"]]
  f <- simulate_case(case)
  summary <- data.frame(case, true_sigma = f$sigma, effect = f$effect, size = f$size,
    elapsed_seconds = NA_real_, status = "initial_fit_failed", error = "",
    fitted_sigma = NA_real_, fitted_poisson = NA, profile_accept = NA,
    point_lr = NA_real_, point_p = NA_real_, mc_lower = NA_real_, mc_upper = NA_real_,
    nominal_accept = NA, mc_ambiguous = NA, exceedances = NA_integer_,
    paired_attempted = 0L, paired_failed = 0L, paired_warned = 0L,
    p_at_cutoff = NA, profile_status = NA_character_, bootstrap_status = NA_character_,
    initial_warnings = "", fit_warnings = "", inference_warnings = "")
  audit <- list(fit_starts = NULL, inference = NULL)
  done <- function() {
    summary$elapsed_seconds <- proc.time()[["elapsed"]] - clock
    list(summary = summary, audit = audit)
  }
  initial <- capture(stats::glm(y ~ hab + z + offset(log(effort)), data = f$data,
    family = MASS::negative.binomial(2)))
  summary$initial_warnings <- paste(initial$warnings, collapse = " | ")
  if (inherits(initial$value, "error") || !isTRUE(initial$value$converged)) {
    summary$error <- if (inherits(initial$value, "error"))
      conditionMessage(initial$value) else "The initial fixed-size GLM did not converge."
    return(done())
  }
  summary$status <- "point_fit_failed"
  opt <- capture(multiScale_optim(initial$value, f$inputs, backend = "joint_nb2",
    n_cores = 1, verbose = FALSE,
    profile_control = list(starts = matrix(c(.075, .20, .43, .71), ncol = 1L))))
  summary$fit_warnings <- paste(opt$warnings, collapse = " | ")
  if (inherits(opt$value, "error")) {
    summary$error <- conditionMessage(opt$value)
    audit$fit_starts <- opt$value$start_results
    return(done())
  }
  fit <- opt$value
  audit$fit_starts <- fit$profile_fit$start_results
  summary$fitted_sigma <- fit$scale_est$Mean[1L]
  summary$fitted_poisson <- fit$dispersion_est$exact_poisson
  summary$status <- "profile_or_bootstrap_failed"
  result <- capture(scale_inference(fit, sigma_values = f$sigma,
    level = level, bootstrap = B, bootstrap_values = f$sigma,
    seed = case$bootstrap_seed, n_starts = 2))
  summary$inference_warnings <- paste(result$warnings, collapse = " | ")
  if (inherits(result$value, "error")) {
    summary$error <- conditionMessage(result$value)
    return(done())
  }
  result <- result$value
  audit$inference <- result[c("grid", "reference", "profile_audit",
    "diagnostics", "bootstrap_records", "bootstrap_set")]
  truth <- result$grid[result$grid$sigma == f$sigma, , drop = FALSE]
  pairs <- result$bootstrap_records
  stopifnot(nrow(truth) == 1L, nrow(pairs) == B,
    all(pairs$generating_sigma == f$sigma), is.null(result$bootstrap_set))
  summary$paired_attempted <- nrow(pairs)
  summary$paired_failed <- sum(!pairs$ok)
  summary$paired_warned <- sum(nzchar(pairs$warnings))
  summary$profile_status <- result$diagnostics$profile_status
  summary$bootstrap_status <- result$diagnostics$bootstrap_status
  summary$profile_accept <- truth$profile_accept
  summary$point_lr <- truth$lr
  summary$point_p <- truth$bootstrap_p
  summary$mc_lower <- truth$mc_lower
  summary$mc_upper <- truth$mc_upper
  summary$nominal_accept <- truth$bootstrap_accept
  summary$mc_ambiguous <- truth$mc_ambiguous
  if (summary$paired_failed) {
    summary$status <- "paired_refit_failed"
    summary$error <- paste(unique(pairs$error[!pairs$ok]), collapse = " | ")
    stopifnot(is.na(truth$bootstrap_p), is.na(truth$bootstrap_accept))
    return(done())
  }

  # Independently recover the tail, adjusted Monte Carlo limits, and exact
  # decision at the 10 percent cutoff without using a floating complement.
  exceed <- sum(pairs$lr >= truth$lr - 1e-8)
  tail <- (1 + exceed) / (B + 1)
  mc <- (1 + B * ci95(exceed, B)) / (B + 1)
  accept <- 10L * (1L + exceed) > B + 1L
  ambiguous <- mc[1L] <= .1 && mc[2L] > .1
  stopifnot(abs(truth$bootstrap_p - tail) <= 1e-14,
    max(abs(c(truth$mc_lower, truth$mc_upper) - mc)) <= 1e-14,
    identical(truth$bootstrap_accept, accept), identical(truth$mc_ambiguous, ambiguous))
  summary$status <- if (accept) "accepted" else "rejected"
  summary$exceedances <- exceed
  summary$p_at_cutoff <- 10L * (1L + exceed) == B + 1L
  done()
}

summarize <- function(cases) {
  do.call(rbind, lapply(scenarios, function(scenario) {
    x <- cases[cases$scenario == scenario & cases$status != "pending", , drop = FALSE]
    n <- nrow(x)
    valid <- !is.na(x$nominal_accept)
    accepted <- sum(x$nominal_accept[valid])
    clear <- valid & !is.na(x$mc_ambiguous) & !x$mc_ambiguous
    certain_accept <- sum(x$nominal_accept[clear])
    ambiguous <- sum(x$mc_ambiguous[valid], na.rm = TRUE)
    nominal_ci <- ci95(accepted, sum(valid))
    operational_ci <- ci95(accepted, n)
    data.frame(scenario, planned = n_each, attempted = n, available = sum(valid),
      accepted, rejected = sum(valid) - accepted, unavailable = n - sum(valid),
      mc_ambiguous = ambiguous, p_at_cutoff = sum(x$p_at_cutoff, na.rm = TRUE),
      nominal_rate = if (sum(valid)) accepted / sum(valid) else NA_real_,
      nominal_lower = nominal_ci[1L], nominal_upper = nominal_ci[2L],
      returned_acceptance_rate = if (n) accepted / n else NA_real_,
      returned_acceptance_lower = operational_ci[1L], returned_acceptance_upper = operational_ci[2L],
      certain_accept, certain_reject = sum(clear) - certain_accept,
      sensitivity_min = if (n) certain_accept / n else NA_real_,
      sensitivity_max = if (n) 1 - (sum(clear) - certain_accept) / n else NA_real_,
      paired_attempted = sum(x$paired_attempted), paired_failed = sum(x$paired_failed),
      seconds_total = sum(x$elapsed_seconds),
      seconds_median = if (n) stats::median(x$elapsed_seconds) else NA_real_)
  }))
}

main <- function() {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  case_dir <- file.path(output, "cases")
  dir.create(case_dir, showWarnings = FALSE)
  lock <- file.path(output, "run.lock")
  if (!dir.create(lock, showWarnings = FALSE)) stop("The output directory is already in use.")
  stopifnot(identical(dirname(normalizePath(lock, winslash = "/")),
    normalizePath(output, winslash = "/")))
  on.exit(unlink(lock, recursive = TRUE), add = TRUE)
  inputs <- hash_files()
  corrected <- readLines(file.path(source_root, "R", "scale_inference.R"), warn = FALSE)
  stopifnot(any(grepl('grid$bootstrap_accept[i] <- 1 - tail["p"] < level', corrected, fixed = TRUE)))
  pkgload::load_all(path = source_root, compile = FALSE, quiet = TRUE)
  stopifnot(identical(inputs, hash_files()))
  config <- list(mode = mode, n_each = n_each, B = B, level = level, plan = plan,
    R = R.version.string, package_files = inputs,
    script_hash = digest::digest(script, file = TRUE, algo = "sha256"))
  config_path <- file.path(output, "configuration.rds")
  if (file.exists(config_path)) stopifnot(identical(readRDS(config_path), config)) else {
    atomic(config, config_path)
    utils::write.csv(plan, file.path(output, "plan.csv"), row.names = FALSE)
    writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"))
  }
  paths <- file.path(case_dir, paste0(plan$id, ".rds"))
  for (i in seq_len(nrow(plan))) {
    if (file.exists(paths[i])) next
    if (file.exists(file.path(output, "STOP"))) break
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "starting", plan$id[i], "\n")
    flush.console()
    result <- fit_case(plan[i, , drop = FALSE])
    atomic(result, paths[i])
    z <- result$summary
    cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "finished", z$id,
      z$status, "p =", z$point_p, "ambiguous =", z$mc_ambiguous,
      "paired failures =", z$paired_failed, "seconds =", z$elapsed_seconds, "\n")
    flush.console()
  }
  stopifnot(identical(inputs, hash_files()))
  outcomes <- lapply(paths[file.exists(paths)], readRDS)
  complete <- do.call(rbind, lapply(outcomes, "[[", "summary"))
  if (is.null(complete)) stop("No case completed.")
  observed <- complete[match(plan$id, complete$id), setdiff(names(complete), names(plan))]
  cases <- cbind(plan, observed, row.names = NULL)
  cases$status[is.na(cases$status)] <- "pending"
  pairs <- do.call(rbind, lapply(outcomes, function(x) {
    rows <- x$audit$inference$bootstrap_records
    if (!is.null(rows) && nrow(rows)) cbind(id = x$summary$id, rows)
  }))
  result <- summarize(cases)
  utils::write.csv(cases, file.path(output, "datasets.csv"), row.names = FALSE)
  utils::write.csv(pairs, file.path(output, "pairs.csv"), row.names = FALSE)
  utils::write.csv(result, file.path(output, "scenario_summary.csv"), row.names = FALSE)
  utils::write.csv(data.frame(path = names(inputs), sha256 = unname(inputs)),
    file.path(output, "source_sha256.csv"), row.names = FALSE)
  print(result, row.names = FALSE)
}

main()
