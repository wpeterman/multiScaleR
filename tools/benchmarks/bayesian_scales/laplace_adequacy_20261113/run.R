# Each shard fits independent datasets serially. Assessment is read-only
# with respect to the immutable, compact fit and audit records.
args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1L] %in% c("fit", "assess") || length(args) > 2L)
  stop("Usage: Rscript run.R fit 1|2 OR Rscript run.R assess")
root <- "tools/benchmarks/bayesian_scales/laplace_adequacy_20261113"
archive <- file.path(root, "records")
dir.create(archive, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "design.R"))
methods <- c("nuts", "laplace", "marginal")
budget <- 65536L
paths <- c(file.path(root, c("design.R", "run.R")), "DESCRIPTION",
  sort(list.files("R", full.names = TRUE, pattern = "[.]R$")),
  sort(list.files("src", full.names = TRUE, pattern = "[.](c|cpp|h)$")))
hash <- vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1))
code_hash <- digest::digest(hash)
configuration <- list(hash = hash, session = R.version.string, locale = Sys.getlocale(),
  versions = vapply(c("TMB", "RTMB", "tmbstan", "rstan", "posterior"),
    function(p) as.character(utils::packageVersion(p)), character(1)),
  plan = plan, methods = methods, budget = budget)
save_once <- function(x, file) {
  temp <- tempfile(tmpdir = dirname(file))
  on.exit(unlink(temp), add = TRUE)
  saveRDS(x, temp, compress = "gzip")
  if (!file.rename(temp, file)) stop("Cannot publish ", file)
}
record_file <- function(case, method) file.path(archive, paste(case$case, method, "rds", sep = "."))
if (args[1L] == "fit") {
  pkgload::load_all(compile = FALSE, quiet = TRUE)
  shard <- if (length(args) == 2L) as.integer(args[2L]) else NA_integer_
  if (!shard %in% 1:2) stop("Choose fit shard 1 or 2.")
  config_file <- file.path(archive, paste0("configuration-", shard, ".rds"))
  if (file.exists(config_file)) {
    if (!identical(readRDS(config_file), configuration)) stop("Study configuration changed.")
  } else save_once(configuration, config_file)

  fit_case <- function(case, method, input) {
    warnings <- character(); stage <- "simulation"; fit <- NULL; diagnostics <- NULL
    start <- proc.time()[["elapsed"]]
    ans <- tryCatch(withCallingHandlers({
      if (inherits(input, "error")) stop(input)
      stage <- "fit"
      control <- if (method == "nuts") list(method = "sample", chains = 4L, iter = 4000L,
        warmup = 2000L, seed = case$seed + 1L,
        control = list(adapt_delta = .99, max_treedepth = 12L)) else
          list(method = method, laplace_draws = 8000L,
            seed = case$seed + match(method, methods))
      fit <- do.call(multiScale_bayes, c(list(formula = input$formula, inputs = input$inputs,
        data = input$data, scales = input$scales, family = case$family, priors = input$priors,
        row_id = "key"), control))
      diagnostics <- fit$diagnostics
      stage <- "assessment"
      spec <- fit$specification
      stopifnot(isTRUE(all.equal(unname(c(spec$center, spec$spread)), unname(input$reference),
        tolerance = 1e-10)))
      thresholds <- list(habitat = c(150, 250, 350))
      if (method == "nuts") {
        derived <- diagnose_bayes_scales(fit)
        probability <- diagnose_bayes_probabilities(fit, thresholds)$probabilities
        status <- if (isTRUE(derived$passed)) "eligible" else "rejected"
        list(status = status, raw = diagnostics, derived = derived$parameters,
          probability = probability)
      } else {
        audit <- diagnose_bayes_laplace(fit, thresholds, candidates = budget,
          seed = case$seed + 9L)
        status <- if (isTRUE(diagnostics$passed)) "eligible" else "rejected"
        list(status = status, raw = diagnostics, audit = audit)
      }
    }, warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = function(e) list(status = "failed", stage = stage,
      error = conditionMessage(e), raw = diagnostics, attempts = e$starts))
    c(as.list(case[c("case", "scenario", "family", "replicate", "seed")]),
      list(method = method, hash = code_hash, elapsed = proc.time()[["elapsed"]] - start,
        warnings = unique(warnings)), ans)
  }

  for (i in which(seq_len(nrow(plan)) %% 2L == shard %% 2L)) {
    case <- plan[i, ]
    input <- tryCatch(simulate_case(case), error = identity)
    for (method in methods) {
      if (method == "marginal" && case$years == 1L) next
      file <- record_file(case, method)
      if (file.exists(file)) {
        if (!identical(readRDS(file)$hash, code_hash)) stop("Existing fit hash mismatch.")
        next
      }
      message(case$case, " ", method)
      result <- fit_case(case, method, input)
      save_once(result, file)
      message("  ", result$status, if (!is.null(result$error)) paste0(" (", result$error, ")"))
    }
  }
} else {
  read_record <- function(case, method) {
    file <- record_file(case, method)
    if (!file.exists(file)) return(c(as.list(case[c("case", "scenario", "family", "replicate", "seed")]),
      list(method = method, status = if (method == "marginal" && case$years == 1L) "outside_scope" else "unrun")))
    result <- readRDS(file)
    if (!identical(result$hash, code_hash)) stop("Fit hash mismatch.")
    result
  }
  fits <- unlist(lapply(seq_len(nrow(plan)), function(i)
    lapply(methods, function(method) read_record(plan[i, ], method))), recursive = FALSE)
  value_or <- function(x, y) if (is.null(x)) y else x
  index <- do.call(rbind, lapply(fits, function(x) data.frame(case = x$case, scenario = x$scenario,
    family = x$family, replicate = x$replicate, seed = x$seed, method = x$method,
    status = x$status, error = value_or(x$error, NA_character_),
    warnings = paste(x$warnings, collapse = " | "), elapsed = value_or(x$elapsed, NA_real_),
    audit_available = value_or(x$audit$available, NA),
    audit_passed = value_or(x$audit$overlap$passed, NA),
    pareto_k = value_or(x$audit$overlap$pareto_k, NA_real_),
    importance_ess = value_or(x$audit$overlap$ess, NA_real_),
    audit_reason = value_or(x$audit$reason, NA_character_))))
  write.csv(index, file.path(root, "fits.csv"), row.names = FALSE)

  comparisons <- list()
  for (i in seq_len(nrow(plan))) {
    case <- plan[i, ]
    ref <- read_record(case, "nuts")
    for (method in methods[-1L]) {
      alt <- read_record(case, method)
      a <- alt$audit$probabilities
      if (is.null(a) || is.null(ref$probability)) next
      names(a)[names(a) == "status"] <- "audit_status"
      b <- ref$probability[, c("covariate", "group", "row", "threshold", "estimate", "mcse")]
      p <- merge(a, b, by = c("covariate", "group", "row", "threshold"), sort = FALSE)
      p$case <- case$case
      p$scenario <- case$scenario
      p$method <- method
      p$fit_status <- alt$status
      p$nuts_status <- ref$status
      p$nuts_estimate <- p$estimate
      p$nuts_mcse <- p$mcse
      p$estimate <- p$mcse <- NULL
      p$abs_error <- abs(p$approximation - p$nuts_estimate)
      p$reference_abs_error <- abs(p$reference - p$nuts_estimate)
      # Independent reference classification includes both NUTS and
      # approximation simulation uncertainty. An ambiguous event stays so.
      p$oracle_status <- "unresolved"
      n <- 8000L; z <- 3
      center <- (p$approximation + z^2 / (2 * n)) / (1 + z^2 / n)
      radius <- z * sqrt(p$approximation * (1 - p$approximation) / n +
        z^2 / (4 * n^2)) / (1 + z^2 / n)
      lo <- center - radius - p$nuts_estimate - z * p$nuts_mcse
      hi <- center + radius - p$nuts_estimate + z * p$nuts_mcse
      eligible <- alt$status == "eligible" & ref$status == "eligible" & is.finite(lo) & is.finite(hi)
      p$oracle_status[which(eligible & (lo > .05 | hi < -.05))] <- "discrepancy"
      p$oracle_status[which(eligible & lo >= -.05 & hi <= .05)] <- "within_tolerance"
      comparisons[[length(comparisons) + 1L]] <- p
    }
  }
  if (length(comparisons)) {
    p <- do.call(rbind, comparisons)
    write.csv(p, file.path(root, "comparisons.csv"), row.names = FALSE)
    print(with(p, table(method, audit_status, oracle_status)))
  }
  print(with(index, table(scenario, method, status)))
}
