# Run each family in one single-thread process; assess exports all planned cases.
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L, args %in% c("gaussian", "nbinom2", "assess", "smoke"))
study <- "tools/benchmarks/bayesian_scales/cv_refits_20261113"
source(file.path(study, "design.R"))
source(file.path(study, "reference.R"))

atomic <- function(x, path) {
  if (file.exists(path)) stop("A record already exists at ", path)
  tmp <- tempfile(tmpdir = dirname(path))
  on.exit(unlink(tmp))
  saveRDS(x, tmp, compress = "gzip")
  if (!file.rename(tmp, path)) stop("Could not publish ", path)
}

if (args != "assess") {
  pkgload::load_all(compile = FALSE, quiet = TRUE)
  stopifnot(packageVersion("RTMB") == "2.0", packageVersion("tmbstan") == "1.2.1")
}

scope <- c(list.files("R", full.names = TRUE), list.files("src",
  pattern = "[.](c|cpp|h)$", full.names = TRUE), "DESCRIPTION", "NAMESPACE",
  file.path(study, c("design.R", "reference.R", "run.R")))
# Use content hashes of lines to make source verification portable across
# Git checkouts with LF and CRLF. Save the native DLL's byte hash separately.
code <- setNames(vapply(scope, function(path) digest::digest(
  readLines(path, warn = FALSE), algo = "sha256"), character(1L)), scope)

models <- c("pooled", "annual")
fold_names <- c("whole_year", "seen_year_new_site")
fit_args <- list(chains = 4L, iter = 3000L, warmup = 1500L, cores = 1L,
  control = list(adapt_delta = .99, max_treedepth = 12L))
budget <- c(1L, 8L, 32L)
quad_order <- c(8L, 16L)
n_stream <- 4L
equal_tol <- 1e-8
block_names <- c("observed", "stress")

attempt <- function(case, input, model, smoke = FALSE) {
  start <- proc.time()[["elapsed"]]
  status <- setNames(rep("failed", length(fold_names)), fold_names)
  if (inherits(input, "error")) return(list(case = case, model = model,
    status = status, error = conditionMessage(input), elapsed = 0))
  fit_setup <- fit_args
  if (smoke) fit_setup <- list(chains = 2L, iter = 80L, warmup = 40L,
    cores = 1L, control = list(adapt_delta = .95))
  warnings <- character()
  cv <- tryCatch(withCallingHandlers(cross_validate_bayes(
    y ~ habitat + offset(offset) + (1 | site), input$data, input$inputs, input$folds,
    scales = if (model == "annual") list(habitat = ~1 + (1 | year)) else NULL,
    family = case$family, row_id = "key", priors = input$priors,
    method = "sample", fit_args = fit_setup, keep_fits = TRUE, seed = case$seed + 1L),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }), error = identity)
  if (inherits(cv, "error")) return(list(case = case, model = model, status = status,
    error = conditionMessage(cv), warnings = warnings, elapsed = proc.time()[["elapsed"]] - start))
  res <- setNames(lapply(seq_along(fold_names), function(k) {
    fold <- fold_names[k]
    base <- list(fold = fold, status = cv$fold_results$status[k],
      diagnostics = cv$diagnostics[[fold]], errors = cv$errors[[fold]],
      warnings = cv$warnings[[fold]], contract = cv$prediction_contracts[[fold]],
      train = cv$training[[fold]], score = cv$joint_scores$neg_log_score[k])
    if (is.null(cv$fits[[fold]])) return(base)
    result <- tryCatch({
      fit <- fit_predictor(cv$fits[[fold]], input, fold, model)
      stopifnot(isTRUE(all.equal(base$train$bin_distances, fit$ref$distances,
        tolerance = 1e-12)))
      pred <- cv$predictions[cv$predictions$fold == fold, ]
      if (is.finite(base$score)) stopifnot(
        identical(pred$new_scale_level, fit$new_scale),
        identical(pred$new_response_level, fit$new_response))
      set.seed(case$seed + 4L + k)
      uniforms <- runif(length(fit$y))
      coarse <- reference_score(fit, uniforms, order = quad_order[1L])
      fine <- reference_score(fit, uniforms, order = quad_order[2L])
      original <- augmented_scores(fit, base$train$score_seed, budgets = 1L)
      streams <- if (fit$uses_latent)
        do.call(rbind, lapply(seq_len(n_stream), function(j)
          augmented_scores(fit, case$seed + 6L + 5L * k + j, budgets = budget))) else
        original
      oracle <- oracle_score(input, fold, fit$y_stress, uniforms)
      fit_error <- abs(c(coarse$observed - fine$observed, coarse$stress - fine$stress,
        max(abs(coarse$pit - fine$pit))))
      score_delta <- original$observed - base$score
      list(reference = fine, reference_refine = setNames(fit_error, c(block_names, "pit")),
        oracle = oracle, original = original, streams = streams,
        package_replay_diff = score_delta, replay_passed = is.finite(score_delta) &&
          abs(score_delta) <= equal_tol,
        truth = input$truth, outcomes = data.frame(visit = pred$row_id, observed = fit$y,
          stress = fit$y_stress, randomized_u = uniforms),
        posterior_n = fit$n, uses_latent = fit$uses_latent)
    }, error = function(e) list(reference_error = conditionMessage(e)))
    c(base, result)
  }), fold_names)
  list(case = case, model = model, status = setNames(vapply(res, "[[", "", "status"),
    fold_names), result = res, warnings = warnings, elapsed = proc.time()[["elapsed"]] - start)
}

if (args == "smoke") {
  for (f in c("gaussian", "nbinom2")) {
    case <- plan[plan$family == f & plan$replicate == 2L, ]
    for (model in models) {
      res <- attempt(case, simulate_cv_case(case), model, smoke = TRUE)
      print(res$status)
      print(res$error)
      print(lapply(res$result, function(r)
        r[c("errors", "reference_error", "score", "package_replay_diff",
          "reference_refine", "uses_latent")]))
      if (!length(res$result) || any(vapply(res$result,
          function(r) !isTRUE(r$replay_passed), logical(1L)))) stop("Smoke replay failed")
    }
  }
  quit(save = "no")
}

records <- file.path(study, "records")
if (args != "assess") {
  dir.create(records, showWarnings = FALSE)
  config <- list(base = "8f7a26a8ead71c29c89e725e40e61829fa11f0e7", code = code,
    plan = plan, model = models, folds = fold_names, fit_args = fit_args,
    budgets = budget, quad_order = quad_order, n_stream = n_stream, replay_tolerance = equal_tol,
    R = R.version.string, locale = Sys.getlocale(), versions = vapply(
      c("RTMB", "TMB", "tmbstan", "rstan", "posterior", "statmod"),
      function(p) as.character(packageVersion(p)), character(1L)),
    dll = digest::digest("src/multiScaleR.dll", file = TRUE, algo = "sha256"))
  conf_file <- file.path(records, paste0("config-", args, ".rds"))
  if (file.exists(conf_file)) stopifnot(identical(config, readRDS(conf_file))) else
    atomic(config, conf_file)
  for (i in which(plan$family == args)) {
    case <- plan[i, , drop = FALSE]
    input <- tryCatch(simulate_cv_case(case), error = identity)
    for (model in models) {
      path <- file.path(records, paste(case$case, model, "rds", sep = "."))
      if (file.exists(path)) next
      res <- attempt(case, input, model)
      atomic(res, path)
      cat(case$case, model, paste(res$status, collapse = " "),
        res$error, round(res$elapsed, 1), "sec\n")
      for (r in res$result) {
        cat(" ", r$fold, r$errors$message, r$reference_error, "replay =", r$package_replay_diff,
          "quadrature =", r$reference_refine, "\n")
      }
      flush.console()
    }
  }
  quit(save = "no")
}

for (f in unique(plan$family)) {
  conf_file <- file.path(records, paste0("config-", f, ".rds"))
  if (file.exists(conf_file)) stopifnot(identical(readRDS(conf_file)$code, code))
}
attempts <- refs <- augmentation <- pits <- list()
for (i in seq_len(nrow(plan))) for (model in models) {
  case <- plan[i, , drop = FALSE]
  path <- file.path(records, paste(case$case, model, "rds", sep = "."))
  a <- if (file.exists(path)) readRDS(path) else list()
  if (length(a)) stopifnot(identical(a$case, case), identical(a$model, model))
  for (fold in fold_names) {
    r <- a$result[[fold]]
    key <- data.frame(case, model = model, fold = fold, row.names = NULL)
    val <- function(x) if (length(x)) x else NA
    attempts[[length(attempts) + 1L]] <- cbind(key, status = if (!length(a)) "pending" else
      a$status[[fold]], stage = val(r$errors$stage),
      error = val(c(r$errors$message, a$error)), reference_error = val(r$reference_error),
      warnings = paste(unique(c(a$warnings, r$warnings)), collapse = " | "),
      package_nll = val(r$score), n_analysis = length(r$train$analysis),
      n_assessment = if (length(r$outcomes$visit)) length(r$outcomes$visit) else NA_integer_,
      replay_error = val(r$package_replay_diff), replay_passed = val(r$replay_passed),
      divergences = val(r$diagnostics$raw$divergences), depth = val(r$diagnostics$raw$max_treedepth_hits),
      fit_elapsed = val(a$elapsed))
    if (!length(r$reference)) next
    for (event in block_names) {
      refs[[length(refs) + 1L]] <- cbind(key, status = r$status, event = event,
        reference_nll = r$reference[[event]], oracle_nll = r$oracle$score[[event]],
        package_seed_nll = r$original[[event]], n = nrow(r$outcomes),
        refinement = r$reference_refine[[event]], oracle_refinement = r$oracle$refine[[event]],
        density_draws = r$reference$effective_density_draws[[event]], posterior_n = r$posterior_n,
        uses_latent = r$uses_latent)
      augmentation[[length(augmentation) + 1L]] <- cbind(key, event = event,
        r$streams[c("budget", "seed")], nll = r$streams[[event]],
        delta = r$streams[[event]] - r$reference[[event]], uses_latent = r$uses_latent)
    }
    pits[[length(pits) + 1L]] <- cbind(key, status = r$status, r$outcomes,
      fit_pit = r$reference$pit, oracle_pit = r$oracle$score$pit,
      fit_refinement = r$reference_refine[["pit"]], oracle_refinement = r$oracle$refine[["pit"]])
  }
}
out <- list(attempts = attempts, references = refs, augmentation = augmentation, pits = pits)
for (nm in names(out)) if (length(out[[nm]]))
  write.csv(do.call(rbind, out[[nm]]), file.path(study, paste0(nm, ".csv")), row.names = FALSE)
