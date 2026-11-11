# Workers rebuild their RTMB tapes. The controller sends only the frozen model
# specification and ordinary R values. Files avoid a local socket requirement.
.msr_bayes_pool <- function(cores, initialize, ...) {
  if (cores <= 1L) return(NULL)
  if (!requireNamespace("processx", quietly = TRUE))
    stop("Install the optional package 'processx' to use cores > 1.", call. = FALSE)
  package_path <- normalizePath(getNamespaceInfo(asNamespace("multiScaleR"), "path"),
    winslash = "/", mustWork = TRUE)
  source <- !file.exists(file.path(package_path, "Meta", "package.rds"))
  script <- file.path(package_path, if (source) "inst" else "", "parallel", "bayes_worker.R")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  path <- tempfile("multiScaleR-bayes-")
  if (!dir.create(path)) stop("Cannot create a temporary directory for Bayesian workers.", call. = FALSE)
  pool <- list(path = path, workers = list())
  ready <- FALSE
  on.exit(if (!ready) .msr_bayes_pool_stop(pool))
  saveRDS(list(paths = .libPaths(), path = package_path, source = source),
    file.path(path, "setup.rds"), compress = FALSE)
  saveRDS(list(initialize = initialize, args = list(...)),
    file.path(path, "context.rds"), compress = FALSE)
  for (i in seq_len(cores)) {
    directory <- file.path(path, as.character(i))
    if (!dir.create(directory)) stop("Cannot create the Bayesian worker directory.", call. = FALSE)
    log <- file.path(directory, "worker.log")
    pool$workers[[i]] <- list(directory = directory, log = log,
      process = processx::process$new(rscript, c("--vanilla", script, path, as.character(i)),
        stdout = log, stderr = "2>&1", cleanup_tree = TRUE, supervise = TRUE))
  }
  .msr_bayes_pool_read(pool$workers, "ready.rds")
  ready <- TRUE
  pool
}

.msr_bayes_pool_stop <- function(pool) {
  for (worker in pool$workers)
    try(worker$process$kill(), silent = TRUE)
  unlink(pool$path, recursive = TRUE)
  invisible(NULL)
}

# Publish a separate marker after saveRDS closes the connection. Readers
# never open data while the other process is still writing it.
.msr_bayes_pool_write <- function(value, file) {
  saveRDS(value, file, compress = FALSE)
  if (!file.create(paste0(file, ".ready")))
    stop("Cannot publish Bayesian worker data at ", basename(file), ".", call. = FALSE)
}

.msr_bayes_pool_take <- function(file) {
  value <- readRDS(file)
  if (unlink(c(file, paste0(file, ".ready"))) != 0L)
    stop("Cannot clear Bayesian worker data at ", basename(file), ".", call. = FALSE)
  value
}

.msr_bayes_pool_read <- function(workers, filename) {
  files <- vapply(workers, function(w) file.path(w$directory, filename), character(1))
  markers <- paste0(files, ".ready")
  results <- vector("list", length(workers))
  pending <- seq_along(workers)
  while (length(pending)) {
    done <- pending[file.exists(markers[pending])]
    for (i in done) {
      reply <- .msr_bayes_pool_take(files[i])
      if (!is.null(reply$error))
        stop("Bayesian worker ", i, " failed: ", reply$error, call. = FALSE)
      for (message in reply$warnings) warning(message, call. = FALSE)
      results[i] <- list(reply$value)
    }
    pending <- setdiff(pending, done)
    for (i in pending) {
      worker <- workers[[i]]
      if (!worker$process$is_alive() && !file.exists(markers[i])) {
        message <- if (file.exists(worker$log)) paste(utils::tail(readLines(worker$log,
          warn = FALSE), 8L), collapse = "\n") else ""
        stop("Bayesian worker ", i, " exited before returning its result. ", message, call. = FALSE)
      }
    }
    if (length(pending)) Sys.sleep(.01)
  }
  results
}

.msr_bayes_pool_map <- function(pool, jobs, task, ...) {
  if (!length(jobs)) return(list())
  n <- min(length(pool$workers), length(jobs))
  workers <- pool$workers[seq_len(n)]
  batches <- parallel::splitIndices(length(jobs), n)
  for (i in seq_len(n))
    .msr_bayes_pool_write(list(jobs = jobs[batches[[i]]], task = task, args = list(...)),
      file.path(workers[[i]]$directory, "job.rds"))
  do.call(c, .msr_bayes_pool_read(workers, "result.rds"))
}

.msr_bayes_worker_main <- function(path, directory) {
  capture <- function(code) {
    warnings <- character()
    reply <- tryCatch(withCallingHandlers(list(value = force(code)),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }), error = function(e) list(error = conditionMessage(e)))
    c(reply, list(warnings = warnings))
  }
  ready <- capture({
    setup <- readRDS(file.path(path, "context.rds"))
    context <- do.call(setup$initialize, setup$args)
    NULL
  })
  .msr_bayes_pool_write(ready, file.path(directory, "ready.rds"))
  if (!is.null(ready$error)) return(invisible(NULL))
  input <- file.path(directory, "job.rds")
  marker <- paste0(input, ".ready")
  repeat {
    if (!file.exists(marker)) {
      Sys.sleep(.01)
      next
    }
    job <- .msr_bayes_pool_take(input)
    reply <- capture(lapply(job$jobs, function(x)
      do.call(job$task, c(list(x, context), job$args))))
    .msr_bayes_pool_write(reply, file.path(directory, "result.rds"))
    rm(job, reply)
  }
}

.msr_bayes_ad_object <- function(spec) {
  RTMB::MakeADFun(.msr_bayes_objective(spec, ad = TRUE), spec$parameters, silent = TRUE)
}

.msr_bayes_nuts_chain <- function(id, obj, iter, warmup, control, seed) {
  # Stan's seed and chain ID determine the trajectory. Seed R separately
  # because RStan uses its RNG to permute the stored post-warmup draws.
  set.seed((as.double(seed) + id - 2) %% .Machine$integer.max + 1)
  initial <- unname(obj$par) + .05 * sin(seq_along(obj$par) + id)
  tmbstan::tmbstan(obj, chains = 1L, chain_id = id, iter = iter,
    warmup = warmup, control = control, seed = seed, init = list(initial),
    cores = 1L, refresh = 0L, laplace = FALSE, silent = TRUE)
}
