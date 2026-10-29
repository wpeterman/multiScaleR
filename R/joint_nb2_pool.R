# Compare numerical helpers before dispatch: development installs can share a
# version number while containing different solvers. No files are written.
.msr_joint_pool_signature <- function() {
  ns <- asNamespace("multiScaleR")
  nms <- sort(ls(ns, all.names = TRUE))
  nms <- nms[startsWith(nms, ".msr_")]
  nms <- nms[vapply(nms, function(n) is.function(get(n, ns)), logical(1))]
  vapply(nms, function(n) paste(deparse(get(n, ns)), collapse = "\n"), character(1))
}

.msr_joint_pool_chunk <- function(ids, payload, bounds, control) {
  obj <- do.call(.msr_joint_nb2_objective, payload)
  bounds$starts <- bounds$starts[ids, , drop = FALSE]
  fit <- .msr_joint_nb2_multistart(obj, bounds, control)
  fit$records$start <- ids
  list(ids = ids, fits = fit$fits, records = fit$records,
       coefficient_audit = if (is.function(obj$coefficient_audit)) obj$coefficient_audit() else NULL)
}

.msr_joint_nb2_pool <- function(obj, bounds, control, pool) {
  payload <- list(engine = obj$engine, y = obj$y, family = obj$family,
    X_fixed = obj$X_fixed, offset = obj$offset, weights = obj$weights,
    trials = obj$trials, coefficient_start = obj$coefficient_start)
  .msr_joint_pool_ram(payload, length(pool))
  signature <- .msr_joint_pool_signature()
  # Do not stop or reconfigure a caller-owned cluster, including on errors.
  bootstrap <- function(paths, expected) {
    .libPaths(paths)
    ns <- loadNamespace("multiScaleR")
    f <- get0(".msr_joint_pool_signature", ns, inherits = FALSE)
    is.function(f) && identical(f(), expected)
  }
  # A namespace-bound closure would load the worker's default installation
  # during deserialization, before .libPaths() can select the caller's build.
  environment(bootstrap) <- baseenv()
  compatible <- parallel::clusterCall(pool, bootstrap, .libPaths(), signature)
  if (!all(unlist(compatible)))
    stop("Worker package code differs from the controller. Install the same multiScaleR build and recreate worker_pool.", call. = FALSE)
  ids <- split(seq_len(nrow(bounds$starts)),
    rep(seq_len(min(length(pool), nrow(bounds$starts))), length.out = nrow(bounds$starts)))
  parts <- parallel::parLapply(pool, ids, .msr_joint_pool_chunk, payload, bounds, control)
  # Reestablish serial ordering before ties are resolved or diagnostics computed.
  fits <- vector("list", nrow(bounds$starts))
  for (part in parts) fits[part$ids] <- part$fits
  records <- do.call(rbind, lapply(parts, `[[`, "records"))
  records <- records[order(records$start), , drop = FALSE]
  rownames(records) <- NULL
  good <- which(records$eligible)
  winner <- if (length(good)) good[which.min(records$value[good])] else NA_integer_
  list(best = if (is.na(winner)) NULL else fits[[winner]], winner = winner,
       fits = fits, records = records, starts = bounds$starts,
       worker_coefficient_audits = lapply(parts, `[[`, "coefficient_audit"))
}

# Include binned matrices even for legacy list inputs that cannot use the
# public kernel_prep RAM estimator. This is a planning guard, not a RAM quota.
.msr_joint_pool_ram <- function(payload, workers) {
  total <- .msr_total_ram_bytes()
  estimate <- (workers + 1) * (4 * as.numeric(utils::object.size(payload)) + 1024^3)
  if (is.finite(total) && estimate > .75 * total)
    stop("worker_pool payload exceeds the conservative RAM budget; use fewer workers or serial fitting.", call. = FALSE)
  invisible(estimate)
}
