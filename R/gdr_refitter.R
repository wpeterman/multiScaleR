# Read the evaluated C settings, never expressions in a mutable original call.
.msr_gdr_C_settings <- function(model) {
  .msr_gdr_native_guard(model)
  if (!methods::is(model, "unmarkedFitGDR") || !is.null(model@TMB))
    stop("Audited GDR refits currently require the native C engine.", call. = FALSE)
  if (as.character(utils::packageVersion("unmarked")) != "1.5.2")
    stop("This native GDR refit adapter requires separate qualification for this unmarked version.", call. = FALSE)
  ns <- asNamespace("unmarked")
  definitions <- list()
  find_nll <- function(expr) {
    if (!is.call(expr)) return(invisible(NULL))
    if (identical(expr[[1L]], as.name("<-")) && identical(expr[[2L]], as.name("nll")) &&
        is.call(expr[[3L]]) && identical(expr[[3L]][[1L]], as.name("function")))
      definitions[[length(definitions) + 1L]] <<- expr[[3L]]
    for (part in as.list(expr)[-1L]) if (is.call(part)) find_nll(part)
    invisible(NULL)
  }
  find_nll(body(get("gdistremoval", ns)))
  if (length(definitions) != 1L)
    stop("Could not identify the qualified native C likelihood definition.", call. = FALSE)
  reference <- eval(definitions[[1L]], ns)
  if (!identical(formals(model@nllFun), formals(reference)) ||
      !identical(body(model@nllFun), body(reference)))
    stop("The model does not retain the qualified native C likelihood.", call. = FALSE)
  env <- environment(model@nllFun)
  get_binding <- function(name) {
    if (!exists(name, envir = env, inherits = FALSE))
      stop("Original native GDR setting is unavailable: ", name, ".", call. = FALSE)
    get(name, envir = env, inherits = FALSE)
  }
  if (!identical(get_binding("data"), model@data) ||
      !identical(get_binding("formulas"), model@formlist))
    stop("Stored native GDR data or formulas do not match the original likelihood.", call. = FALSE)
  for (name in c("K", "keyfun", "mixture", "output", "unitsOut"))
    if (!identical(as.vector(get_binding(name)), as.vector(methods::slot(model, name))))
      stop("Stored native GDR setting differs from its likelihood: ", name, ".", call. = FALSE)
  out <- stats::setNames(lapply(c("method", "se", "engine", "threads"), get_binding),
    c("method", "se", "engine", "threads"))
  if (!identical(out$engine, "C") || !is.logical(out$se) || length(out$se) != 1L ||
      is.na(out$se) || !out$method %in% c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B", "Brent"))
    stop("Original native C optimizer settings are not qualified for deterministic refits.", call. = FALSE)
  validate_scalar_numeric(out$threads, "original native threads", positive = TRUE, integerish = TRUE)
  # Build the expression explicitly: these dots belong to the saved likelihood
  # environment, rather than this inspector's formal arguments.
  out$dots <- eval(call("list", as.name("...")), envir = env)
  if (length(out$dots) && (is.null(names(out$dots)) || anyNA(names(out$dots)) ||
      anyDuplicated(names(out$dots)) || any(!names(out$dots) %in% c("control", "lower", "upper"))))
    stop("Original native optimizer extras need separate data-binding qualification; supported extras are control, lower and upper.", call. = FALSE)
  out$K <- model@K
  out$coefficient_names <- names(unmarked::coef(model))
  out$design <- get_binding("gd")[c("X_lambda", "X_phi", "X_dist", "X_rem")]
  out
}

.msr_gdr_hash <- function(object)
  digest::digest(object, algo = "sha256", serializeVersion = 2)

.msr_gdr_adapter_identity <- function() {
  functions <- c("make_gdr_refitter", ".msr_gdr_C_settings", ".msr_gdr_variance_status",
    ".msr_gdr_atomic_audit", ".msr_gdr_hash", ".msr_gdr_design_hash", ".msr_gdr_adapter_identity")
  .msr_gdr_hash(lapply(functions, function(name) {
    fun <- get(name, envir = environment(.msr_gdr_adapter_identity), inherits = FALSE)
    list(name = name, formals = deparse(formals(fun)), body = deparse(body(fun)))
  }))
}

.msr_gdr_design_hash <- function(data, formulas) {
  design <- get("getDesign", asNamespace("unmarked"))(data, formulas)
  .msr_gdr_hash(design[c("X_lambda", "X_phi", "X_dist", "X_rem")])
}

.msr_gdr_variance_status <- function(model, requested) {
  if (!requested) return(list(requested = FALSE, status = "not_requested"))
  n <- length(unmarked::coef(model))
  check <- .msr_refit_hessian_se(model@opt$hessian, rep(TRUE, n), rep(1, n))
  native <- vapply(model@estimates@estimates, function(estimate) {
    v <- diag(estimate@covMat)
    length(v) == length(estimate@estimates) && all(is.finite(v)) && all(v > 0)
  }, logical(1))
  list(requested = TRUE, status = if (check$diagnostic$status == "available" && all(native))
    "available" else "unusable", hessian_status = check$diagnostic$status,
    native_component_covariance_finite_positive = all(native),
    uncertainty = "native_coefficients_conditional_on_spatial_scales")
}

.msr_gdr_atomic_audit <- function(value, path) {
  temporary <- paste0(path, ".tmp")
  abort <- function(message) stop(structure(list(message = message, call = NULL,
    audit = value, audit_file = path), class = c("msr_gdr_audit_error", "error", "condition")))
  if (file.exists(path) || file.exists(temporary))
    abort(paste0("Retain the existing GDR audit record: ", path))
  tryCatch(saveRDS(value, temporary, version = 2), error = function(e)
    abort(paste0("Could not write the GDR audit record: ", temporary, ". ", conditionMessage(e))))
  if (!file.rename(temporary, path))
    abort(paste0("Could not finish the atomic GDR audit record: ", temporary))
  path
}

#' Create an audited native C refitter for combined distance-removal models
#'
#' @description Retain a small declared set of coefficient starts while
#'   refitting a GDR model during ecological-scale optimization. Intermediate
#'   fits omit coefficient Hessians; the returned final fit uses the original
#'   evaluated standard-error setting.
#' @param model An original \code{unmarkedFitGDR} fitted with the native C engine
#'   in unmarked 1.5.2. Data, formulas, coefficient identities, actual integration
#'   cutoff and evaluated optimizer settings must remain recoverable. The TMB
#'   engine, active offsets and stochastic optimizer controls need separate
#'   qualification. Existing generic refit interfaces remain available.
#' @param cold_starts A nonempty named list of finite numeric vectors on the
#'   native coefficient link scales. Every vector must have exactly the unique
#'   names in \code{unmarked::coef(model)}. Their order is matched by name.
#'   Provide independent starts distinct from the fitted coefficients. These
#'   are coefficient starts, not ecological kernel scales or native cutoff K.
#' @param audit_dir A new directory, or an existing empty directory, for compact
#'   audit records. Each refit writes atomic request, attempt and summary RDS
#'   files with all attempted starts and diagnostics. Raw input data, models and likelihood environments
#'   are hashed in memory or omitted; they are not checkpointed. Keep this
#'   directory to inspect failures, including scalar objective failures.
#' @param provenance Optional named list of caller-supplied source, data or
#'   archive identities. Use simple atomic values and lists. This information
#'   supplements the automatically recorded input and provider fingerprints.
#' @return A function accepting \code{model}, \code{data} and \code{context},
#'   suitable for \code{refit_fn} in \code{\link{multiScale_optim}}. It returns
#'   one native fitted model with a compact \code{msr_gdr_refit_audit} attribute.
#'   Failed calls raise a condition carrying their audit record and its path.
#' @details
#' The factory binds to the supplied original model. The callback permits only
#' updates to static spatial columns named by \code{context$covs}. All other
#' native frame fields and nonspatial columns remain fixed. Native engine,
#' actual cutoff K, mixture, distance key, density units, formulas, method,
#' threads, bounds and optimizer controls are retained from evaluated settings,
#' rather than re-evaluating the original call's potentially changed variables.
#'
#' Every call compares the original fitted coefficients with every declared
#' cold start. It retains errors, warnings, native convergence, likelihood,
#' available gradients, coefficient identities and uncertainty status. Only
#' one winning model is retained. The highest finite likelihood must come from
#' a warning-free converged fit. A better but unusable attempt causes an explicit
#' failure rather than an inferior fit being labeled the winner.
#'
#' When \code{context$return_model = FALSE}, coefficient Hessians are omitted.
#' Final calls polish the selected coefficients using the original SE request.
#' A final likelihood loss greater than \code{1e-5}, inconsistent parameter
#' count, or unusable requested native covariance fails with its audit retained.
#' Kernel-conditional coefficient SEs do not include uncertainty in estimated
#' spatial scales. This helper does not add Bayesian priors, bootstrap intervals
#' or a guarantee that unvisited likelihood modes have been found.
#' An explicitly supplied optimizer iteration budget must be positive;
#' a zero-iteration evaluation is not accepted as a converged refit.
#'
#' Audit input fingerprints use SHA256 in memory through digest. The records
#' contain coefficient vectors and summaries, not candidate-model arrays.
#' Restored callbacks verify their adapter and native function identities.
#' Native R and DLL files are hashed once again in each new worker process;
#' a changed provider fails before fitting and leaves an audit record.
#' The native integration cutoff K is distinct from the K column in information
#' criterion tables, which counts fitted parameters. Reusing one cutoff does
#' not replace a separate cutoff-sensitivity analysis.
#' @seealso \code{\link{multiScale_optim}}, \code{\link{validate_unmarked_data}},
#'   \code{\link{aic_tab}}
#' @examples
#' if (as.character(utils::packageVersion("unmarked")) == "1.5.2") {
#'   yd <- rbind(c(0, 0), c(1, 0), c(0, 2), c(2, 1), c(0, 0))
#'   yr <- rbind(c(0, 0, 0), c(0, 1, 0), c(1, 0, 1), c(1, 1, 1), c(0, 0, 0))
#'   frame <- unmarked::unmarkedFrameGDR(yd, yr, dist.breaks = c(0, 25, 75),
#'     unitsIn = "m", period.lengths = c(1, 2, 1),
#'     siteCovs = data.frame(hab = seq(-1, 1, length.out = 5)))
#'   # Fit the original model before declaring independent refit starts.
#'   model <- unmarked::gdistremoval(lambdaformula = ~hab, data = frame,
#'     output = "abund", K = 40, se = FALSE, engine = "C",
#'     starts = c(log(8), 0, log(40), qlogis(0.3)), control = list(maxit = 2000))
#'   cold <- unmarked::coef(model)
#'   cold[] <- 0
#'   cold["dist(Int)"] <- log(30)
#'   callback <- make_gdr_refitter(model, list(cold = cold),
#'     audit_dir = tempfile("msr-gdr-audit-"))
#'   is.function(callback)
#'   # Supply callback as refit_fn in multiScale_optim().
#' }
#' @export
make_gdr_refitter <- function(model, cold_starts, audit_dir, provenance = list()) {
  settings <- .msr_gdr_C_settings(model)
  # Native optim can report code zero for a point-only maxit = 0 evaluation.
  maxit <- settings$dots$control$maxit
  if (!is.null(maxit) && (length(maxit) != 1L || !is.numeric(maxit) ||
      !is.finite(maxit) || maxit <= 0 || maxit != floor(maxit)))
    stop("The original GDR optimizer must permit a positive integer maxit.", call. = FALSE)
  coefficients <- unmarked::coef(model)
  coefficient_names <- settings$coefficient_names
  if (!length(coefficients) || any(!is.finite(coefficients)) || is.null(coefficient_names) ||
      anyNA(coefficient_names) || anyDuplicated(coefficient_names))
    stop("The original GDR model must have finite uniquely named coefficients.", call. = FALSE)
  if (!is.list(cold_starts) || !length(cold_starts) || is.null(names(cold_starts)) ||
      anyNA(names(cold_starts)) || any(!nzchar(names(cold_starts))) ||
      anyDuplicated(names(cold_starts)) || any(names(cold_starts) %in% c("original", "final")))
    stop("cold_starts must be a nonempty uniquely named list; original and final are reserved names.", call. = FALSE)
  cold_starts <- lapply(cold_starts, function(start) {
    if (!is.numeric(start) || !is.null(dim(start)) || any(!is.finite(start)) ||
        is.null(names(start)) || anyNA(names(start)) || anyDuplicated(names(start)) ||
        !setequal(names(start), coefficient_names))
      stop("Each cold start must have exactly the finite named native coefficients.", call. = FALSE)
    start <- start[coefficient_names]
    if (identical(as.numeric(start), as.numeric(coefficients)))
      stop("A cold start must be independent of the original fitted coefficient vector.", call. = FALSE)
    start
  })
  if (anyDuplicated(vapply(cold_starts, .msr_gdr_hash, character(1))))
    stop("Cold coefficient starts must be distinct from one another.", call. = FALSE)
  simple <- function(value) {
    if (is.list(value) && !is.object(value)) return(all(vapply(value, simple, logical(1))))
    is.null(value) || (is.atomic(value) && !is.object(value) && is.null(dim(value)) && length(value) <= 64L)
  }
  if (!is.list(provenance) || (length(provenance) && (is.null(names(provenance)) ||
      anyNA(names(provenance)) || any(!nzchar(names(provenance))) || anyDuplicated(names(provenance)))) ||
      !simple(provenance) || length(serialize(provenance, NULL, version = 2)) > 65536L)
    stop("provenance must be a compact named list of simple atomic values and lists (at most 64 values per vector and 64 KiB total).", call. = FALSE)
  validate_character_scalar(audit_dir, "audit_dir")
  if (dir.exists(audit_dir) && length(list.files(audit_dir, all.files = TRUE, no.. = TRUE)))
    stop("audit_dir must be new or empty; retain existing records.", call. = FALSE)
  if (!dir.exists(audit_dir) && !dir.create(audit_dir, recursive = TRUE))
    stop("Could not create audit_dir.", call. = FALSE)
  audit_dir <- normalizePath(audit_dir, winslash = "/", mustWork = TRUE)
  formulas <- model@formlist
  original_data <- model@data
  metadata <- lapply(c("mixture", "keyfun", "output", "unitsOut", "K"),
    function(name) methods::slot(model, name))
  names(metadata) <- c("mixture", "keyfun", "output", "unitsOut", "K")
  design_hash <- .msr_gdr_hash(settings$design)
  if (!identical(.msr_gdr_design_hash(original_data, formulas), design_hash))
    stop("Current training transformations differ from the retained native design; restore the original training settings.", call. = FALSE)
  original_key <- .msr_gdr_hash(list(data = original_data, formulas = lapply(formulas, deparse),
    coefficients = as.list(coefficients), metadata = metadata))
  source <- list(unmarked_version = as.character(utils::packageVersion("unmarked")),
    multiScaleR_version = as.character(utils::packageVersion("multiScaleR")),
    adapter_sha256 = .msr_gdr_adapter_identity(),
    native_function_sha256 = .msr_gdr_hash(list(formals = deparse(formals(unmarked::gdistremoval)),
      body = deparse(body(unmarked::gdistremoval)))),
    native_R_sha256 = digest::digest(file = file.path(find.package("unmarked"), "R", "unmarked.rdb"), algo = "sha256"),
    native_DLL_sha256 = digest::digest(file = getLoadedDLLs()[["unmarked"]][["path"]], algo = "sha256"))
  # Full matrices are needed only to verify the original design at creation.
  settings$design <- NULL
  header <- list(schema = 1L, original_model_sha256 = original_key, provider = source,
    metadata = metadata, formulas = lapply(formulas, deparse), controls = settings,
    cold_starts = cold_starts, original_coefficients = coefficients,
    provenance = provenance, models_or_data_checkpointed = FALSE)
  .msr_gdr_atomic_audit(header, file.path(audit_dir, "contract.rds"))
  rm(model)
  counter <- new.env(parent = emptyenv())
  counter$n <- 0L
  counter$provider_pid <- Sys.getpid()
  function(model, data, context = list(return_model = TRUE)) {
    counter$n <- counter$n + 1L
    stamp <- gsub("[^0-9]", "", format(Sys.time(), "%Y%m%d%H%M%OS6", tz = "UTC"))
    path <- file.path(audit_dir, paste0("refit-", Sys.getpid(), "-", stamp, "-", counter$n, ".rds"))
    audit <- list(schema = 1L, original_model_sha256 = original_key, provider = source,
      metadata = metadata, provenance = provenance, attempts = list(), selected = NULL,
      final_SE_requested = settings$se, objective_only = NULL,
      parameter_count = length(coefficients), input_sha256 = NULL, status = "started")
    fail <- function(message) {
      audit$status <- "failed"
      audit$error <- message
      audit$path <- path
      .msr_gdr_atomic_audit(audit, path)
      stop(structure(list(message = paste0(message, " Audit: ", path), call = NULL,
        audit = audit, audit_file = path), class = c("msr_gdr_refit_error", "error", "condition")))
    }
    .msr_gdr_atomic_audit(audit, sub("[.]rds$", "-request.rds", path))
    tryCatch({
      if (!is.list(context) || (!is.null(context$return_model) &&
          (!is.logical(context$return_model) || length(context$return_model) != 1L || is.na(context$return_model))))
        stop("context must contain a scalar logical return_model when supplied.")
      audit$objective_only <- identical(context$return_model, FALSE)
      audit$provider_check <- list(unmarked_version = as.character(utils::packageVersion("unmarked")),
        adapter_sha256 = .msr_gdr_adapter_identity(),
        native_function_sha256 = .msr_gdr_hash(list(formals = deparse(formals(unmarked::gdistremoval)),
          body = deparse(body(unmarked::gdistremoval)))), files_rechecked = FALSE)
      if (!identical(audit$provider_check$unmarked_version, source$unmarked_version) ||
          !identical(audit$provider_check$adapter_sha256, source$adapter_sha256) ||
          !identical(audit$provider_check$native_function_sha256, source$native_function_sha256))
        stop("The GDR refitter no longer matches its qualified provider or adapter.")
      # Recheck binaries on a restored process, then reuse the verified result.
      # Function identities above are inexpensive and checked on every call.
      if (!identical(counter$provider_pid, Sys.getpid())) {
        audit$provider_check$files_rechecked <- TRUE
        audit$provider_check$native_R_sha256 <- digest::digest(file =
          file.path(find.package("unmarked"), "R", "unmarked.rdb"), algo = "sha256")
        audit$provider_check$native_DLL_sha256 <- digest::digest(file =
          getLoadedDLLs()[["unmarked"]][["path"]], algo = "sha256")
        if (!identical(audit$provider_check$native_R_sha256, source$native_R_sha256) ||
            !identical(audit$provider_check$native_DLL_sha256, source$native_DLL_sha256))
          stop("The restored GDR native files differ from the qualified provider.")
        counter$provider_pid <- Sys.getpid()
      }
      .msr_gdr_native_guard(model)
      if (!methods::is(model, "unmarkedFitGDR") || !methods::is(data, "unmarkedFrameGDR"))
        stop("Supply the original GDR model and a native GDR frame.")
      current_metadata <- stats::setNames(lapply(names(metadata), function(name) methods::slot(model, name)), names(metadata))
      current_key <- .msr_gdr_hash(list(data = model@data, formulas = lapply(model@formlist, deparse),
        coefficients = as.list(unmarked::coef(model)), metadata = current_metadata))
      if (!identical(current_key, original_key) || !is.null(model@TMB))
        stop("The callback is bound to a different original model; create a new refitter.")
      covs <- context$covs
      if (is.null(covs)) covs <- character()
      if (!is.character(covs) || anyNA(covs) || anyDuplicated(covs) ||
          any(!covs %in% names(original_data@siteCovs)))
        stop("context$covs must name original static site columns.")
      .msr_gdr_spatial_guard(model, covs)
      if (!identical(names(data@siteCovs), names(original_data@siteCovs)) ||
          !identical(rownames(data@siteCovs), rownames(original_data@siteCovs)) ||
          !identical(data@siteCovs[setdiff(names(data@siteCovs), covs)],
            original_data@siteCovs[setdiff(names(original_data@siteCovs), covs)]))
        stop("A refit changed nonspatial site columns or site order.")
      check_data <- data
      check_data@siteCovs <- original_data@siteCovs
      if (!identical(check_data, original_data))
        stop("A refit changed histories, protocol, yearly/observation data or frame layout.")
      if (!identical(.msr_gdr_design_hash(original_data, formulas), design_hash))
        stop("The original training transformations changed after the refitter was created.")
      audit$input_sha256 <- .msr_gdr_hash(data)
      audit$spatial_columns <- covs
    }, error = function(e) fail(conditionMessage(e)))

    run <- function(start, name, se) {
      record <- list(name = name, start = start, SE_requested = se, warnings = character(),
        error = NULL, converged = FALSE, logLik = NA_real_, coefficients = NULL,
        gradient = NULL, gradient_status = "not_provided_by_native_optimizer",
        covariance = list(requested = se, status = "unavailable"), eligible = FALSE)
      fit <- tryCatch(withCallingHandlers({
        args <- c(list(lambdaformula = formulas$lambda, phiformula = formulas$phi,
          distanceformula = formulas$dist, removalformula = formulas$rem, data = data,
          keyfun = metadata$keyfun, output = metadata$output, unitsOut = metadata$unitsOut,
          mixture = metadata$mixture, K = metadata$K, starts = unname(start),
          method = settings$method, se = se, engine = "C", threads = settings$threads), settings$dots)
        value <- do.call(unmarked::gdistremoval, args)
        if (!identical(names(unmarked::coef(value)), coefficient_names) ||
            .msr_parameter_count(value) != length(coefficients) ||
            !identical(as.numeric(value@K), as.numeric(metadata$K)) || !is.null(value@TMB))
          stop("A native refit changed coefficient identities, parameter count, engine or cutoff.")
        record$logLik <- as.numeric(unmarked::logLik(value))
        record$convergence <- value@opt$convergence
        record$converged <- identical(as.integer(value@opt$convergence), 0L)
        record$coefficients <- unmarked::coef(value)
        record$gradient <- value@opt$gradient
        if (!is.null(record$gradient)) record$gradient_status <- "provided_by_native_optimizer"
        record$covariance <- .msr_gdr_variance_status(value, se)
        value
      }, warning = function(w) record$warnings <<- c(record$warnings, conditionMessage(w))),
        error = function(e) { record$error <<- conditionMessage(e); NULL })
      record$eligible <- !is.null(fit) && is.finite(record$logLik) && record$converged &&
        !length(record$warnings) && (!se || identical(record$covariance$status, "available"))
      record$audit_file <- sub("[.]rds$", paste0("-attempt-", length(audit$attempts) + 1L, ".rds"), path)
      .msr_gdr_atomic_audit(record, record$audit_file)
      list(fit = fit, record = record)
    }
    starts <- c(list(original = coefficients), cold_starts)
    best <- NULL
    best_value <- -Inf
    best_index <- NA_integer_
    for (i in seq_along(starts)) {
      attempt <- run(starts[[i]], names(starts)[i], FALSE)
      audit$attempts[[i]] <- attempt$record
      if (!is.null(attempt$fit) && is.finite(attempt$record$logLik) && attempt$record$logLik > best_value) {
        best <- attempt$fit
        best_value <- attempt$record$logLik
        best_index <- i
      }
      attempt <- NULL
    }
    if (is.null(best)) fail("Every declared native coefficient start failed.")
    audit$selected <- best_index
    if (!audit$attempts[[best_index]]$eligible)
      fail("The highest declared native likelihood is not an eligible converged warning-free fit.")
    if (!audit$objective_only) {
      final <- run(unmarked::coef(best), "final", settings$se)
      audit$attempts[[length(audit$attempts) + 1L]] <- final$record
      if (!final$record$eligible || !is.finite(final$record$logLik) ||
          final$record$logLik < best_value - 1e-5)
        fail("The final native refit lost likelihood or has unusable requested uncertainty.")
      best <- final$fit
    }
    audit$status <- "returned"
    audit$path <- path
    .msr_gdr_atomic_audit(audit, path)
    attr(best, "msr_gdr_refit_audit") <- audit
    best
  }
}
