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
    ".msr_gdr_atomic_audit", ".msr_gdr_hash", ".msr_gdr_design_hash", ".msr_gdr_adapter_identity",
    ".msr_gdr_optimizer_checkpoint", "gdr_start_values", ".msr_gdr_exact_start_map",
    ".msr_gdr_TMB_settings", ".msr_gdr_TMB_provider", ".msr_gdr_TMB_diagnostics",
    ".msr_gdr_TMB_refit", ".msr_gdr_refit_start", ".msr_gdr_TMB_preparation", ".msr_gdr_TMB_fit",
    ".msr_gdr_TMB_layout", ".msr_gdr_TMB_policy", ".msr_gdr_TMB_adfun",
    ".msr_gdr_TMB_phase", ".msr_gdr_TMB_point")
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
  if (!is.null(model@TMB)) {
    value <- attr(model, "msr_gdr_tmb_diagnostics")
    if (is.null(value)) return(list(requested = requested, status = "unavailable"))
    value$requested <- requested
    return(value)
  }
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

# Retain only native optim's small numeric/character fields, including failures.
# A checkpoint is evidence of an attempt, never a new fit-eligibility decision.
.msr_gdr_optimizer_checkpoint <- function(optimizer) {
  if (!is.list(optimizer) || is.object(optimizer) || is.null(names(optimizer)) ||
      anyNA(names(optimizer)) || anyDuplicated(names(optimizer)) ||
      !all(c("par", "value", "convergence") %in% names(optimizer)))
    stop("Native optimizer checkpoint metadata is incomplete.", call. = FALSE)
  out <- optimizer[intersect(c("par", "value", "counts", "convergence", "message"), names(optimizer))]
  plain <- function(x) is.null(x) || (!is.object(x) && is.null(dim(x)) &&
    all(names(attributes(x)) %in% "names"))
  if (!all(vapply(out, plain, logical(1))) || !is.numeric(out$par) ||
      !length(out$par) || !is.numeric(out$value) || length(out$value) != 1L ||
      !is.numeric(out$convergence) || length(out$convergence) != 1L ||
      (!is.null(out$counts) && (!is.numeric(out$counts) || length(out$counts) > 64L)) ||
      (!is.null(out$message) && (!is.character(out$message) || length(out$message) != 1L)))
    stop("Native optimizer checkpoints require plain numeric vectors and a scalar message.", call. = FALSE)
  out
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

#' Create an audited native refitter for combined distance-removal models
#'
#' @description Retain a small declared set of coefficient starts while
#'   refitting a GDR model during ecological-scale optimization. The C engine
#'   omits intermediate coefficient Hessians. The TMB engine retains native
#'   joint covariance checks; both engines return an actual native fitted model.
#' @param model An original \code{unmarkedFitGDR} fitted with the native C or
#'   TMB engine in unmarked 1.5.2 (TMB 1.9.25 for TMB fits). Data, formulas,
#'   coefficient identities, actual integration cutoff and evaluated optimizer
#'   settings must remain available. Active offsets, TMB ZIP and stochastic
#'   optimizer controls remain unsupported. No engine is switched silently.
#' @param cold_starts A nonempty named list of finite numeric vectors on the
#'   native link scales. For C models, every vector must have exactly the unique
#'   names in \code{unmarked::coef(model)} and is matched by name. For TMB models,
#'   supply unnamed constructor-order vectors from \code{gdr_start_values(model)$starts},
#'   including every estimated log-SD and excluding latent modes. Alternatively,
#'   name each TMB coordinate as \code{block[within_block]}, using that snapshot's
#'   map (for example, \code{beta_lambda[1]} or \code{lsigma_lambda[1]}).
#'   Named TMB vectors are matched by those unique names. Provide independent
#'   starts distinct from the fitted coordinates. These are neither ecological
#'   kernel scales nor the native integration cutoff K.
#' @param audit_dir A new directory, or an existing empty directory, for compact
#'   audit records. Each refit writes atomic request, validated-input, attempt
#'   and summary RDS files with all attempted starts and diagnostics. An attempt
#'   request is written before entering the native constructor; an interrupted
#'   request is retained without inventing a completion. Raw input data, models and likelihood environments
#'   are hashed in memory or omitted; they are not checkpointed. Keep this
#'   directory to inspect failures, including scalar objective failures.
#' @param provenance Optional named list of caller-supplied source, data or
#'   archive identities. Use simple atomic values and lists. This information
#'   supplements the automatically recorded input and provider fingerprints.
#' @param tmb_settings For a TMB model, an explicit named list containing
#'   \code{method}, \code{se}, \code{threads} and \code{control}. Use a supported
#'   deterministic \code{optim()} method (\code{"BFGS"}, \code{"Nelder-Mead"},
#'   \code{"CG"} or \code{"L-BFGS-B"}), a scalar logical final SE request, \code{threads = 1}
#'   and a named list of finite numeric controls with positive integer
#'   \code{maxit}. Declare the intended evaluated settings; they cannot be
#'   recovered reliably from a saved TMB call. Optional \code{lower}/\code{upper}
#'   bounds require L-BFGS-B and use scalar or native fitted-optimizer order,
#'   which can differ from constructor-start order. Leave \code{NULL} for C models,
#'   whose original evaluated settings are recovered automatically.
#'   Optional \code{random.start} is \code{"cached_best"} (the unchanged native
#'   default) or \code{"zero"}. Optional \code{inner.control} explicitly supplies
#'   positive \code{maxit}, \code{tol}, \code{grad.tol} and \code{step.tol} for
#'   native conditional-mode Newton iterations. These numerical policies are
#'   bound to the audit contract; they do not change priors or prove a rescue.
#'   Optional \code{record_memory = TRUE} requires \pkg{ps} and records process
#'   memory at phase boundaries; its default is \code{FALSE}.
#' @return A function accepting \code{model}, \code{data} and \code{context},
#'   suitable for \code{refit_fn} in \code{\link{multiScale_optim}}. It returns
#'   one native fitted model with a compact \code{msr_gdr_refit_audit} attribute.
#'   Failed calls raise a condition carrying their audit record and its path.
#'   The audit adds `contract_sha256`, identifying the original `contract.rds`,
#'   and `input_context_sha256`, binding the actual input data/design hashes and
#'   authorized spatial columns to that contract. Attempt records retain these
#'   identities, their request/completion paths and the native optimizer's
#'   `par`, `value`, `counts`, `convergence` and `message` when supplied. Values
#'   stay on native link scales; failed/nonfinite diagnostics remain visible.
#'   Each TMB attempt also identifies an `optimizer_point_file`, written before
#'   covariance computation. Its `msr_gdr_optimizer_point` packet contains native
#'   fixed optimizer positions, repeated display labels, data/layout hashes,
#'   finite-cutoff value and storage type, declared policy and phase records.
#'   It is always inference-ineligible and contains no fitted S4 object, latent
#'   states or covariance. A completed point survives a later covariance error;
#'   it does not provide coefficient SEs or permission to resume inference.
#' @details
#' The factory binds to the supplied original model. The callback permits only
#' updates to static spatial columns named by \code{context$covs}. All other
#' native frame fields and nonspatial columns remain fixed. Native engine,
#' actual cutoff K, mixture, distance key, density units, formulas, method,
#' threads and optimizer controls are bound to evaluated settings rather than
#' re-evaluating the original call's potentially changed variables. C bounds
#' retain their existing contract. TMB settings are supplied explicitly because
#' the native saved objective does not retain the original evaluated controls.
#' The unmarked 1.5.2 TMB constructor does not forward its optimizer extras;
#' this adapter forwards the declared controls through a local fitting helper
#' without changing the installed provider or its namespace. Native data,
#' random-effect preparation, compiled likelihood and S4 construction are used.
#' Every TMB attempt creates a fresh tape with zero incoming latent coordinates.
#' Fixed parameters and log-SDs use the declared original/cold starts. Native
#' \code{sdreport()} is evaluated explicitly at the returned optimizer coordinates
#' so that estimates, covariance and optimizer metadata describe the same point.
#'
#' Returned native calls retain the declared controls without capturing the
#' temporary constructor environment. Use this callback for audited refits;
#' ordinary \code{unmarked::update()} retains the provider's TMB control handling.
#'
#' Every call compares the original fitted coefficients with every declared
#' cold start. It retains errors, warnings, native convergence, likelihood,
#' available gradients, coefficient identities and uncertainty status. Only
#' one winning model is retained. The highest finite likelihood must come from
#' a warning-free converged fit. A better but unusable attempt causes an explicit
#' failure rather than an inferior fit being labeled the winner.
#'
#' When \code{context$return_model = FALSE}, C coefficient Hessians are omitted.
#' TMB still computes native \code{sdreport()} for every attempt, including
#' \code{se = FALSE}; this is not a faster likelihood-only TMB interface. TMB
#' eligibility also requires finite scores, positive-definite fixed covariance
#' and, when latent effects are present, positive-definite full joint precision.
#' These checks retain cross-component covariance. Unusable covariance remains
#' a failure even when the final SE request is false. Final calls polish the
#' selected coordinates using the original C or explicitly declared TMB SE request.
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
#' The contract retains original data and evaluated-design fingerprints,
#' numeric constructor-coordinate metadata and the original optimizer summary.
#' Its SHA256 covers the entire contract before its own fingerprint is added,
#' including evaluated settings and native provider identity. This is stronger
#' than the coordinate-only fingerprint from `gdr_start_values()`. Caller
#' provenance is recorded as supplied; it is not independently verified.
#' Input-design fingerprints identify the actual permitted spatial update.
#' None of these records restores an executable objective or evaluates a new
#' coefficient score. Existing fitting, curvature and uncertainty gates still
#' apply, and no Hessians, raw data or design matrices are checkpointed.
#' Restored callbacks verify their adapter and native function identities.
#' Native R and DLL files are hashed once again in each new worker process;
#' a changed provider fails before fitting and leaves an audit record.
#' Identified native GDR response-field models retain their fixed exponential
#' correlation, frozen component designs and full joint covariance through this
#' callback. Supply the same explicit tmb_settings used by fit_gdr_response_field().
#' Ordinary iid update, field simulation, profile and spatial-CV routes remain
#' guarded; the private RTMB fitting admission is unchanged.
#'
#' The native integration cutoff K is distinct from the K column in information
#' criterion tables, which counts fitted parameters. Reusing one cutoff does
#' not replace a separate cutoff-sensitivity analysis.
#'
#' TMB diagnostics evaluate the direct joint gradient and current latent
#' precision at the retained complete mode. They record its maximum latent
#' gradient and precision-scaled Newton decrement, each required below 1e-3,
#' and align the covariance's random coordinates with that mode. The returned
#' marginal value must also agree within 1e-7. These are local consistency checks,
#' not evidence of a unique/global mode or accurate high-dimensional integration.
#' A zero inner start or tighter tolerances may reach another numerical branch;
#' retain every failure and qualify each policy independently. No automatic
#' policy rescue, population-prior change or uncertainty adjustment is applied.
#'
#' Attempt records separate native constructor preparation from TMB tape
#' construction, outer optimization, returned-value checking, covariance and
#' diagnostic qualification. Elapsed times include the declared observer work.
#' With `record_memory = TRUE`, `ps` supplies whole-process RSS and lifetime
#' maximum RSS in bytes at phase boundaries. These include other objects loaded
#' in that R process. They are not isolated tape costs, phase-specific maxima or
#' operating-system memory caps. An unavailable observation retains its reason.
#' Ordinary unmarked constructor replay does not forward these policies in the
#' qualified provider; use this audited refitter to retain the declared controls.
#' Returned native models still retain full native coefficient/latent uncertainty
#' conditional on ecological scales. The numeric point does not replace them.
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
#' @md
#' @export
make_gdr_refitter <- function(model, cold_starts, audit_dir, provenance = list(),
                              tmb_settings = NULL) {
  if (methods::is(model, "unmarkedFitGDRField"))
    return(.msr_gdr_field_refitter(model, cold_starts, audit_dir, provenance, tmb_settings))
  tmb <- methods::is(model, "unmarkedFitGDR") && !is.null(model@TMB)
  if (!tmb && !is.null(tmb_settings))
    stop("tmb_settings apply only to native TMB models.", call. = FALSE)
  settings <- if (tmb) .msr_gdr_TMB_settings(model, tmb_settings) else .msr_gdr_C_settings(model)
  # Native optim can report code zero for a point-only maxit = 0 evaluation.
  maxit <- settings$dots$control$maxit
  if (!is.null(maxit) && (length(maxit) != 1L || !is.numeric(maxit) ||
      !is.finite(maxit) || maxit <= 0 || maxit != floor(maxit)))
    stop("The original GDR optimizer must permit a positive integer maxit.", call. = FALSE)
  coefficients <- .msr_gdr_refit_start(model, settings)
  coefficient_names <- settings$coefficient_names
  start_names <- names(coefficients)
  if (!length(coefficients) || any(!is.finite(coefficients)) || is.null(start_names) ||
      anyNA(start_names) || anyDuplicated(start_names))
    stop("The original GDR model must have finite uniquely named coefficients.", call. = FALSE)
  if (!is.list(cold_starts) || !length(cold_starts) || is.null(names(cold_starts)) ||
      anyNA(names(cold_starts)) || any(!nzchar(names(cold_starts))) ||
      anyDuplicated(names(cold_starts)) || any(names(cold_starts) %in% c("original", "final")))
    stop("cold_starts must be a nonempty uniquely named list; original and final are reserved names.", call. = FALSE)
  cold_starts <- lapply(cold_starts, function(start) {
    if (tmb && is.numeric(start) && is.null(names(start)) && length(start) == length(start_names))
      names(start) <- start_names
    if (!is.numeric(start) || !is.null(dim(start)) || any(!is.finite(start)) ||
        is.null(names(start)) || anyNA(names(start)) || anyDuplicated(names(start)) ||
        !setequal(names(start), start_names))
      stop("Each cold start must have exactly the finite named native coefficients.", call. = FALSE)
    start <- start[start_names]
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
  if (tmb) source_TMB <- .msr_gdr_TMB_provider()
  source <- list(unmarked_version = as.character(utils::packageVersion("unmarked")),
    multiScaleR_version = as.character(utils::packageVersion("multiScaleR")),
    adapter_sha256 = .msr_gdr_adapter_identity(),
    native_function_sha256 = .msr_gdr_hash(list(formals = deparse(formals(unmarked::gdistremoval)),
      body = deparse(body(unmarked::gdistremoval)))),
    native_R_sha256 = digest::digest(file = file.path(find.package("unmarked"), "R", "unmarked.rdb"), algo = "sha256"),
    native_DLL_sha256 = digest::digest(file = getLoadedDLLs()[["unmarked"]][["path"]], algo = "sha256"))
  if (tmb) source$TMB <- source_TMB
  # Full matrices are needed only to verify the original design at creation.
  settings$design <- NULL
  header <- list(schema = 1L, original_model_sha256 = original_key, provider = source,
    metadata = metadata, formulas = lapply(formulas, deparse), controls = settings,
    cold_starts = cold_starts, original_coefficients = coefficients,
    provenance = provenance, models_or_data_checkpointed = FALSE,
    original_data_sha256 = .msr_gdr_hash(original_data), original_design_sha256 = design_hash,
    original_coordinates = gdr_start_values(model),
    original_optimizer = .msr_gdr_optimizer_checkpoint(model@opt))
  contract_sha256 <- .msr_gdr_hash(header)
  header$contract_sha256 <- contract_sha256
  contract_file <- file.path(audit_dir, "contract.rds")
  .msr_gdr_atomic_audit(header, contract_file)
  contract_file_sha256 <- digest::digest(file = contract_file, algo = "sha256")
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
      parameter_count = if (tmb) settings$parameter_count else length(coefficients), input_sha256 = NULL, status = "started",
      contract_sha256 = contract_sha256, contract_file = contract_file,
      contract_file_sha256 = contract_file_sha256,
      input_design_sha256 = NULL, input_context_sha256 = NULL)
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
      }
      if (tmb) {
        actual <- .msr_gdr_TMB_provider(files = audit$provider_check$files_rechecked)
        audit$provider_check$TMB <- actual
        if (!identical(actual, source$TMB[names(actual)]))
          stop("The TMB refitter provider or native TMB files changed.")
      }
      counter$provider_pid <- Sys.getpid()
      .msr_gdr_native_guard(model)
      if (!methods::is(model, "unmarkedFitGDR") || !methods::is(data, "unmarkedFrameGDR"))
        stop("Supply the original GDR model and a native GDR frame.")
      current_metadata <- stats::setNames(lapply(names(metadata), function(name) methods::slot(model, name)), names(metadata))
      current_key <- .msr_gdr_hash(list(data = model@data, formulas = lapply(model@formlist, deparse),
        coefficients = as.list(.msr_gdr_refit_start(model, settings)), metadata = current_metadata))
      if (!identical(current_key, original_key) || !identical(!is.null(model@TMB), tmb) ||
          (tmb && !identical(.msr_gdr_hash(model@TMB$env$data), settings$native_data_sha256)))
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
      if (tmb && !identical(.msr_gdr_TMB_layout(get("getDesign", asNamespace("unmarked"))(data, formulas)),
          settings$design_layout))
        stop("A TMB refit changed native design columns, rows or random-effect layout.")
      if (!identical(settings, header$controls) || !identical(metadata, header$metadata) ||
          !identical(lapply(formulas, deparse), header$formulas) ||
          !identical(cold_starts, header$cold_starts) ||
          !identical(coefficients, header$original_coefficients) ||
          !identical(source, header$provider) || !identical(provenance, header$provenance) ||
          !identical(.msr_gdr_hash(original_data), header$original_data_sha256) ||
          !identical(design_hash, header$original_design_sha256) ||
          !file.exists(contract_file) ||
          !identical(digest::digest(file = contract_file, algo = "sha256"), contract_file_sha256))
        stop("The retained GDR contract or evaluated settings changed; create a new refitter.")
      audit$input_sha256 <- .msr_gdr_hash(data)
      audit$spatial_columns <- covs
      audit$input_design_sha256 <- .msr_gdr_design_hash(data, formulas)
      audit$input_context_sha256 <- .msr_gdr_hash(list(contract_sha256 = contract_sha256,
        input_sha256 = audit$input_sha256, input_design_sha256 = audit$input_design_sha256,
        spatial_columns = covs))
      .msr_gdr_atomic_audit(audit, sub("[.]rds$", "-input.rds", path))
    }, error = function(e) fail(conditionMessage(e)))

    run <- function(start, name, se) {
      record <- list(name = name, start = start, SE_requested = se, warnings = character(),
        error = NULL, converged = FALSE, logLik = NA_real_, coefficients = NULL,
        gradient = NULL, gradient_status = "not_provided_by_native_optimizer",
        covariance = list(requested = se, status = "unavailable"), eligible = FALSE,
        contract_sha256 = contract_sha256, contract_file = contract_file,
        input_context_sha256 = audit$input_context_sha256, input_sha256 = audit$input_sha256,
        input_design_sha256 = audit$input_design_sha256, optimizer = NULL, status = "started")
      record$audit_file <- sub("[.]rds$", paste0("-attempt-", length(audit$attempts) + 1L, ".rds"), path)
      record$request_file <- sub("[.]rds$", "-request.rds", record$audit_file)
      if (tmb) record$native_preparation_started <- .msr_gdr_TMB_phase("native_preparation_started",
        proc.time()[["elapsed"]], settings$policy)
      .msr_gdr_atomic_audit(record, record$request_file)
      fit <- tryCatch(withCallingHandlers({
        args <- c(list(lambdaformula = formulas$lambda, phiformula = formulas$phi,
          distanceformula = formulas$dist, removalformula = formulas$rem, data = data,
          keyfun = metadata$keyfun, output = metadata$output, unitsOut = metadata$unitsOut,
          mixture = metadata$mixture, K = metadata$K, starts = unname(start),
          method = settings$method, se = se, engine = settings$engine, threads = settings$threads),
          if (!tmb) settings$dots)
        checkpoint <- function(point) {
          point$native_preparation <- list(started = record$native_preparation_started,
            elapsed_seconds = point$phases$construction_started$process_elapsed_seconds -
              record$native_preparation_started$process_elapsed_seconds,
            scope = "native_constructor_preparation_before_TMB_tape_helper")
          point$binding <- list(contract_sha256 = contract_sha256, provider = header$provider,
            input_context_sha256 = audit$input_context_sha256, input_sha256 = audit$input_sha256,
            input_design_sha256 = audit$input_design_sha256, start_name = name)
          .msr_gdr_atomic_audit(point, sub("[.]rds$", "-point.rds", record$audit_file))
        }
        value <- if (tmb) .msr_gdr_TMB_refit(args, settings, checkpoint = checkpoint) else do.call(unmarked::gdistremoval, args)
        record$optimizer <- .msr_gdr_optimizer_checkpoint(value@opt)
        if (!identical(names(unmarked::coef(value)), coefficient_names) ||
            .msr_parameter_count(value) != audit$parameter_count ||
            !identical(as.numeric(value@K), as.numeric(metadata$K)) ||
            !identical(!is.null(value@TMB), tmb) || !identical(value@sitesRemoved, model@sitesRemoved))
          stop("A native refit changed coefficient identities, parameter count, engine or cutoff.")
        record$logLik <- as.numeric(unmarked::logLik(value))
        record$convergence <- value@opt$convergence
        record$converged <- identical(as.integer(value@opt$convergence), 0L)
        record$coefficients <- unmarked::coef(value)
        if (tmb) record$constructor_coordinates <- .msr_gdr_refit_start(value, settings)
        record$gradient <- value@opt$gradient
        if (!is.null(record$gradient)) record$gradient_status <- "provided_by_native_optimizer"
        record$covariance <- .msr_gdr_variance_status(value, se)
        if (tmb) record$optimizer_point_file <- record$covariance$workflow$point_file
        if (tmb) record$native_preparation_elapsed_seconds <-
          record$covariance$workflow$phases$construction_started$process_elapsed_seconds -
          record$native_preparation_started$process_elapsed_seconds
        if (tmb) {
          record$gradient <- record$covariance$gradient
          record$gradient_status <- "native_TMB_marginal_score_at_optimizer"
        }
        value
      }, warning = function(w) record$warnings <<- c(record$warnings, conditionMessage(w))),
        error = function(e) {
          if (inherits(e, "msr_gdr_audit_error")) stop(e)
          record$error <<- conditionMessage(e)
          if (inherits(e, "msr_gdr_tmb_error") && !is.null(e$optimizer)) {
            record$optimizer_point_file <<- e$point_file
            record$workflow <<- e$workflow
            record$optimizer <<- e$optimizer
            record$logLik <<- -e$optimizer$value
            record$convergence <<- e$optimizer$convergence
            record$converged <<- identical(as.integer(e$optimizer$convergence), 0L)
            if (!is.null(e$covariance)) record$covariance <<- e$covariance
          }
          NULL
        })
      record$eligible <- !is.null(fit) && is.finite(record$logLik) && record$converged &&
        !length(record$warnings) && ((!se && !tmb) || identical(record$covariance$status, "available"))
      record$status <- if (is.null(record$error)) "returned" else "failed"
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
      if (is.finite(attempt$record$logLik) && attempt$record$logLik > best_value) {
        best <- attempt$fit
        best_value <- attempt$record$logLik
        best_index <- i
      }
      attempt <- NULL
    }
    if (is.na(best_index)) fail("Every declared native coefficient start failed.")
    audit$selected <- best_index
    if (is.null(best) || !audit$attempts[[best_index]]$eligible)
      fail("The highest declared native likelihood is not an eligible converged warning-free fit.")
    if (!audit$objective_only) {
      final <- run(.msr_gdr_refit_start(best, settings), "final", settings$se)
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
