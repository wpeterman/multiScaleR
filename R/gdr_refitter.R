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

#' Create an Audited Refit Function for Distance-Removal Models
#'
#' Build a refit function for `multiScale_optim(refit_fn = ...)` that refits
#' an `unmarked` distance-removal (GDR) model at each candidate scale of effect
#' from its original estimates and from extra starting values you declare. It
#' keeps the best clean fit and writes a record of every attempt to an audit
#' directory. A distance-removal likelihood can have more than one local
#' optimum, and several starts make it less likely that a poor optimum is
#' accepted unnoticed.
#'
#' @param model The starting model. One of:
#'   * an `unmarkedFitGDR` fitted by `unmarked::gdistremoval()` with
#'     `engine = "C"`;
#'   * an `unmarkedFitGDR` fitted with `engine = "TMB"`, including models with
#'     random intercepts;
#'   * an [unmarkedFitGDRField][unmarkedFitGDRField-class] from
#'     [fit_gdr_response_field()].
#'
#'   The refit function is tied to this model. It reuses the model's data,
#'   formulas, cutoff `K`, mixture, distance function, output units and
#'   optimizer settings, and it never switches engine. The model needs
#'   `unmarked` 1.5.2, plus `TMB` 1.9.25 for the TMB engine. Formulas with
#'   `offset()` terms and zero-inflated TMB models (`mixture = "ZIP"` with
#'   `engine = "TMB"`) are not supported.
#' @param cold_starts Named list of one or more extra starting vectors on the
#'   link scales of the model, for example `list(cold = cold)`. For a C-engine
#'   model, name each vector exactly like `unmarked::coef(model)`, such as
#'   `"lambda(Int)"` and `"dist(Int)"`. For a TMB-engine or field model, supply
#'   one value for every fixed coefficient and log SD, either unnamed in the
#'   order of `gdr_start_values(model)$starts`, or named by parameter block and
#'   position, such as `"beta_lambda[1]"` or `"lsigma_lambda[1]"`. Random
#'   effects are not included. Each vector must differ from the model's
#'   estimates and from the other vectors, and the names `"original"` and
#'   `"final"` are reserved. These values are model coefficients, not kernel
#'   scales of effect and not the cutoff `K`.
#' @param audit_dir Character string naming the directory for audit records.
#'   It must be new or empty, and it is created if needed. (For a field model,
#'   an existing directory with records is also accepted.) Keep it, because the
#'   records are the only trace of failed attempts. See "Read the audit
#'   directory".
#' @param provenance Optional named list of your own identifiers, such as a
#'   data source, a script version or an archive DOI. Use simple values: atomic
#'   vectors of at most 64 elements, or nested lists of them, with at most
#'   64 KiB in total. It is copied into every record as supplied and is not
#'   checked.
#' @param tmb_settings For a TMB-engine or field model, a named list of the
#'   optimizer settings to use for every refit. Leave it `NULL` (the default)
#'   for a C-engine model, whose settings are read from the fitted model.
#'   Required elements:
#'   * `method`: one of `"BFGS"`, `"Nelder-Mead"`, `"CG"` or `"L-BFGS-B"`.
#'   * `se`: logical. Whether the final refit of each call computes standard
#'     errors.
#'   * `threads`: must be 1.
#'   * `control`: named list of finite numeric `optim()` controls with a
#'     positive integer `maxit`, such as
#'     `list(maxit = 2000L, reltol = 1e-10)`.
#'
#'   Optional elements:
#'   * `lower`, `upper`: parameter bounds for `method = "L-BFGS-B"`, either one
#'     number or one value per coordinate in the fitted optimizer order
#'     (`names(gdr_start_values(model)$fitted_values)`), which can differ from
#'     the start order.
#'   * `random.start`: how TMB starts each search for the random-effect values.
#'     `"cached_best"` (the TMB default) starts from the best values found so
#'     far. `"zero"` starts from zero.
#'   * `inner.control`: a list with positive `maxit`, `tol`, `grad.tol` and
#'     `step.tol` for that inner search.
#'   * `record_memory`: `TRUE` records the memory used by the R process at each
#'     phase of a fit. It requires the `ps` package. The default is `FALSE`.
#'
#'   For a field model, pass the identical list used in
#'   [fit_gdr_response_field()].
#' @return A function with arguments `model`, `data` and `context`, for the
#'   `refit_fn` argument of [multiScale_optim()]. The optimizer calls it with
#'   the starting model, a copy of its frame in which the kernel-weighted
#'   covariates have new values, and a context whose `covs` element names
#'   those columns. Each call returns one fitted model of the same class as
#'   `model`, with an attribute `msr_gdr_refit_audit` that summarizes the call:
#'   * `status`: `"returned"` for a successful call.
#'   * `attempts`: one record per start (`"original"`, each cold start, and
#'     `"final"`), each with the `start`, `converged`, the `logLik`,
#'     `coefficients` on the link scale, `warnings`, any `error`, the
#'     `covariance` status, whether it was `eligible`, and the optimizer
#'     output (`par`, `value`, `counts`, `convergence`, `message`). TMB
#'     attempts add `marginal_score` (largest absolute gradient, its 0.001
#'     threshold and whether it passed) and `optimizer_point_file`.
#'   * `selected`: the index of the start with the highest log-likelihood.
#'   * `objective_only`: `TRUE` when the call skipped the final refit.
#'   * `contract_sha256` and `input_context_sha256`: fingerprints that tie the
#'     call to `contract.rds` and to the inputs it received. `path` is the
#'     summary record.
#'
#'   A field model returns the record described in [fit_gdr_response_field()]
#'   instead. A failed call stops with an error of class
#'   `msr_gdr_refit_error`, whose `audit` and `audit_file` elements hold the
#'   record and its path.
#' @details
#' # How each refit works
#'
#' At every call the function fits the model from its original estimates
#' (the start named `"original"`) and from each cold start, then keeps the
#' start with the highest finite log-likelihood. That start must be eligible:
#' converged (code 0) with no warnings and, for TMB models, usable covariance
#' and a largest absolute gradient below 0.001. If it is not eligible, the
#' call fails rather than returning a lower-likelihood fit in its place.
#'
#' When [multiScale_optim()] needs the fitted model rather than only its
#' likelihood, a final fit starts from the selected coefficients with the
#' original standard-error request. It must also be eligible and may lose no
#' more than 1e-5 log-likelihood units. C-engine intermediate fits skip
#' coefficient Hessians, which saves time.
#'
#' Only the columns named in `context$covs` may change between calls, and
#' only static site covariates qualify. A covariate whose name also appears in
#' the yearly or observation covariates is rejected, so give it a distinct
#' name.
#' Changes to the counts, the other covariates, the survey protocol or the
#' row order stop the call.
#'
#' Each TMB attempt builds a new TMB objective and starts its random effects
#' at zero. It always runs `TMB::sdreport()`, even with `se = FALSE`, so TMB
#' refits are not faster without standard errors. Usable covariance requires
#' a positive-definite Hessian, fixed covariance and joint precision. The
#' random effects must sit at a mode with positive-definite curvature and a
#' largest gradient and Newton decrement below 0.001. The objective,
#' recomputed at the optimizer's coordinates, must also agree with the
#' reported value within 1e-7. These are local checks. They
#' do not show that the mode is unique or that TMB's Laplace approximation is
#' accurate.
#'
#' # Declare TMB settings
#'
#' A C-engine fit keeps the optimizer settings it used, so they are read
#' automatically. A TMB fit in `unmarked` 1.5.2 does not keep them, and its
#' constructor does not pass optimizer controls through. Declare them in
#' `tmb_settings` instead. The refit function passes them to `optim()`
#' directly, without changing `unmarked`. A typical list is:
#'
#' ```
#' list(method = "BFGS", se = TRUE, threads = 1L,
#'      control = list(maxit = 2500L, reltol = 1e-10))
#' ```
#'
#' `random.start` and `inner.control` change how TMB searches for the random
#' effects. A different setting can reach a different numerical solution, so
#' treat a change as a new analysis and check its results. No setting is
#' changed automatically after a failure. With `record_memory = TRUE`, each
#' phase record stores the resident memory of the whole R process and its
#' lifetime peak, in bytes. These numbers include everything else loaded in
#' the session, so they are not the memory used by one phase.
#'
#' # Read the audit directory
#'
#' | File | Written | Contents |
#' |:---|:---|:---|
#' | `contract.rds` | When the refit function is created | Fingerprints of the original data and design, the settings, the starts, the original coefficients and optimizer summary, and the identity of `unmarked`, TMB and `multiScaleR`, with a SHA-256 of the whole contract |
#' | `refit-*-request.rds` | At the start of each call | The call's opening record |
#' | `refit-*-input.rds` | After the inputs pass their checks | Fingerprints of the new data and the changed columns |
#' | `refit-*-attempt-*-request.rds` | Before each start enters `unmarked` | The start and its name |
#' | `refit-*-attempt-*-point.rds` | TMB only, after optimization and before covariance | The optimizer result |
#' | `refit-*-attempt-*.rds` | After each start | The completed attempt record |
#' | `refit-*.rds` | At the end of each call | The summary, also attached to the returned model |
#'
#' A request without a matching completed record documents an attempt that
#' never finished, for example after an interruption. Keep it. Records contain
#' numbers, settings and SHA-256 fingerprints. They contain no raw data,
#' design matrices, Hessians or fitted models, so they stay small.
#'
#' The optimizer point file (class `msr_gdr_optimizer_point`) is written
#' before the covariance step, so a later covariance failure still leaves
#' evidence that the optimization finished. It holds the optimizer's `par`,
#' `value`, `counts`, `convergence` and `message`, the coordinate labels, data
#' and layout fingerprints, the cutoff `K`, the declared settings and the
#' elapsed time of each phase. It is always marked
#' `inference_eligible = FALSE`. It is not a fitted model, it gives no
#' standard errors, and a fit cannot be resumed from it.
#'
#' At each call the refit function also confirms that `unmarked`, TMB and
#' this function's own code have not changed since it was created. A change,
#' including a different installation in a new worker process, stops the call
#' before fitting.
#'
#' # Refit a response field
#'
#' For a model from [fit_gdr_response_field()], each call refits the field
#' with all of that function's checks, with the range held fixed. Use the
#' same `cold_starts` and the identical `tmb_settings`, and run
#' [multiScale_optim()] with `n_cores = 1`. The ordinary `update()` method
#' stops for field models, so this refit function is the way to refit one.
#'
#' # What the checks do not establish
#'
#' Coefficient standard errors from the returned model treat the selected
#' scales of effect as known. The function adds no priors, bootstrap
#' intervals or kernel-scale uncertainty. Agreement among the declared starts
#' does not guarantee that a better optimum does not exist elsewhere.
#'
#' The cutoff `K` (the largest abundance the likelihood considers at a point)
#' is unrelated to the `K` column of [aic_tab()], which counts parameters.
#' Reusing one cutoff does not replace a check that results are insensitive to
#' it.
#'
#' See the "Refit from several starting values" section of
#' `vignette("v12-distance-removal", package = "multiScaleR")` for a worked
#' example.
#' @seealso [multiScale_optim()], [gdr_start_values()],
#'   [fit_gdr_response_field()], [validate_unmarked_data()], [aic_tab()]
#' @examples
#' if (as.character(utils::packageVersion("unmarked")) == "1.5.2") {
#'   # Five point counts: birds by distance band and by time interval
#'   yd <- rbind(c(0, 0), c(1, 0), c(0, 2), c(2, 1), c(0, 0))
#'   yr <- rbind(c(0, 0, 0), c(0, 1, 0), c(1, 0, 1), c(1, 1, 1), c(0, 0, 0))
#'   frame <- unmarked::unmarkedFrameGDR(yd, yr, dist.breaks = c(0, 25, 75),
#'     unitsIn = "m", period.lengths = c(1, 2, 1),
#'     siteCovs = data.frame(hab = seq(-1, 1, length.out = 5)))
#'   # Fit the starting model (C engine) before declaring extra starts
#'   model <- unmarked::gdistremoval(lambdaformula = ~hab, data = frame,
#'     output = "abund", K = 40, se = FALSE, engine = "C",
#'     starts = c(log(8), 0, log(40), qlogis(0.3)), control = list(maxit = 2000))
#'   # One extra start, named like the model coefficients (link scales)
#'   cold <- unmarked::coef(model)
#'   cold[] <- 0
#'   cold["dist(Int)"] <- log(30)
#'   audit_dir <- tempfile("msr-gdr-audit-")
#'   refit <- make_gdr_refitter(model, list(cold = cold), audit_dir = audit_dir)
#'
#'   # multiScale_optim() calls refit() at each candidate scale. Calling it
#'   # once here, on the original data, shows what each call records.
#'   fit <- refit(model, frame, context = list(covs = "hab"))
#'   audit <- attr(fit, "msr_gdr_refit_audit")
#'   # Log-likelihood and eligibility of each start, then the selected start
#'   print(sapply(audit$attempts, function(a)
#'     c(logLik = a$logLik, eligible = a$eligible)))
#'   print(audit$attempts[[audit$selected]]$name)
#'   # contract.rds, plus request, input, attempt and summary records
#'   print(list.files(audit_dir))
#' }
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
      if (tmb) {
        score <- record$gradient
        maximum <- if (is.numeric(score) && length(score) == audit$parameter_count &&
          all(is.finite(score))) max(abs(score)) else Inf
        record$marginal_score <- list(maximum_absolute_gradient = maximum,
          threshold = 1e-3, eligible = is.finite(maximum) && maximum < 1e-3)
      }
      record$eligible <- !is.null(fit) && is.finite(record$logLik) && record$converged &&
        !length(record$warnings) && ((!se && !tmb) || identical(record$covariance$status, "available")) &&
        (!tmb || isTRUE(record$marginal_score$eligible))
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
