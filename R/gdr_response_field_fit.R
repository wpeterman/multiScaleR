.msr_gdr_field_source <- function() {
  functions <- c(".msr_gdr_field_source", "fit_gdr_response_field", ".msr_gdr_field_fit", ".msr_gdr_field_guard",
    ".msr_gdr_field_predict", ".msr_gdr_field_prediction_function",
    ".msr_gdr_field_refitter", ".msr_gdr_field_map", ".msr_gdr_field_geometry_hash", ".msr_gdr_field_provider",
    ".msr_gdr_response_field_context", ".msr_gdr_gaussian_data",
    ".msr_unmarked_prediction_basis", ".msr_unmarked_fixed_design",
    ".msr_unmarked_curve_basis", ".msr_cv_prediction_terms",
    ".msr_unmarked_random_design", ".msr_response_field_crs",
    ".msr_gdr_TMB_fit", ".msr_gdr_TMB_diagnostics", ".msr_gdr_field_configuration_hash",
    ".msr_gdr_field_binding", ".msr_gdr_field_geometry", ".msr_gdr_field_contract",
    ".msr_gdr_field_prediction_budget", ".msr_gdr_field_raster_eval",
    ".add_site_covariate_rasters")
  # Runtime byte compilation is not a source edit. Bind canonical formals
  # and code, as in the existing GDR adapter identity, rather than bytecode.
  ordinary <- lapply(functions, function(name) {
    fn <- get(name, environment(.msr_gdr_field_source), inherits = FALSE)
    list(name, deparse(formals(fn)), deparse(body(fn)))
  })
  registered <- lapply(c("predict", "vcov", "update", "simulate", "profile", "ranef", "confint"),
    function(name) {
      fn <- methods::selectMethod(name, "unmarkedFitGDRField",
        fdef = methods::getGeneric(name, where = asNamespace("unmarked")))
      list(name, deparse(formals(fn)), deparse(body(fn)))
    })
  .msr_gdr_hash(list(ordinary, registered))
}

.msr_gdr_field_provider <- function() {
  value <- .msr_gdr_TMB_provider()
  value$unmarked_R_sha256 <- digest::digest(file = file.path(find.package("unmarked"),
    "R", "unmarked.rdb"), algo = "sha256")
  value$constructor_sha256 <- .msr_gdr_hash(list(formals(unmarked::gdistremoval),
    body(unmarked::gdistremoval)))
  value
}

.msr_gdr_field_map <- function(Z, rows) {
  if (!ncol(Z)) return(integer())
  z <- Matrix::summary(Z)
  if (nrow(z) != rows || any(z$x != 1) || anyDuplicated(z$i) ||
      !setequal(z$i, seq_len(rows)))
    stop("The GDR field requires one original random intercept per retained row.", call. = FALSE)
  as.integer(z$j[match(seq_len(rows), z$i)])
}

#' Fit a Fixed-Range Spatial Abundance Field to GDR
#'
#' Replace one native iid abundance intercept with an exponential Gaussian
#' response field while retaining the original distance/removal observation
#' likelihood. Correlation range is supplied, not estimated. Ecological kernel
#' scales remain separate parameters in the usual fixed-reference machinery.
#' @param model An actual native TMB `unmarkedFitGDR` used as a template. The
#'   qualified domain is density per square kilometer, Poisson abundance,
#'   half-normal distance detection, one
#'   primary period, complete canonical 75-meter point-circle bins, one
#'   abundance random intercept and optionally one independent observer
#'   intercept in distance detection. Existing offset and other native guards
#'   apply. Each abundance state must have an observed location.
#' @param coords Numeric matrix or data frame with two projected coordinates
#'   per retained survey. Row names must equal the original site identities.
#'   Repeated surveys for one abundance state must share exactly one location.
#' @param crs Known projected two-dimensional CRS, supplied as an EPSG code or
#'   another input accepted by `sf::st_crs()`. No coordinate conversion occurs.
#' @param range Positive fixed exponential correlation distance in CRS units.
#'   It adds no estimated parameter. This is a response-correlation distance,
#'   not an ecological kernel scale.
#' @param cold_starts Nonempty named list of independent finite coefficient/SD
#'   starts in the constructor order returned by [gdr_start_values()]. Original
#'   fitted starts and zero latent states are also checked. Failed attempts and
#'   disagreements stop rather than returning a preferred passing fit.
#' @param audit_dir Directory for compact atomic request/attempt/result records.
#'   Records retain failures, controls and identities, not models or tapes.
#' @param tmb_settings Explicit native settings as in [make_gdr_refitter()].
#'   This initial field route requires BFGS, `se = TRUE`, one thread and positive
#'   `control$maxit`. Full fixed and joint curvature remain mandatory.
#' @param max_memory_mb Positive finite conservative planning budget in MiB.
#'   Default 512. Dense field/covariance and finite-K tape planning must fit it.
#'   This estimate is not an operating-system memory cap.
#' @return An `unmarkedFitGDRField` with the native maximized Laplace likelihood,
#'   coefficient/SD estimates and identified field metadata in `response_field`.
#'   Abundance random coefficients and their component covariance are expressed
#'   as physical field values; native TMB uses equivalent correlation-whitened
#'   coordinates internally. Fixed range is not included in the AIC penalty.
#'   Use `predict(fit, type = "lambda", ...)` for conditional density predictions
#'   and [make_gdr_refitter()] for fixed-reference ecological-scale refits.
#' @details
#' The finite integration cutoff, sampled area, distance/removal constants and
#'   normalized Gaussian SD are unchanged. No year field or population prior is
#'   introduced. Both declared starts must return warning-free converged fits
#'   with usable full covariance, positive current latent precision and a small
#'   marginal score. This is an opt-in native fitting route; private RTMB fitting
#'   admission and all unrelated guards remain unchanged.
#'
#' Predictions retain full fitted coefficient/field cross-covariance. At new
#'   locations, conditional field residual variance is identified separately and
#'   included in the latent predictor SE. Back-transformed predictions are the
#'   inverse-link conditional predictor, not an average over a lognormal field.
#'   No observation noise or ecological-scale uncertainty is added. Setting
#'   `re.form = NA` supplies a zero-effect reference, not field integration.
#'
#' Ecological-scale fitting uses the refit backend with n_cores = 1 and
#' PSOCK = FALSE. Native fit/tape objects are not sent to fit workers.
#'
#' Estimated ranges, other GDR mixtures/keys, annual availability, spatial CV,
#'   field simulation and profile intervals require separate qualification.
#'   Ordinary update is guarded; retain the audited refitter. See the inference
#'   contracts vignette for the supported scope and uncertainty interpretation.
#' @seealso [make_gdr_refitter()], [gdr_start_values()], [response_field_geometry()]
#' @examples
#' # coords must be in exactly the retained survey order, with a known CRS.
#' # Fit the ordinary native GDR template first, then supply its original
#' # independent starts and explicit native controls to this opt-in constructor.
#' @md
#' @export
fit_gdr_response_field <- function(model, coords, crs, range, cold_starts,
    audit_dir, tmb_settings, max_memory_mb = 512) {
  if (methods::is(model, "unmarkedFitGDRField"))
    stop("Use make_gdr_refitter() to retain an existing response field.", call. = FALSE)
  .msr_gdr_selection_data(model, seq_len(nrow(model@data@yDistance)))
  settings <- .msr_gdr_TMB_settings(model, tmb_settings)
  if (!identical(settings$method, "BFGS") || !isTRUE(settings$se))
    stop("The qualified GDR field route requires BFGS and se = TRUE.", call. = FALSE)
  native <- .msr_gdr_TMB_preparation(model)
  pars <- model@TMB$env$parameters
  rows <- nrow(native$X_lambda)
  if (!identical(model@mixture, "P") || !identical(model@keyfun, "halfnorm") ||
      !identical(model@output, "density") || !identical(model@unitsOut, "kmsq") ||
      model@data@numPrimary != 1L || length(model@sitesRemoved) ||
      !identical(as.numeric(model@data@dist.breaks), as.numeric(native$db)))
    stop("The GDR field requires complete single-period Poisson/half-normal rows with output = 'density' and unitsOut = 'kmsq'.", call. = FALSE)
  ids <- .msr_unmarked_observations(model)$ids
  numeric <- list(Xlambda = native$X_lambda, Xdist = native$X_dist, Xrem = native$X_rem,
    yd = unname(model@data@yDistance), yr = unname(model@data@yRemoval),
    block = .msr_gdr_field_map(native$Z_lambda, rows),
    observer = .msr_gdr_field_map(native$Z_dist, rows),
    periods = as.numeric(model@data@period.lengths), breaks = as.numeric(native$db),
    area = native$A, K = as.numeric(model@K), annulus_area_m2 = native$a,
    annulus_weights = native$u, row_ids = list(survey = ids,
      removal_survey = rep(ids, each = ncol(model@data@yRemoval)),
      removal_interval = rep(seq_len(ncol(model@data@yRemoval)), rows)))
  context <- .msr_gdr_response_field_context(numeric, pars, coords, crs, range,
    max_memory_mb, .geometry_only = TRUE)
  bytes <- 8 * (64 * as.double(length(pars$b_lambda) + length(pars$b_dist) +
    length(model@opt$par))^2 + 32 * rows * length(pars$b_lambda)) +
    512 * as.double(rows) * (as.double(model@K) + 1)
  if (!is.finite(bytes) || bytes > max_memory_mb * 2^20)
    stop("Native GDR field covariance/tape planning exceeds max_memory_mb.", call. = FALSE)
  components <- names(.msr_gdr_active_formulas(model))
  bases <- stats::setNames(lapply(components, function(type)
    .msr_unmarked_prediction_basis(model, type)), components)
  raw <- stats::setNames(lapply(components, function(type)
    get("get_covariates", asNamespace("unmarked"))(model, type)), components)
  for (type in components) {
    bases[[type]]$terms <- .msr_cv_prediction_terms(bases[[type]]$terms,
      raw[[type]], qualify = TRUE, allow_learning = FALSE)
    environment(bases[[type]]$terms) <- baseenv()
    original <- native[[paste0("X_", type)]]
    if (!isTRUE(all.equal(.msr_unmarked_fixed_design(bases[[type]], raw[[type]])$X,
        original, tolerance = 1e-12, check.attributes = FALSE)))
      stop("Isolated field prediction recipes must reproduce every original design.", call. = FALSE)
  }
  forms <- lapply(model@formlist, function(f) { environment(f) <- baseenv(); f })
  # State order follows the original native random design, not coordinate sort.
  locations <- unname(as.matrix(coords))[match(seq_along(pars$b_lambda), numeric$block), , drop = FALSE]
  lower <- t(chol(exp(-as.matrix(stats::dist(locations))/range)))
  field <- list(definition = context$field, project = context$project,
    field_id = context$field_id, coords = unname(as.matrix(coords)), crs = context$field$crs,
    row_ids = ids, block = numeric$block, lower = lower, bases = bases, raw = raw, forms = forms,
    original_Z = native[c("Z_lambda", "Z_dist")], settings = settings,
    provider = .msr_gdr_field_provider(), source = .msr_gdr_field_source(),
    max_memory_mb = max_memory_mb, planned_bytes = bytes, format_version = 1L)
  output <- .msr_gdr_field_fit(model, model@data, cold_starts, audit_dir, field)
  output@call <- match.call()
  output
}

.msr_gdr_field_fit <- function(template, frame, cold_starts, audit_dir, field,
    .objective_only = FALSE, .provenance = list()) {
  settings <- field$settings
  snapshot <- gdr_start_values(template); starts <- snapshot$starts
  labels <- paste0(snapshot$map$block, "[", snapshot$map$within_block, "]")
  if (!is.list(cold_starts) || !length(cold_starts) || is.null(names(cold_starts)) ||
      anyNA(names(cold_starts)) || any(!nzchar(names(cold_starts))) ||
      anyDuplicated(names(cold_starts)) || any(names(cold_starts) %in% c("original", "final")))
    stop("Supply distinct named independent cold coefficient/SD starts.", call. = FALSE)
  cold_starts <- lapply(cold_starts, function(p) {
    if (!is.numeric(p) || !is.null(dim(p)) || length(p) != length(starts) || any(!is.finite(p)))
      stop("Every cold start must match the finite original constructor coordinates.", call. = FALSE)
    if (!is.null(names(p))) {
      if (anyNA(names(p)) || anyDuplicated(names(p)) || !setequal(names(p), labels))
        stop("Named field starts must use exact constructor block/position names.", call. = FALSE)
      p <- p[labels]
    }
    if (identical(as.numeric(p), as.numeric(starts)))
      stop("A cold field start must be independent of original fitted coordinates.", call. = FALSE)
    unname(p)
  })
  if (anyDuplicated(vapply(cold_starts, .msr_gdr_hash, character(1))))
    stop("Cold field starts must differ from one another.", call. = FALSE)
  if (!is.character(audit_dir) || length(audit_dir) != 1L || is.na(audit_dir) || !nzchar(audit_dir))
    stop("audit_dir must name one explicit directory.", call. = FALSE)
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  audit_dir <- normalizePath(audit_dir, winslash = "/", mustWork = TRUE)
  record <- list(status = "started", target = "native_P_finite_K_fixed_range_exponential_Gaussian",
    field_id = field$field_id, source = field$source, provider = field$provider,
    controls = settings[c("method", "se", "threads", "dots")], attempts = list(),
    K = as.numeric(template@K), parameter_count = snapshot$coordinate_count,
    range_estimated = FALSE, annual_field_coordinates = 0L,
    objective_only = .objective_only, provenance = .provenance,
    input_sha256 = .msr_gdr_hash(frame))
  audit <- tempfile("field-fit-", tmpdir = audit_dir, fileext = ".rds")
  record$path <- audit
  on.exit({
    if (identical(record$status, "started")) record$status <- "failed_before_return"
    .msr_gdr_atomic_audit(record, audit)
  })
  .msr_gdr_atomic_audit(record, sub("[.]rds$", "-request.rds", audit))
  changed <- template; changed@data <- frame
  field$raw <- stats::setNames(lapply(names(field$bases), function(type)
    get("get_covariates", asNamespace("unmarked"))(changed, type)), names(field$bases))
  matrices <- lapply(names(field$bases), function(type) {
    x <- .msr_unmarked_fixed_design(field$bases[[type]], field$raw[[type]])
    if (!is.null(x$offset) && any(x$offset != 0))
      stop("Active GDR offsets remain unqualified.", call. = FALSE)
    x$X
  }); names(matrices) <- paste0("X_", names(field$bases))
  record$input_design_sha256 <- .msr_gdr_hash(matrices)
  .msr_gdr_atomic_audit(record, sub("[.]rds$", "-input.rds", audit))
  run <- function(start, name, index) {
    attempt <- list(name = name, start = unname(start), status = "started", warnings = character(), error = NULL)
    path <- sub("[.]rds$", paste0("-attempt-", index, ".rds"), audit)
    on.exit(.msr_gdr_atomic_audit(attempt, path))
    constructor <- unmarked::gdistremoval
    scope <- new.env(parent = environment(constructor)); result <- NULL
    scope$fit_TMB <- function(model, data, params, random, starts, method, ...) {
      if (!isTRUE(all.equal(data$Z_lambda, field$original_Z$Z_lambda, tolerance = 0)) ||
          !isTRUE(all.equal(data$Z_dist, field$original_Z$Z_dist, tolerance = 0)))
        stop("A GDR field refit changed original random-state rows or columns.", call. = FALSE)
      for (n in names(matrices)) {
        if (!identical(dim(data[[n]]), dim(matrices[[n]])))
          stop("A field refit changed original fixed-design dimensions.", call. = FALSE)
        data[[n]] <- matrices[[n]]
      }
      data$Z_lambda <- Matrix::Matrix(data$Z_lambda %*% field$lower, sparse = TRUE)
      old_threads <- TMB::openmp(DLL = "unmarked_TMBExports")
      on.exit(TMB::openmp(old_threads, DLL = "unmarked_TMBExports"))
      TMB::openmp(1L, DLL = "unmarked_TMBExports")
      result <<- .msr_gdr_TMB_fit(model, data, params, random, starts, method, settings$dots)
      result
    }
    environment(constructor) <- scope
    forms <- template@formlist
    fit <- tryCatch(withCallingHandlers(constructor(lambdaformula = forms$lambda,
      phiformula = forms$phi, distanceformula = forms$dist, removalformula = forms$rem,
      data = frame, K = template@K, keyfun = template@keyfun, mixture = template@mixture,
      output = template@output, unitsOut = template@unitsOut, engine = "TMB", threads = 1L,
      starts = unname(start), method = settings$method),
      warning = function(w) attempt$warnings <<- c(attempt$warnings, conditionMessage(w))),
      error = function(e) {
        attempt$status <<- "failed"; attempt$error <<- conditionMessage(e)
        if (inherits(e, "msr_gdr_tmb_error")) {
          attempt$optimizer <<- e$optimizer; attempt$covariance <<- e$covariance
        }
        NULL
      })
    if (is.null(fit)) return(list(fit = NULL, attempt = attempt))
    attempt$optimizer <- .msr_gdr_optimizer_checkpoint(result$opt)
    attempt$logLik <- -result$opt$value
    attempt$convergence <- result$opt$convergence
    attempt$converged <- identical(as.integer(result$opt$convergence), 0L)
    attempt$SE_requested <- TRUE
    attempt$covariance <- .msr_gdr_TMB_diagnostics(result)
    attempt$gradient <- attempt$covariance$gradient
    attempt$eligible <- identical(as.integer(result$opt$convergence), 0L) &&
      !length(attempt$warnings) && identical(attempt$covariance$status, "available") &&
      attempt$covariance$maximum_absolute_gradient < 1e-3
    attempt$status <- "returned"
    list(fit = fit, result = result, attempt = attempt)
  }
  declared <- c(list(original = starts), cold_starts)
  values <- rep(NA_real_, length(declared)); first <- NULL; difference <- 0
  eligible <- logical(length(declared)); selected <- NULL; selected_index <- NA_integer_
  for (i in seq_along(declared)) {
    value <- run(declared[[i]], names(declared)[i], i); record$attempts[[i]] <- value$attempt
    eligible[i] <- isTRUE(value$attempt$eligible)
    if (!is.null(value$result)) {
      values[i] <- value$result$opt$value
      if (is.null(first)) first <- value$result$opt$par
      difference <- max(difference, abs(value$result$opt$par - first))
      if (is.null(selected) || values[i] < selected$result$opt$value) {
        selected <- value; selected_index <- i
      }
    }
    .msr_gdr_atomic_audit(record, sub("[.]rds$", paste0("-checkpoint-", i, ".rds"), audit))
    rm(value); gc(verbose = FALSE)
  }
  if (!all(eligible)) {
    record$status <- "failed"
    stop("Every declared GDR field start must pass the unchanged fit and curvature gates; inspect audit_dir.", call. = FALSE)
  }
  if (diff(range(values)) > 1e-5 || difference > 1e-3) {
    record$status <- "start_disagreement"
    stop("Declared GDR field starts disagree; retain the failed reference gate.", call. = FALSE)
  }
  fit <- selected$fit; result <- selected$result
  result$TMB$gr(result$opt$par); point <- result$TMB$env$last.par
  joint_covariance <- as.matrix(Matrix::solve(result$sdr$jointPrecision))
  fixed <- names(point) == "beta_lambda"; latent <- names(point) == "b_lambda"
  indices <- c(which(fixed), which(latent))
  transform <- as.matrix(Matrix::bdiag(diag(sum(fixed)), field$lower))
  estimate <- fit@estimates@estimates$lambda
  estimate@estimates[] <- drop(transform %*% point[indices])
  covariance_names <- dimnames(estimate@covMat)
  estimate@covMat <- transform %*% joint_covariance[indices, indices, drop = FALSE] %*% t(transform)
  dimnames(estimate@covMat) <- covariance_names
  fit@estimates@estimates$lambda <- estimate
  field$point <- point; field$joint_covariance <- joint_covariance
  field$diagnostics <- selected$attempt$covariance; field$audit_file <- audit
  field$configuration_sha256 <- .msr_gdr_field_configuration_hash(field)
  field$native_data_sha256 <- .msr_gdr_hash(fit@TMB$env$data)
  field$response_sha256 <- .msr_gdr_hash(frame)
  field$prediction_sha256 <- .msr_gdr_hash(list(point, joint_covariance, field$bases, field$raw, field$forms, field$field_id))
  field$geometry_sha256 <- .msr_gdr_field_geometry_hash(field)
  field$summary_sha256 <- .msr_gdr_hash(list(fit@estimates, fit@opt$par,
    fit@opt$value, fit@opt$convergence, fit@negLogLike, fit@AIC, lapply(fit@formlist, deparse),
    fit@K, fit@mixture, fit@keyfun, fit@output, fit@unitsOut))
  record$selected <- selected_index
  record$maximum_start_value_difference <- diff(range(values))
  fit@call[[1L]] <- quote(multiScaleR::fit_gdr_response_field)
  attr(fit, "msr_gdr_tmb_diagnostics") <- selected$attempt$covariance
  output <- methods::new("unmarkedFitGDRField", fit, response_field = field)
  tryCatch(.msr_gdr_field_guard(output), error = function(e) {
    record$error <<- conditionMessage(e); stop(e)
  })
  record$status <- "returned"
  attr(output, "msr_gdr_refit_audit") <- record
  output
}

.msr_gdr_field_guard <- function(model) {
  field <- model@response_field
  .msr_gdr_native_guard(model)
  checks <- c(
    source = identical(field$source, .msr_gdr_field_source()),
    provider = identical(field$provider, .msr_gdr_field_provider()),
    likelihood = identical(model@nllFun, model@TMB$fn),
    configuration = identical(field$configuration_sha256, .msr_gdr_field_configuration_hash(field)),
    response = identical(field$response_sha256, .msr_gdr_hash(model@data)),
    native_data = identical(field$native_data_sha256, .msr_gdr_hash(model@TMB$env$data)),
    geometry = identical(field$geometry_sha256, .msr_gdr_field_geometry_hash(field)),
    summary = identical(field$summary_sha256, .msr_gdr_hash(list(model@estimates,
      model@opt$par, model@opt$value, model@opt$convergence, model@negLogLike,
      model@AIC, lapply(model@formlist, deparse), model@K, model@mixture, model@keyfun,
      model@output, model@unitsOut))),
    prediction = identical(field$prediction_sha256,
      .msr_gdr_hash(list(field$point, field$joint_covariance, field$bases, field$raw, field$forms, field$field_id))))
  if (!all(checks))
    stop("The retained GDR response-field provider, target or prediction context changed: ",
      paste(names(checks)[!checks], collapse = ", "), ".", call. = FALSE)
  invisible(field)
}

# Bind every captured numeric projection input and both function definitions.
# A first prediction can compile a closure in place; its code and inputs are
# unchanged. Serializing that closure would incorrectly treat this as drift.
.msr_gdr_field_geometry_hash <- function(field) {
  project <- field$project; scope <- environment(project)
  expected <- c("m", "n", "max_memory_mb", "declared", "locations", "range", "lower", "field_id", "budget")
  if (!is.function(project) || !is.environment(scope) ||
      !identical(parent.env(scope), environment(.msr_gdr_response_field_context)) ||
      !setequal(ls(scope, all.names = TRUE), expected))
    stop("The retained GDR field projection context changed.", call. = FALSE)
  budget <- get("budget", scope, inherits = FALSE)
  if (!is.function(budget) || !identical(environment(budget), scope))
    stop("The retained GDR field projection budget context changed.", call. = FALSE)
  .msr_gdr_hash(list(field$definition, field$coords, field$block, field$lower,
    field$max_memory_mb, deparse(formals(project)), deparse(body(project)),
    deparse(formals(budget)), deparse(body(budget)),
    mget(setdiff(expected, "budget"), scope, inherits = FALSE)))
}

# Controls, native layout and diagnostics are part of this exact fitted target.
# They cannot be replaced on an existing object to admit a different refitter.
.msr_gdr_field_configuration_hash <- function(field) {
  .msr_gdr_hash(field[c("settings", "original_Z", "row_ids", "crs", "planned_bytes",
    "format_version", "diagnostics")])
}
