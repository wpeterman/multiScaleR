# TMB's saved objective does not retain evaluated optim controls. Require those
# values explicitly, and bind both the native input layout and provider bytes.
.msr_gdr_TMB_settings <- function(model, settings) {
  .msr_gdr_native_guard(model)
  if (!methods::is(model, "unmarkedFitGDR") || is.null(model@TMB) ||
      !is.environment(model@TMB$env) || !is.function(model@TMB$fn))
    stop("Supply an actual native TMB GDR model.", call. = FALSE)
  if (as.character(utils::packageVersion("unmarked")) != "1.5.2" ||
      as.character(utils::packageVersion("TMB")) != "1.9.25")
    stop("The TMB GDR adapter requires separate qualification for these provider versions.", call. = FALSE)
  if (!is.list(settings) || is.object(settings) || is.null(names(settings)) ||
      anyNA(names(settings)) || anyDuplicated(names(settings)) ||
      !all(c("method", "se", "threads", "control") %in% names(settings)) ||
      any(!names(settings) %in% c("method", "se", "threads", "control", "lower", "upper", "random.start", "inner.control", "record_memory")))
    stop("tmb_settings must explicitly contain method, se, threads and control.", call. = FALSE)
  if (!is.character(settings$method) || length(settings$method) != 1L ||
      is.na(settings$method) || !settings$method %in% c("BFGS", "Nelder-Mead", "CG", "L-BFGS-B") ||
      !is.logical(settings$se) || length(settings$se) != 1L || is.na(settings$se) ||
      !is.numeric(settings$threads) || length(settings$threads) != 1L ||
      is.na(settings$threads) || settings$threads != 1)
    stop("TMB refits require a deterministic native optimizer, scalar logical se and threads = 1.", call. = FALSE)
  control <- settings$control
  allowed <- c("maxit", "reltol", "abstol", "trace", "REPORT", "parscale", "ndeps",
    "alpha", "beta", "gamma", "type", "lmm", "factr", "pgtol", "warn.1d.NelderMead")
  if (!is.list(control) || is.object(control) || !length(control) ||
      is.null(names(control)) || anyNA(names(control)) || anyDuplicated(names(control)) ||
      any(!names(control) %in% allowed) || !all(vapply(control, function(x)
        is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) > 0L && all(is.finite(x)), logical(1))) ||
      is.null(control$maxit) || length(control$maxit) != 1L ||
      control$maxit <= 0 || control$maxit != floor(control$maxit) ||
      (!is.null(control$parscale) && any(control$parscale <= 0)) ||
      (!is.null(control$ndeps) && any(control$ndeps <= 0)))
    stop("TMB control must contain a positive integer maxit and supported finite numeric optim controls.", call. = FALSE)
  policy <- .msr_gdr_TMB_policy(settings)
  snapshot <- gdr_start_values(model)
  bounds <- settings[intersect(c("lower", "upper"), names(settings))]
  if (length(bounds) && (settings$method != "L-BFGS-B" || !all(vapply(bounds,
      function(x) is.numeric(x) && !is.object(x) && is.null(dim(x)) &&
        length(x) %in% c(1L, snapshot$coordinate_count) && !anyNA(x), logical(1)))))
    stop("TMB bounds require L-BFGS-B and scalar or optimizer-order numeric bounds.", call. = FALSE)
  if (!identical(model@nllFun, model@TMB$fn))
    stop("Stored TMB likelihood and native objective differ.", call. = FALSE)
  design <- get("getDesign", asNamespace("unmarked"))(model@data, model@formlist)
  native <- model@TMB$env$data
  prepared <- .msr_gdr_TMB_preparation(model)
  if (!isTRUE(all.equal(native, prepared, check.attributes = FALSE, tolerance = 0)))
    stop("The original TMB target differs from its stored frame, protocol or cutoff.", call. = FALSE)
  for (component in c("lambda", "phi", "dist", "rem")) {
    for (prefix in c("X", "Z")) {
      name <- paste0(prefix, "_", component)
      if (!isTRUE(all.equal(native[[name]], design[[name]], check.attributes = FALSE, tolerance = 0)))
        stop("The original TMB design differs from the current native preparation: ", name, ".", call. = FALSE)
    }
  }
  list(method = settings$method, se = settings$se, engine = "TMB", threads = 1L,
    dots = c(list(control = control), bounds), policy = policy, K = model@K,
    coefficient_names = names(unmarked::coef(model)),
    design = design[c("X_lambda", "X_phi", "X_dist", "X_rem")],
    constructor_names = paste0(snapshot$map$block, "[", snapshot$map$within_block, "]"),
    native_data_sha256 = .msr_gdr_hash(native),
    coordinate_map = snapshot$map, parameter_count = .msr_parameter_count(model),
    design_layout = .msr_gdr_TMB_layout(design))
}

.msr_gdr_TMB_layout <- function(design) {
  lapply(design[paste0(rep(c("X", "Z"), each = 4L), "_", c("lambda", "phi", "dist", "rem"))],
    function(x) list(dim = dim(x), dimnames_sha256 = .msr_gdr_hash(dimnames(x))))
}

# Stop the native constructor before tape creation. This reuses its actual data
# preparation, without evaluating an expired call or executing the saved tape.
.msr_gdr_TMB_preparation <- function(model) {
  constructor <- unmarked::gdistremoval
  local <- new.env(parent = environment(constructor))
  local$fit_TMB <- function(model, data, ...) stop(structure(
    list(message = "native preparation captured", call = NULL,
      data = c(list(model = model), data)),
    class = c("msr_gdr_preparation", "error", "condition")))
  environment(constructor) <- local
  formulas <- model@formlist
  captured <- tryCatch(constructor(lambdaformula = formulas$lambda,
    phiformula = formulas$phi, distanceformula = formulas$dist,
    removalformula = formulas$rem, data = model@data, K = model@K,
    keyfun = model@keyfun, mixture = model@mixture, output = model@output,
    unitsOut = model@unitsOut, engine = "TMB", threads = 1L),
    msr_gdr_preparation = identity)
  if (!inherits(captured, "msr_gdr_preparation"))
    stop("Could not capture the qualified native TMB preparation.", call. = FALSE)
  captured$data
}

.msr_gdr_TMB_provider <- function(files = TRUE) {
  native <- get("fit_TMB", asNamespace("unmarked"))
  out <- list(TMB_version = as.character(utils::packageVersion("TMB")),
    fit_TMB_sha256 = .msr_gdr_hash(list(formals = deparse(formals(native)), body = deparse(body(native)))))
  if (files) {
    out$TMB_R_sha256 <- digest::digest(file = file.path(find.package("TMB"), "R", "TMB.rdb"), algo = "sha256")
    out$TMB_DLL_sha256 <- digest::digest(file = getLoadedDLLs()[["TMB"]][["path"]], algo = "sha256")
    out$native_TMB_DLL_sha256 <- digest::digest(file = getLoadedDLLs()[["unmarked_TMBExports"]][["path"]], algo = "sha256")
  }
  out
}

# Keep covariance diagnostics compact. No covariance, joint precision or tape is
# serialized. Joint precision checks retain fixed/latent cross-component terms.
.msr_gdr_TMB_diagnostics <- function(result) {
  started <- proc.time()[["elapsed"]]
  sdr <- result$sdr
  covariance <- sdr$cov.fixed
  positive <- function(x) {
    if (is.null(x) || nrow(x) != ncol(x) || !nrow(x) || any(!is.finite(x))) return(FALSE)
    isTRUE(tryCatch({
      if (inherits(x, "sparseMatrix")) Matrix::Cholesky(x, LDL = FALSE) else chol(x)
      TRUE
    }, error = function(e) FALSE))
  }
  returned_value_error <- NULL
  returned_value <- tryCatch(result$TMB$fn(result$opt$par), error = function(e) {
    returned_value_error <<- conditionMessage(e); NA_real_
  })
  returned_value_gap <- abs(returned_value - result$opt$value)
  value_aligned <- length(returned_value_gap) == 1L && is.finite(returned_value_gap) && returned_value_gap < 1e-7
  gradient <- as.numeric(result$TMB$gr(result$opt$par))
  current <- result$TMB$env$last.par
  random <- result$TMB$env$random
  latent <- length(random)
  current_fixed <- setdiff(seq_along(current), random)
  current_fixed_aligned <- isTRUE(all.equal(as.numeric(current[current_fixed]),
    as.numeric(result$opt$par), tolerance = 1e-7))
  random_aligned <- !latent || isTRUE(all.equal(as.numeric(sdr$par.random),
    as.numeric(current[random]), tolerance = 1e-7))
  joint <- sdr$jointPrecision
  fixed_ok <- positive(covariance) && nrow(covariance) == length(result$opt$par)
  joint_ok <- if (latent) positive(joint) &&
    nrow(joint) == length(result$opt$par) + latent else TRUE
  conditional_ok <- TRUE; maximum <- decrement <- 0
  if (latent) {
    # Always evaluate the current complete point, never spHess's initial default.
    precision <- result$TMB$env$spHess(current, random = TRUE)
    conditional_ok <- positive(precision) && nrow(precision) == latent
    stationarity <- tryCatch({
      score <- as.numeric(result$TMB$env$f(current, order = 1L))
      if (length(score) != length(current) || any(!is.finite(score)))
        stop("The direct native joint gradient is unavailable or nonfinite.")
      score <- score[random]
      list(maximum = max(abs(score)), decrement = if (conditional_ok)
        sqrt(max(0, sum(score * as.numeric(solve(precision, score))))) else NA_real_)
    }, error = function(e) list(maximum = NA_real_, decrement = NA_real_, reason = conditionMessage(e)))
    maximum <- stationarity$maximum; decrement <- stationarity$decrement
  }
  stationary <- is.finite(maximum) && is.finite(decrement) && maximum < 1e-3 && decrement < 1e-3
  aligned <- isTRUE(all.equal(as.numeric(sdr$par.fixed), as.numeric(result$opt$par), tolerance = 1e-7))
  workflow <- result$workflow
  if (is.null(workflow)) workflow <- list(policy = .msr_gdr_TMB_policy(), phases = list())
  workflow$phases$qualification_complete <- .msr_gdr_TMB_phase("qualification_complete", started, workflow$policy)
  list(requested = TRUE, status = if (isTRUE(sdr$pdHess) && fixed_ok && joint_ok && conditional_ok && aligned && current_fixed_aligned && random_aligned && value_aligned && stationary &&
      length(gradient) == length(result$opt$par) && all(is.finite(gradient))) "available" else "unusable",
    pdHess = isTRUE(sdr$pdHess), fixed_covariance_positive_definite = fixed_ok,
    joint_precision_positive_definite = joint_ok, covariance_at_optimizer = aligned,
    covariance_at_current_latent_mode = random_aligned, returned_value_gap = returned_value_gap,
    current_fixed_coordinates_at_optimizer = current_fixed_aligned,
    conditional_latent_precision_positive_definite = conditional_ok,
    conditional_curvature_point = if (latent) "native_current_mode" else "no_latent_coordinates",
    conditional_latent_stationary = stationary, conditional_latent_maximum_absolute_gradient = maximum,
    conditional_latent_Newton_decrement = decrement, conditional_stationarity_threshold = 1e-3,
    conditional_stationarity_error = if (latent) stationarity$reason else NULL,
    returned_value_error = returned_value_error,
    conditional_stationarity_scope = "local_native_mode_not_global_uniqueness_or_integral_accuracy",
    fixed_coordinates = length(result$opt$par), latent_coordinates = latent,
    gradient = gradient, maximum_absolute_gradient = max(abs(gradient)),
    sdreport_computed = TRUE, coefficient_hessians_omitted = FALSE, workflow = workflow,
    uncertainty = "native_joint_fixed_and_latent_conditional_on_spatial_scales")
}

.msr_gdr_TMB_refit <- function(args, settings, checkpoint = NULL) {
  constructor <- unmarked::gdistremoval
  local <- new.env(parent = environment(constructor))
  diagnostic <- NULL
  optimizer <- NULL
  # gdistremoval 1.5.2 does not forward ... to fit_TMB. This local binding
  # forwards the explicitly evaluated controls without patching its namespace.
  local$fit_TMB <- function(model, data, params, random, starts, method, ...) {
    old_threads <- TMB::openmp(DLL = "unmarked_TMBExports")
    on.exit(TMB::openmp(old_threads, DLL = "unmarked_TMBExports"))
    TMB::openmp(1L, DLL = "unmarked_TMBExports")
    result <- .msr_gdr_TMB_fit(model, data, params, random, starts, method, settings$dots,
      policy = settings$policy, checkpoint = checkpoint)
    optimizer <<- .msr_gdr_optimizer_checkpoint(result$opt)
    diagnostic <<- .msr_gdr_TMB_diagnostics(result)
    result
  }
  environment(constructor) <- local
  value <- tryCatch(do.call(constructor, args), error = function(e) {
    if (inherits(e, c("msr_gdr_tmb_error", "msr_gdr_audit_error"))) stop(e)
    stop(structure(list(message = conditionMessage(e), call = NULL,
      optimizer = optimizer, covariance = diagnostic),
      class = c("msr_gdr_tmb_error", "error", "condition")))
  })
  # Do not retain the temporary constructor closure and its call environment.
  value@call[[1L]] <- quote(unmarked::gdistremoval)
  for (name in names(settings$dots)) value@call[[name]] <- settings$dots[[name]]
  attr(value, "msr_gdr_tmb_diagnostics") <- diagnostic
  value
}

# The qualified native fitting sequence, with evaluated controls and an explicit
# sdreport point. Native fit_TMB otherwise uses last.par.best, whose values can
# differ from opt$par, leaving returned estimates and optimizer metadata unequal.
# The public constructor supplies zero latent starts; reference diagnostics may use warm states.
.msr_gdr_TMB_fit <- function(model, data, params, random, starts, method, dots,
    policy = NULL, checkpoint = NULL) {
  if (is.null(policy)) policy <- .msr_gdr_TMB_policy()
  if (!is.null(checkpoint) && !is.function(checkpoint)) stop("Use an internal checkpoint callback.")
  started <- proc.time()[["elapsed"]]
  phases <- list(construction_started = .msr_gdr_TMB_phase("construction_started", started, policy))
  fixed <- names(params)[!names(params) %in% random]
  template <- params[fixed]
  count <- length(unlist(template))
  if (length(starts) != count)
    stop("The number of native TMB constructor starts changed.", call. = FALSE)
  params <- replace(params, fixed, utils::relist(unname(starts), template))
  objective <- .msr_gdr_TMB_adfun(model, data, params, random, policy)
  objective$starts_order <- rep(fixed, lengths(template))
  phases$construction_complete <- .msr_gdr_TMB_phase("construction_complete", started, policy)
  optimization_started <- proc.time()[["elapsed"]]
  opt <- do.call(stats::optim, c(list(par = objective$par, fn = objective$fn,
    gr = objective$gr, method = method), dots))
  phases$optimization_complete <- .msr_gdr_TMB_phase("optimization_complete", optimization_started, policy)
  # This write precedes expensive covariance and even the returned-value check.
  # A numeric point is retained evidence, never a fitted S4 or inference result.
  point_file <- if (is.null(checkpoint)) NULL else
    checkpoint(.msr_gdr_TMB_point(opt, model, data, params, random, policy, phases))
  covariance_started <- covariance_request_file <- NULL
  sdr <- tryCatch({
    value_started <- proc.time()[["elapsed"]]
    value <- objective$fn(opt$par)
    if (!is.finite(value) || !is.finite(opt$value) || abs(value - opt$value) > 1e-7)
      stop("Native TMB likelihood differs at the returned optimizer coordinates.", call. = FALSE)
    phases$value_qualification_complete <- .msr_gdr_TMB_phase("value_qualification_complete", value_started, policy)
    covariance_started <- proc.time()[["elapsed"]]
    phases$covariance_started <- .msr_gdr_TMB_phase("covariance_started", covariance_started, policy)
    if (!is.null(point_file)) {
      covariance_request_file <- sub("[.]rds$", "-covariance-started.rds", point_file)
      .msr_gdr_atomic_audit(list(schema_version = 1L, type = "GDR_covariance_started",
        inference_eligible = FALSE, point_file = point_file,
        point_sha256 = digest::digest(file = point_file, algo = "sha256"),
        phase = phases$covariance_started), covariance_request_file)
    }
    value <- TMB::sdreport(objective, par.fixed = opt$par, getJointPrecision = TRUE)
    phases$covariance_complete <- .msr_gdr_TMB_phase("covariance_complete", covariance_started, policy)
    value
  }, error = function(e) {
    if (inherits(e, "msr_gdr_audit_error")) stop(e)
    if (!is.null(covariance_started)) phases$covariance_failed <-
      .msr_gdr_TMB_phase("covariance_failed", covariance_started, policy)
    stop(structure(list(message = conditionMessage(e), call = NULL,
      optimizer = .msr_gdr_optimizer_checkpoint(opt), covariance = NULL, point_file = point_file,
      workflow = list(policy = policy, phases = phases, covariance_request_file = covariance_request_file)),
      class = c("msr_gdr_tmb_error", "error", "condition")))
  })
  sdr$par <- objective$par
  list(opt = opt, TMB = objective, sdr = sdr, AIC = 2 * opt$value + 2 * count,
    workflow = list(policy = policy, phases = phases, point_file = point_file,
      covariance_request_file = covariance_request_file))
}

.msr_gdr_refit_start <- function(model, settings) {
  if (settings$engine == "C") return(unmarked::coef(model))
  snapshot <- gdr_start_values(model)
  stats::setNames(snapshot$starts, paste0(snapshot$map$block, "[", snapshot$map$within_block, "]"))
}
