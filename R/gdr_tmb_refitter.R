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
      any(!names(settings) %in% c("method", "se", "threads", "control", "lower", "upper")))
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
    dots = c(list(control = control), bounds), K = model@K,
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
  sdr <- result$sdr
  covariance <- sdr$cov.fixed
  positive <- function(x) {
    if (is.null(x) || nrow(x) != ncol(x) || !nrow(x) || any(!is.finite(x))) return(FALSE)
    isTRUE(tryCatch({
      if (inherits(x, "sparseMatrix")) Matrix::Cholesky(x, LDL = FALSE) else chol(x)
      TRUE
    }, error = function(e) FALSE))
  }
  gradient <- as.numeric(result$TMB$gr(result$opt$par))
  latent <- length(result$TMB$env$random)
  joint <- sdr$jointPrecision
  fixed_ok <- positive(covariance) && nrow(covariance) == length(result$opt$par)
  joint_ok <- if (latent) positive(joint) &&
    nrow(joint) == length(result$opt$par) + latent else TRUE
  conditional_ok <- if (latent) {
    # spHess() otherwise defaults to the tape's INITIAL joint coordinates.
    # gr(opt$par) above updates last.par to the actual conditional mode, as
    # used by native sdreport. Always supply that complete point explicitly.
    precision <- result$TMB$env$spHess(result$TMB$env$last.par, random = TRUE)
    positive(precision) && nrow(precision) == latent
  } else TRUE
  aligned <- isTRUE(all.equal(as.numeric(sdr$par.fixed), as.numeric(result$opt$par),
    tolerance = 1e-7))
  list(requested = TRUE, status = if (isTRUE(sdr$pdHess) && fixed_ok && joint_ok && conditional_ok && aligned &&
      length(gradient) == length(result$opt$par) && all(is.finite(gradient))) "available" else "unusable",
    pdHess = isTRUE(sdr$pdHess), fixed_covariance_positive_definite = fixed_ok,
    joint_precision_positive_definite = joint_ok, covariance_at_optimizer = aligned,
    conditional_latent_precision_positive_definite = conditional_ok,
    conditional_curvature_point = if (latent) "native_current_mode" else "no_latent_coordinates",
    fixed_coordinates = length(result$opt$par), latent_coordinates = latent,
    gradient = gradient, maximum_absolute_gradient = max(abs(gradient)),
    sdreport_computed = TRUE, coefficient_hessians_omitted = FALSE,
    uncertainty = "native_joint_fixed_and_latent_conditional_on_spatial_scales")
}

.msr_gdr_TMB_refit <- function(args, settings) {
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
    result <- .msr_gdr_TMB_fit(model, data, params, random, starts, method, settings$dots)
    optimizer <<- .msr_gdr_optimizer_checkpoint(result$opt)
    diagnostic <<- .msr_gdr_TMB_diagnostics(result)
    result
  }
  environment(constructor) <- local
  value <- tryCatch(do.call(constructor, args), error = function(e) {
    if (inherits(e, "msr_gdr_tmb_error")) stop(e)
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
.msr_gdr_TMB_fit <- function(model, data, params, random, starts, method, dots) {
  fixed <- names(params)[!names(params) %in% random]
  template <- params[fixed]
  count <- length(unlist(template))
  if (length(starts) != count)
    stop("The number of native TMB constructor starts changed.", call. = FALSE)
  params <- replace(params, fixed, utils::relist(unname(starts), template))
  objective <- TMB::MakeADFun(data = c(model = model, data), parameters = params,
    random = random, silent = TRUE, DLL = "unmarked_TMBExports")
  objective$starts_order <- rep(fixed, lengths(template))
  opt <- do.call(stats::optim, c(list(par = objective$par, fn = objective$fn,
    gr = objective$gr, method = method), dots))
  sdr <- tryCatch({
    value <- objective$fn(opt$par)
    if (!is.finite(value) || !is.finite(opt$value) || abs(value - opt$value) > 1e-7)
      stop("Native TMB likelihood differs at the returned optimizer coordinates.", call. = FALSE)
    TMB::sdreport(objective, par.fixed = opt$par, getJointPrecision = TRUE)
  }, error = function(e) stop(structure(list(message = conditionMessage(e), call = NULL,
      optimizer = .msr_gdr_optimizer_checkpoint(opt), covariance = NULL),
      class = c("msr_gdr_tmb_error", "error", "condition"))))
  sdr$par <- objective$par
  list(opt = opt, TMB = objective, sdr = sdr, AIC = 2 * opt$value + 2 * count)
}

.msr_gdr_refit_start <- function(model, settings) {
  if (settings$engine == "C") return(unmarked::coef(model))
  snapshot <- gdr_start_values(model)
  stats::setNames(snapshot$starts, paste0(snapshot$map$block, "[", snapshot$map$within_block, "]"))
}
