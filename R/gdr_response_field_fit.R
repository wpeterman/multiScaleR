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
    ".add_site_covariate_rasters", ".msr_gdr_TMB_policy", ".msr_gdr_TMB_adfun",
    ".msr_gdr_TMB_phase", ".msr_gdr_TMB_point")
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

#' Fit a Fixed-Range Spatial Abundance Field to a Distance-Removal Model
#'
#' Replace the independent random intercepts on abundance in a native TMB
#' distance-removal (GDR) model with a spatially correlated Gaussian field.
#' Nearby survey blocks then tend to have similar density deviations. The
#' correlation between two blocks falls with the distance `d` between them as
#' `exp(-d / range)`. You supply `range`; it is not estimated. Ecological
#' scales of effect remain separate parameters, estimated later with
#' [make_gdr_refitter()] and [multiScale_optim()].
#'
#' @param model An `unmarkedFitGDR` fitted by `unmarked::gdistremoval()` with
#'   `engine = "TMB"`. It serves as the template. Its data, formulas, cutoff
#'   `K` and estimates are reused, and its independent abundance intercepts
#'   become the field. See "Check that your model qualifies" for the supported
#'   model forms.
#' @param coords Numeric matrix or data frame with two projected coordinate
#'   columns (x, y) and one row per survey, in the row order of the model's
#'   frame. Row names must equal the frame's site row names, which are
#'   `rownames(unmarked::siteCovs(model@data))` when the frame has site
#'   covariates. All surveys in one block (one level of the abundance grouping
#'   factor) must carry exactly the same coordinates, usually the block center,
#'   and different blocks must have different coordinates.
#' @param crs The projected coordinate reference system of `coords`, as an EPSG
#'   code (for example `32617`) or any other input accepted by `sf::st_crs()`.
#'   It must have two horizontal axes and declared linear units. Coordinates are
#'   never transformed, so `range` and `coords` must use these units.
#' @param range Positive number. The fixed correlation range of the
#'   exponential field, in the linear units of `crs`. Two blocks `range` apart
#'   have correlation `exp(-1)`, about 0.37, and blocks three ranges apart have
#'   correlation about 0.05. This is a distance over which density deviations
#'   are shared. It is not a kernel scale of effect, and it is not counted as a
#'   parameter.
#' @param cold_starts Named list of one or more extra starting vectors (for
#'   example `list(cold = cold)`). Each vector holds one value for every fixed
#'   coefficient and log standard deviation, on the link scales listed by
#'   [gdr_start_values()]. Supply it either unnamed, in the order of
#'   `gdr_start_values(model)$starts`, or named by parameter block and
#'   position, such as `"beta_lambda[1]"` (log density intercept),
#'   `"beta_dist[1]"` (log distance scale), `"beta_rem[1]"` (logit removal
#'   probability) and `"lsigma_lambda[1]"` (log SD of the abundance
#'   intercepts). Each vector must differ from the template's estimates and
#'   from the other vectors. The names `"original"` and `"final"` are reserved.
#' @param audit_dir Character string naming a directory for audit records. It
#'   is created if needed. Every fit writes small `.rds` records there: the
#'   request, the validated inputs, each attempt (including failures), the
#'   optimizer result saved before covariance work, and the final outcome.
#'   Records hold numbers, settings and identities, not fitted models. Keep the
#'   directory to inspect a failed fit.
#' @param tmb_settings Named list of optimizer settings with elements `method`,
#'   `se`, `threads` and `control`. This function requires `method = "BFGS"`,
#'   `se = TRUE`, `threads = 1` and a positive integer `control$maxit`, for
#'   example `list(method = "BFGS", se = TRUE, threads = 1L, control =
#'   list(maxit = 2500L, reltol = 1e-10))`. The optional elements
#'   `random.start`, `inner.control` and `record_memory` are described in
#'   [make_gdr_refitter()]. Pass the identical list to [make_gdr_refitter()]
#'   when you later estimate scales of effect with the field.
#' @param max_memory_mb Positive number, default 512. A conservative planning
#'   budget, in MiB, for the dense field covariance and the TMB tape. The fit
#'   stops before optimizing when its estimate exceeds the budget. The same
#'   budget limits later `predict()` calls on the fitted object. It is an
#'   estimate, not an operating-system memory cap.
#' @return An object of class [unmarkedFitGDRField][unmarkedFitGDRField-class].
#'   It extends the `unmarkedFitGDR` class of `unmarked`, so `coef()`, `logLik()`,
#'   `unmarked::sigma()` and `summary()` work as usual, while `predict()`,
#'   `vcov()` and `confint()` use the retained field covariance. The table
#'   under "Read the fitted model" says where each parameter lives and on what
#'   scale. The log-likelihood is the native Laplace approximation of the
#'   marginal likelihood, with the field integrated out. The AIC counts the
#'   same parameters as the template model, because the range is fixed.
#'   The attribute `msr_gdr_refit_audit` summarizes the fit: `status`, the
#'   `attempts` from every start with their diagnostics, the index of the
#'   `selected` start, `maximum_start_value_difference` (the spread of the
#'   negative log-likelihoods among starts) and the audit `path`.
#' @details
#' A random intercept for each survey block lets density differ among blocks
#' beyond what the covariates explain. In the template model those block
#' deviations are independent. This function keeps the same distance and
#' removal observation model, the same cutoff `K` and the same SD of the
#' deviations, but makes the deviations of nearby blocks correlated. See the
#' "Model spatial correlation among survey blocks" section of
#' `vignette("v12-distance-removal", package = "multiScaleR")` for a worked
#' example.
#'
#' # Check that your model qualifies
#'
#' The template must have all of the following:
#'
#' * Poisson abundance (`mixture = "P"`) and a half-normal distance function
#'   (`keyfun = "halfnorm"`).
#' * Density output in birds per square kilometer (`output = "density"`,
#'   `unitsOut = "kmsq"`).
#' * One visit per point (one primary period) and complete counts for every
#'   survey, with no sites dropped by `unmarked`.
#' * Distance bands in meters with boundaries `c(0, 25, 75)` or
#'   `c(0, 25, 50, 75)`, so the survey circle has a 75 m radius.
#' * Exactly one random intercept in the abundance formula, such as
#'   `~ land1 + (1 | block)`. Every level of the grouping factor must have at
#'   least one survey.
#' * Optionally, one independent random intercept in the distance formula,
#'   such as an observer effect. It stays independent. The removal formula
#'   takes no random effects.
#' * No `offset()` terms.
#' * `unmarked` 1.5.2 and `TMB` 1.9.25.
#'
#' The field correlation matrix must also be numerically positive definite.
#' Blocks at nearly identical locations relative to `range` can fail this
#' check. No numerical jitter is added.
#'
#' # How the fit is checked
#'
#' The field model is fitted from the template's own estimates (the start
#' named `"original"`) and from every vector in `cold_starts`. Each fit starts
#' its block deviations at zero. Every start must return a converged,
#' warning-free fit with usable covariance and a largest absolute gradient
#' below 0.001. All starts must also agree, with negative log-likelihoods
#' within 1e-5 and coefficients and log SDs within 0.001 of one another. If
#' any start fails or the starts disagree, the function stops with an error
#' and leaves the records in `audit_dir`. When every check passes, it returns
#' the fit with the highest likelihood. Agreement among a few starts does not
#' prove that no better optimum exists. If the starts disagree, compare the
#' attempts recorded in `audit_dir` before trying other starting values.
#'
#' # Read the fitted model
#'
#' | Quantity | Where to find it | Scale | How to read it |
#' |:---|:---|:---|:---|
#' | Density intercept | `coef(fit)["lambda(Int)"]` | log | `exp()` gives density per square km when covariates and the field are zero |
#' | Abundance covariate effect | `coef(fit)["lambda(x)"]` | log | `exp()` is the factor by which density changes per unit of `x` |
#' | Distance scale | `coef(fit)["dist(Int)"]` | log | `exp()` gives the half-normal scale in meters |
#' | Removal probability | `coef(fit)["rem(Int)"]` | logit | `plogis()` gives the probability, per unit of interval time, that a bird present is first detected |
#' | Field SD | `unmarked::sigma(fit)` | log density | Typical size of a block deviation. A block one SD above zero has `exp(SD)` times the density of a block with no deviation |
#' | Block deviations | `coef(fit, type = "lambda", fixedOnly = FALSE)` | log density | The entries labeled `lambda(b_lambda)`, one fitted field value per block, follow the fixed coefficients |
#' | Correlation range | `fit@response_field$definition$range` | CRS units | Supplied by you, not estimated |
#'
#' Internally, TMB works with uncorrelated (whitened) versions of the block
#' deviations. Every value reported by `coef()`, `vcov()` and `predict()` is
#' converted back to the physical field, the deviation at each block's
#' location.
#'
#' # How uncertainty is reported
#'
#' [vcov()][vcov,unmarkedFitGDRField-method] returns the covariance retained
#' from the fit. It includes the covariance between components (abundance,
#' distance, removal) and between the fixed coefficients and the field values.
#' It is computed once during fitting, so `vcov()` does not rerun TMB.
#'
#' [predict()][predict,unmarkedFitGDRField-method] uses that covariance for
#' its standard errors. At a block used in fitting, the SE comes from the
#' fitted coefficients and field values alone. At a new location, the SE adds
#' the field variance that remains after conditioning on the fitted blocks.
#' That extra variance describes how uncertain the field is between blocks.
#' It is not observation noise. Density predictions are `exp()` of the linear
#' predictor at the fitted field values. They are not averages over possible
#' field values, which would be larger. With `re.form = NA`, predictions set
#' the field to zero, which describes a block with no deviation rather than an
#' average over blocks.
#'
#' All of these intervals treat the ecological scales of effect and the
#' correlation range as known. They add no uncertainty for either.
#'
#' # Estimate scales of effect with the field in place
#'
#' Build a refit function from the field fit, using the same starts and the
#' identical `tmb_settings`, and pass it to [multiScale_optim()]:
#'
#' ```
#' refit <- make_gdr_refitter(field_fit, list(cold = cold),
#'   audit_dir = tempfile("gdr-field-refits-"), tmb_settings = settings)
#' opt <- multiScale_optim(field_fit, kernel_inputs, refit_fn = refit,
#'   n_cores = 1)
#' ```
#'
#' Each candidate scale refits the field from every start, with the range held
#' fixed. Expect minutes even for 100 surveys. The search must run on one core
#' with the default `backend = "refit"` and `PSOCK = FALSE`. The ordinary
#' `update()` method stops for this class. Use the refit function instead.
#'
#' # Map density with the field
#'
#' A density map needs the field value at every raster cell. Bind the fitting
#' coordinates to the scale-of-effect result with [response_field_geometry()],
#' then pass the descriptor to [predict.multiScaleR()]:
#'
#' ```
#' geometry <- response_field_geometry(opt, coords = block_xy, crs = 32617)
#' layers <- kernel_scale.raster(land, multiScaleR = opt,
#'   scale_center = TRUE, match_buffer = TRUE)
#' dens <- predict(opt, layers, type = "lambda", se.fit = TRUE,
#'   uncertainty = "conditional",
#'   raster_args = list(field_geometry = geometry))
#' ```
#'
#' The raster CRS must match `crs`, and prediction runs on one core. Each cell
#' gets the field value predicted at its center from the fitted block values.
#' The grouping factor of the field (`block` above) is not needed as a
#' raster. A factor that also enters a fixed part of the model still needs an
#' explicit value through `site_values` in [kernel_scale.raster()] or a raster
#' coded with the fitted levels. To map the zero-field reference instead, omit
#' `field_geometry` and set `re.form = NA`.
#'
#' # Not supported
#'
#' Estimating the range, other abundance mixtures or distance functions,
#' repeat visits with availability, a removal field, `simulate()`,
#' `profile()`, `ranef()` and spatial cross-validation are not available for
#' this model.
#' @seealso [unmarkedFitGDRField-class], [make_gdr_refitter()],
#'   [gdr_start_values()], [response_field_geometry()],
#'   `vignette("v12-distance-removal", package = "multiScaleR")`
#' @examples
#' \donttest{
#' # The field needs unmarked 1.5.2 and TMB 1.9.25.
#' if (requireNamespace("TMB", quietly = TRUE) &&
#'     utils::packageVersion("unmarked") == "1.5.2" &&
#'     utils::packageVersion("TMB") == "1.9.25") {
#'   # Eight survey blocks 1 km apart along a line, six point counts in each
#'   block <- factor(rep(paste0("b", 1:8), each = 6))
#'   hab <- rep(seq(-1, 1, length.out = 6), 8)
#'   shift <- c(-0.4, -0.3, -0.1, 0, 0.1, 0.2, 0.35, 0.4)[block]
#'   sim <- sim_dat_gdr(density = exp(log(400) + 0.5 * hab + shift),
#'     distance_scale = 40, removal = 0.2, dist_breaks = c(0, 25, 50, 75),
#'     interval_lengths = c(2, 3, 5), interval_units = "minutes",
#'     site_covs = data.frame(hab = hab, block = block), seed = 1)
#'   umf <- sim$observed
#'
#'   # Template: an independent random intercept for each block (TMB engine)
#'   template <- unmarked::gdistremoval(lambdaformula = ~ hab + (1 | block),
#'     removalformula = ~ 1, distanceformula = ~ 1, data = umf,
#'     keyfun = "halfnorm", mixture = "P", output = "density",
#'     unitsOut = "kmsq", K = 60, engine = "TMB",
#'     starts = c(log(400), 0, log(40), qlogis(0.2), log(0.3)))
#'
#'   # Every survey carries the projected coordinates (m) of its block center
#'   block_xy <- cbind(x = 500000 + 1000 * (as.integer(block) - 1),
#'                     y = 4400000)
#'   rownames(block_xy) <- rownames(unmarked::siteCovs(umf))
#'
#'   # A second set of starting values, named by parameter block and position
#'   cold <- c("beta_lambda[1]" = log(300), "beta_lambda[2]" = 0,
#'             "beta_dist[1]" = log(30), "beta_rem[1]" = qlogis(0.4),
#'             "lsigma_lambda[1]" = log(0.5))
#'   settings <- list(method = "BFGS", se = TRUE, threads = 1L,
#'                    control = list(maxit = 2500L, reltol = 1e-10))
#'
#'   field_fit <- fit_gdr_response_field(template, coords = block_xy,
#'     crs = 32617, range = 2000, cold_starts = list(cold = cold),
#'     audit_dir = tempfile("gdr-field-"), tmb_settings = settings)
#'
#'   # Log density, log distance scale (m) and logit removal probability
#'   print(unmarked::coef(field_fit))
#'   # SD of the block deviations, on the log density scale
#'   print(unmarked::sigma(field_fit))
#'   # The range adds no parameter, so the AIC penalty matches the template
#'   print(c(template = template@AIC, field = field_fit@AIC))
#' }
#' }
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
    attempt$native_preparation_started <- .msr_gdr_TMB_phase("native_preparation_started",
      proc.time()[["elapsed"]], settings$policy)
    .msr_gdr_atomic_audit(attempt, sub("[.]rds$", "-request.rds", path))
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
      checkpoint <- function(point) {
        point$native_preparation <- list(started = attempt$native_preparation_started,
          elapsed_seconds = point$phases$construction_started$process_elapsed_seconds -
            attempt$native_preparation_started$process_elapsed_seconds,
          scope = "native_constructor_and_field_loading_before_TMB_tape_helper")
        point$binding <- list(field_id = field$field_id, provider = field$provider,
          source = field$source, input_design_sha256 = record$input_design_sha256,
          start_name = name, input_field_configuration_sha256 = .msr_gdr_field_configuration_hash(field))
        .msr_gdr_atomic_audit(point, sub("[.]rds$", "-point.rds", path))
      }
      result <<- .msr_gdr_TMB_fit(model, data, params, random, starts, method, settings$dots,
        policy = settings$policy, checkpoint = checkpoint)
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
        if (inherits(e, "msr_gdr_audit_error")) stop(e)
        attempt$status <<- "failed"; attempt$error <<- conditionMessage(e)
        if (inherits(e, "msr_gdr_tmb_error")) {
          attempt$optimizer_point_file <<- e$point_file; attempt$workflow <<- e$workflow
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
    attempt$optimizer_point_file <- result$workflow$point_file
    attempt$native_preparation_elapsed_seconds <-
      result$workflow$phases$construction_started$process_elapsed_seconds -
      attempt$native_preparation_started$process_elapsed_seconds
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
