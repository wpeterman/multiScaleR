# Find raw offset inputs without evaluating formula environments. Double-sided
# unmarked formulas and glmmTMB component formulas can contain several offsets.
.msr_raster_offset_terms <- function(expr) {
  if (is.list(expr)) return(unlist(lapply(expr, .msr_raster_offset_terms), recursive = FALSE))
  if (!is.call(expr)) return(list())
  fun <- as.character(expr[[1L]])[1L]
  if (fun == "~") {
    left <- if (length(expr) == 3L && is.call(expr[[2L]]) &&
      identical(expr[[2L]][[1L]], as.name("~"))) .msr_raster_offset_terms(expr[[2L]]) else list()
    return(c(left, .msr_raster_offset_terms(expr[[length(expr)]])))
  }
  if (fun == "offset") return(list(expr[[2L]]))
  unlist(lapply(as.list(expr)[-1L], .msr_raster_offset_terms), recursive = FALSE)
}

.msr_site_values <- function(values) {
  if (is.null(values) || !length(values)) return(list())
  if (!is.list(values)) values <- as.list(values)
  if (is.null(names(values)) || anyNA(names(values)) || any(!nzchar(names(values))) ||
      anyDuplicated(names(values)) || any(vapply(values, function(x)
        !is.atomic(x) || length(x) != 1L || anyNA(x) ||
          (is.numeric(x) && !is.finite(x)), logical(1))))
    stop("site_values must name distinct, nonmissing scalar scenario values; use a list for mixed types.", call. = FALSE)
  values
}

# Append explicit nonspatial scenarios while retaining historical numeric-zero
# defaults. Categorical rasters retain all fitted labels, not arbitrary codes.
.add_site_covariate_rasters <- function(smooth_stack, multiScaleR, raster_covs,
                                        offset_values = NULL, site_values = NULL) {
  if (is.null(multiScaleR) || !inherits(multiScaleR, "multiScaleR")) return(smooth_stack)
  mod <- .analysis_model(multiScaleR$opt_mod)
  if (is.null(mod)) return(smooth_stack)
  site_values <- .msr_site_values(site_values)
  if (any(names(site_values) %in% names(offset_values)))
    stop("Supply each offset scenario in site_values or offset_values, not both.", call. = FALSE)
  unmarked <- methods::is(mod, "unmarkedFit")
  mf <- if (unmarked) mod@data@siteCovs else
    tryCatch(stats::model.frame(mod), error = function(e) NULL)
  raw <- if (unmarked) mf else tryCatch(.model_data(mod, effects = "all"),
                                        error = function(e) NULL)
  if (is.null(mf) && is.null(raw) && !length(site_values)) {
    warning("Could not inspect fitted model frame; site-level covariates were not added.", call. = FALSE)
    return(smooth_stack)
  }
  formulas <- if (unmarked) .unmarked_model_formulas(mod) else list(stats::formula(mod))
  if (inherits(mod, "glmmTMB")) formulas <- c(formulas,
    mod$modelInfo$allForm[c("ziformula", "dispformula")])
  offsets <- unlist(lapply(formulas, .msr_raster_offset_terms), recursive = FALSE)
  offset_vars <- character()
  offset_defaults <- numeric()
  available <- c(names(offset_values), names(site_values), raster_covs)
  for (expr in offsets) {
    vars <- all.vars(expr)
    if (!length(vars)) next
    simple_log <- is.call(expr) && identical(expr[[1L]], as.name("log")) &&
      length(expr) == 2L && is.symbol(expr[[2L]])
    simple <- is.symbol(expr) || simple_log
    if (!simple && !all(vars %in% available))
      stop("Offset expression `", deparse(expr),
        "` needs explicit raw values for nonspatial inputs in offset_values or site_values.", call. = FALSE)
    for (v in vars) {
      value <- if (simple_log) 1 else 0
      if (v %in% names(offset_defaults) && offset_defaults[[v]] != value &&
          !v %in% available)
        stop("Conflicting offset reference values for '", v, "'; supply it explicitly.", call. = FALSE)
      offset_defaults[v] <- value
    }
    offset_vars <- union(offset_vars, vars)
  }
  if (length(setdiff(names(offset_values), offset_vars)))
    stop("offset_values contains names not used as formula offsets.", call. = FALSE)
  model_vars <- if (unmarked) intersect(.unmarked_model_predictors(mod),
    union(names(mf), names(site_values))) else union(.model_predictors(mod),
      unique(unlist(lapply(formulas, function(f) if (is.null(f)) character() else
        all.vars(f[[length(f)]])), use.names = FALSE)))
  site_vars <- setdiff(union(model_vars, offset_vars), raster_covs)
  # GDR component maps use explicit scenarios. In particular, do not create
  # zero rasters for unprojected spatial covariates from other components.
  if (methods::is(mod, "unmarkedFitGDR"))
    site_vars <- intersect(site_vars, union(names(site_values), offset_vars))
  if (length(setdiff(names(site_values), site_vars)))
    stop("site_values must name nonraster model predictors; spatial layers cannot be overwritten.", call. = FALSE)
  if (any(names(offset_values) %in% raster_covs))
    stop("offset_values cannot overwrite a computed spatial layer.", call. = FALSE)
  template <- smooth_stack[[1L]]
  for (v in site_vars) {
    x <- if (!is.null(mf) && v %in% names(mf)) mf[[v]] else raw[[v]]
    if (v %in% names(site_values)) {
      value <- site_values[[v]]
      if (is.factor(x) || is.character(x) || is.character(value) || is.factor(value)) {
        fitted_levels <- if (is.list(mod)) mod$xlevels[[v]] else NULL
        if (is.null(fitted_levels)) fitted_levels <- if (is.factor(x)) levels(x) else
          if (is.character(x)) levels(as.factor(x)) else NULL
        if (!length(fitted_levels) || !as.character(value) %in% fitted_levels)
          stop("site_values must select a fitted level for '", v, "'.", call. = FALSE)
        rr <- terra::as.factor(template * 0 + match(as.character(value), fitted_levels))
        levels(rr) <- data.frame(ID = seq_along(fitted_levels), label = fitted_levels)
      } else {
        if (!(is.numeric(value) || is.logical(value)))
          stop("A numeric site predictor needs a finite numeric scenario value.", call. = FALSE)
        rr <- template * 0 + as.numeric(value)
      }
    } else if (v %in% offset_vars) {
      value <- if (v %in% names(offset_values)) offset_values[[v]] else offset_defaults[[v]]
      rr <- template * 0 + value
    } else if (is.factor(x) || is.character(x)) {
      warning(paste0("Site-level covariate '", v, "' is categorical. ",
        "A dummy raster cannot be added safely without explicitly defining a projection value."),
        call. = FALSE)
      next
    } else if (is.numeric(x) || is.integer(x) || is.logical(x)) {
      rr <- template * 0
    } else {
      warning(paste0("Site-level covariate '", v,
        "' is not numeric/integer/logical. Skipping dummy raster creation."), call. = FALSE)
      next
    }
    names(rr) <- v
    smooth_stack <- c(smooth_stack, rr)
  }
  smooth_stack
}
