# Native GDR must describe the fitted likelihood before it can be optimized,
# compared or projected. Qualification failures are not repaired by dropping
# parameters, switching engines or returning apparently valid uncertainty.
.msr_gdr_native_guard <- function(object) {
  if (!methods::is(object, "unmarkedFitGDR")) return(invisible(NULL))
  if (identical(object@mixture, "ZIP") && !is.null(object@TMB))
    stop("GDR TMB ZIP models are not qualified for multiScaleR fitting, selection or prediction: the audited native likelihood ignores the zero-inflation parameter. Use a separately validated model; do not change its parameter penalty or silently switch engines.",
      call. = FALSE)

  formulas <- .msr_gdr_active_formulas(object)
  # terms() recognizes actual R offsets. Walking arbitrary calls would also
  # reject ordinary I(offset(x)) predictors and quoted column names. A zero-row
  # frame supplies dot-expansion names without evaluating training transforms.
  columns <- unique(unlist(lapply(c("siteCovs", "yearlySiteCovs", "obsCovs"),
    function(name) names(methods::slot(object@data, name))), use.names = FALSE))
  data <- as.data.frame(stats::setNames(rep(list(numeric()), length(columns)),
    columns), check.names = FALSE)
  for (component in names(formulas)) {
    formula <- formulas[[component]]
    if (!inherits(formula, "formula"))
      stop("Could not inspect the active GDR ", component, " formula.", call. = FALSE)
    tt <- tryCatch(stats::terms(formula, data = data), error = function(e)
      stop("Could not inspect the active GDR ", component, " formula: ",
        conditionMessage(e), call. = FALSE))
    if (length(attr(tt, "offset")))
      stop("Active GDR ", component, " offsets are not qualified for multiScaleR fitting, selection or prediction: the audited native engine omits their contribution. Supply a separately validated model; do not remove an offset from the intended model.",
        call. = FALSE)
  }
  invisible(NULL)
}

# Availability and distance designs can be inactive in native GDR models.
# They must not create an optimized ecological scale or a fitted predictor.
.msr_gdr_active_formulas <- function(object) {
  active <- c("lambda", if (object@data@numPrimary > 1L) "phi",
    if (!identical(object@keyfun, "uniform")) "dist", "rem")
  stats::setNames(lapply(active, function(component) object@formlist[[component]]), active)
}

# Expand dot terms against the same component's available native columns,
# without evaluating transformations or learning from prediction rows.
.msr_gdr_model_predictors <- function(object) {
  formulas <- .msr_gdr_active_formulas(object)
  site <- names(object@data@siteCovs)
  yearly <- unique(c(names(object@data@yearlySiteCovs), site))
  observation <- unique(c(names(object@data@obsCovs), yearly))
  columns <- list(lambda = site, phi = yearly, dist = yearly, rem = observation)
  as.character(unique(unlist(lapply(names(formulas), function(component) {
    component_names <- columns[[component]]
    data <- as.data.frame(stats::setNames(rep(list(numeric()), length(component_names)),
      component_names), check.names = FALSE)
    all.vars(attr(stats::delete.response(stats::terms(formulas[[component]], data = data)), "variables"))
  }), use.names = FALSE)))
}

# Native GDR combines observation/yearly tables with repeated static columns.
# A duplicated optimized name cannot identify one spatial design unambiguously.
.msr_gdr_spatial_guard <- function(object, covariates) {
  if (!methods::is(object, "unmarkedFitGDR")) return(invisible(NULL))
  shadow <- intersect(covariates, union(names(object@data@yearlySiteCovs),
    names(object@data@obsCovs)))
  if (length(shadow))
    stop("Optimized GDR site covariates are shadowed in yearlySiteCovs or obsCovs: ",
      paste(shadow, collapse = ", "), ". Use distinct names for distinct native inputs.",
      call. = FALSE)
  invisible(NULL)
}
