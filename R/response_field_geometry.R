# A planar two-axis map is required, including when reading an older descriptor.
# Excluding longitude/latitude alone would also admit geocentric/vertical CRS.
.msr_response_field_crs <- function(crs) {
  declared <- sf::st_crs(crs)
  if (!is.na(declared)) {
    wkt <- gsub("[[:space:]]", "", declared$wkt)
    projected <- any(startsWith(wkt, c("PROJCRS[", "PROJCS[",
      "BOUNDCRS[SOURCECRS[PROJCRS[", "BOUNDCRS[SOURCECRS[PROJCS[")))
    if (projected && grepl("CS[Cartesian,2]", wkt, fixed = TRUE) &&
        !is.null(declared$units_gdal) && !is.na(declared$units_gdal) &&
        nzchar(declared$units_gdal)) return(declared)
  }
  stop("Use a known projected CRS with two horizontal dimensions and declared linear units.", call. = FALSE)
}

# Bind a field descriptor to the retained fit without retaining a model or tape.
.msr_response_field_binding <- function(object) {
  model <- .analysis_model(object$opt_mod)
  .msr_gdr_hash(list(frame = list(rows = rownames(model$frame),
    columns = lapply(model$frame, identity)),
    fit = model$fit[c("par", "parfull", "objective", "convergence")],
    fields = .msr_glmmTMB_response_fields(model),
    projection = .msr_projection_signature(object)))
}

#' Bind Projected Coordinates to a Native Response Field
#'
#' Verify the physical coordinates used by a fitted response field before making
#' coordinate-aware raster predictions. The descriptor stores metadata and
#' hashes, without a model, coordinate array or prediction array.
#'
#' @param object A \code{multiScaleR} refit result containing a native glmmTMB
#'   Gaussian, log-link Poisson/NB2/Gamma or logit-binomial model with one shared
#'   exponential, Gaussian or Matern response field and no other random effects.
#' @param coords Numeric matrix or data frame with two finite projected coordinate
#'   columns, one row per retained native fitting row in its original order.
#'   Repeated locations must repeat their coordinate pairs. Use the original
#'   coordinates supplied to \code{glmmTMB::numFactor()}, in the same units.
#' @param crs Known two-dimensional projected horizontal coordinate reference
#'   system, accepted by
#'   \code{sf::st_crs()}. Its linear units must describe \code{coords} exactly.
#'   Prediction rasters must have an equivalent CRS. Coordinates are never
#'   transformed or rescaled automatically. Geographic, geocentric, vertical
#'   and compound definitions are rejected. Choose the two-dimensional map
#'   CRS explicitly before fitting and prediction.
#' @param batch_size Positive integer maximum number of rows passed to each
#'   native prediction call. The default is 128. Smaller batches reduce dense
#'   field expansion but repeat native prediction preparation.
#' @param max_memory_mb Positive finite working-memory budget in MiB, default
#'   512. A conservative dense-array estimate must fit this budget before native
#'   prediction. This planning gate is not an operating-system memory cap and
#'   cannot account exactly for opaque native automatic-differentiation tapes.
#' @return A \code{multiScaleR_field_geometry} descriptor with schema version,
#'   projected \code{crs} WKT, \code{units}, native \code{position} and
#'   \code{group} column names, the fitted \code{group_value}/\code{group_levels},
#'   \code{training_locations}/\code{training_rows} counts, \code{provider_version},
#'   prediction controls \code{batch_size}/\code{max_memory_mb}, and SHA256
#'   \code{reference_sha256}, \code{binding_sha256} and \code{definition_sha256}.
#'   The hashes bind coordinate verification, fitted context and descriptor
#'   contents respectively. Recreate the descriptor after changing its inputs.
#' @details Supply this object as
#'   \code{raster_args = list(field_geometry = geometry)} in
#'   \code{predict.multiScaleR()}. Cells supply their own projected coordinates.
#'   Their native position factor replaces any constant position scenario.
#'   Every cell uses the verified shared field group. An inconsistent group
#'   scenario fails rather than selecting another field. Prediction uses one
#'   worker, with serial native sub-batches. No full-grid covariance is retained.
#'   For previously unobserved positions, explicitly set
#'   \code{native_args = list(allow.new.levels = TRUE)}. A missing or false
#'   setting fails before that batch enters native prediction. Native warnings
#'   from otherwise supported calls pass through unchanged.
#'
#'   Prepare ecological layers with \code{kernel_scale.raster()},
#'   \code{scale_center = TRUE} and \code{match_buffer = TRUE}. Supply ordinary
#'   site and offset scenarios as usual. Position/group scenarios can select
#'   fitted levels during raster preparation. During prediction, the geometry
#'   descriptor supplies cell-specific positions. Do not supply a constant
#'   physical coordinate in place of this binding.
#'
#'   Native mean-prediction SEs remain conditional on selected ecological kernel
#'   scales. They do not add observation noise or kernel-scale uncertainty.
#'   Binomial output keeps the native per-trial probability scale, including
#'   when fitting grouped success/failure counts. Poisson/NB2 output is the
#'   native expected count; Gaussian/Gamma output is the native response mean.
#'   \code{re.form = NA} gives a zero-field reference, not integration over the
#'   random field. Output retains fit and requested SE bands even if all cells
#'   are NA. Multiple rows in \code{raster_args$const} produce one fit/SE band
#'   pair per scenario row, in that order, with names made unique. Scenarios
#'   are evaluated serially without a full-grid coordinate or scenario array.
#'   Explicit \code{wopt$names} are retained. Cross-row covariance requires point predictions. Coordinate
#'   verification does not establish recovery, calibrated intervals, spatial CV
#'   scoring or eligibility for the guarded \code{joint_mixed} backend.
#' @seealso \code{\link{predict.multiScaleR}}, \code{\link{kernel_scale.raster}},
#'   \code{\link{prediction_contract}}
#' @examples
#' # After fitting a native response field, use its retained fitting rows:
#' # geometry <- response_field_geometry(fit, coords = xy, crs = 32617)
#' # map <- predict(fit, layers,
#' #   native_args = list(allow.new.levels = TRUE),
#' #   raster_args = list(field_geometry = geometry))
#' @export
response_field_geometry <- function(object, coords, crs, batch_size = 128L,
                                    max_memory_mb = 512) {
  if (!inherits(object, "multiScaleR") || !identical(object$backend, "refit"))
    stop("Use a native glmmTMB multiScaleR refit result.", call. = FALSE)
  model <- .analysis_model(object$opt_mod)
  fields <- .msr_glmmTMB_response_fields(model)
  if (is.null(fields) || !identical(fields$status, "native_layout_identified") ||
      length(fields$fields) != 1L || length(model$modelInfo$reStruc$condReStruc) != 1L ||
      length(model$modelInfo$reStruc$ziReStruc) ||
      !fields$fields[[1L]]$covariance %in% c("exponential", "gaussian", "matern") ||
      fields$fields[[1L]]$independent_fields != 1L)
    stop("Coordinate-aware maps require one shared native exponential, Gaussian or Matern response field.", call. = FALSE)
  field <- fields$fields[[1L]]
  if (length(field$coordinate_variables) != 1L || length(field$grouping_variables) != 1L)
    stop("The field needs simple position and group column names.", call. = FALSE)
  family <- model$modelInfo$family
  qualified <- family$family == "gaussian" ||
    (family$family %in% c("poisson", "nbinom2", "Gamma") && family$link == "log") ||
    (family$family == "binomial" && family$link == "logit")
  if (!qualified)
    stop("Coordinate-aware field maps require Gaussian, log-Poisson, log-NB2, log-Gamma or logit-binomial models.", call. = FALSE)
  validate_scalar_numeric(batch_size, "batch_size", integerish = TRUE, positive = TRUE)
  validate_scalar_numeric(max_memory_mb, "max_memory_mb", positive = TRUE)
  if (batch_size > .Machine$integer.max || !is.finite(max_memory_mb))
    stop("Use a representable batch size and finite memory budget.", call. = FALSE)
  reference_crs <- .msr_response_field_crs(crs)
  if ((!is.matrix(coords) && !is.data.frame(coords)) || ncol(coords) != 2L ||
      !is.numeric(as.matrix(coords)) || any(!is.finite(as.matrix(coords))) ||
      nrow(coords) != nrow(model$frame))
    stop("coords must contain finite coordinate pairs for every retained native fitting row.", call. = FALSE)
  xy <- unname(as.matrix(coords))
  position <- field$coordinate_variables
  group <- field$grouping_variables
  if (!is.factor(model$frame[[position]]) ||
      !identical(as.character(glmmTMB::numFactor(xy)), as.character(model$frame[[position]])) ||
      nrow(unique(xy)) != field$locations_per_field)
    stop("Coordinates must match native position factors in retained row order and units.", call. = FALSE)
  groups <- model$frame[[group]]
  if (!is.factor(groups) || anyNA(groups) || length(unique(as.character(groups))) != 1L)
    stop("The retained response field must use one nonmissing shared group.", call. = FALSE)
  geometry <- list(schema_version = 1L, crs = reference_crs$wkt,
    units = reference_crs$units_gdal, position = position, group = group,
    group_value = as.character(groups[1L]), group_levels = levels(groups),
    training_locations = field$locations_per_field, training_rows = nrow(model$frame),
    provider_version = fields$provider_version, batch_size = as.integer(batch_size),
    max_memory_mb = max_memory_mb, reference_sha256 = .msr_gdr_hash(xy),
    binding_sha256 = .msr_response_field_binding(object))
  geometry$definition_sha256 <- .msr_gdr_hash(geometry)
  structure(geometry, class = "multiScaleR_field_geometry")
}

.msr_response_field_budget <- function(geometry, model, rows) {
  # Allow substantial dense covariance/factorization/AD workspace. The opaque
  # native tape is not measurable by object.size(), so do not call this a cap.
  m <- as.double(geometry$training_locations) + rows
  p <- length(model$fit$par)
  bytes <- 8 * (128 * m^2 + 16 * (geometry$training_rows + rows) * (m + p) + 8 * p^2) +
    2 * as.numeric(utils::object.size(model))
  if (!is.finite(bytes) || bytes > geometry$max_memory_mb * 2^20)
    stop("Native response-field expansion exceeds max_memory_mb; reduce batch_size or use a separately validated field representation.", call. = FALSE)
  bytes
}

.msr_response_field_preflight <- function(geometry, object, raster, cores, const) {
  if (!inherits(geometry, "multiScaleR_field_geometry") ||
      !identical(geometry$definition_sha256,
        .msr_gdr_hash(unclass(geometry)[setdiff(names(geometry), "definition_sha256")])))
    stop("Recreate the unchanged response_field_geometry() descriptor.", call. = FALSE)
  if (!identical(geometry$binding_sha256, .msr_response_field_binding(object)))
    stop("Field geometry must match this retained fitted context.", call. = FALSE)
  if (!isTRUE(sf::st_crs(terra::crs(raster)) == .msr_response_field_crs(geometry$crs)))
    stop("Prediction raster CRS and linear units must match field geometry; transform explicitly.", call. = FALSE)
  if (!is.numeric(cores) || length(cores) != 1L || !identical(as.numeric(cores), 1))
    stop("Coordinate-aware field maps require cores = 1 and serial native batches.", call. = FALSE)
  internal <- c(".msr_field_x", ".msr_field_y")
  if (any(internal %in% c(names(raster), names(const))))
    stop("Reserved field coordinate names conflict with raster predictors or constants.", call. = FALSE)
  .msr_response_field_budget(geometry, .analysis_model(object$opt_mod),
    min(geometry$batch_size, terra::ncell(raster)))
}

.msr_response_field_raster_eval <- function(model, data) {
  constants <- model$field_constants
  if (!is.null(constants)) {
    model$field_constants <- NULL
    bands <- if (model$se.fit) 2L else 1L
    output <- matrix(NA_real_, nrow(data), bands * nrow(constants))
    for (i in seq_len(nrow(constants))) {
      if (isTRUE(model$field_na_rm) && anyNA(constants[i, , drop = FALSE])) next
      scenario <- data
      for (nm in names(constants)) scenario[[nm]] <- rep(constants[[nm]][i], nrow(data))
      value <- .msr_response_field_raster_eval(model, scenario)
      output[, (i - 1L) * bands + seq_len(bands)] <- as.matrix(value)
    }
    colnames(output) <- make.names(rep(if (model$se.fit) c("fit", "se.fit") else "fit",
      nrow(constants)), unique = TRUE)
    return(output)
  }
  geometry <- model$field_geometry
  model$field_geometry <- NULL
  xy_names <- c(".msr_field_x", ".msr_field_y")
  if (!all(xy_names %in% names(data)) ||
      any(!is.finite(as.matrix(data[xy_names]))))
    stop("Every predicted cell needs finite projected field coordinates.", call. = FALSE)
  if (geometry$group %in% names(data) && (anyNA(data[[geometry$group]]) ||
      any(as.character(data[[geometry$group]]) != geometry$group_value)))
    stop("Raster scenarios must retain the verified shared field group.", call. = FALSE)
  if (!nrow(data)) return(if (model$se.fit)
    data.frame(fit = numeric(), se.fit = numeric()) else numeric())
  indices <- split(seq_len(nrow(data)), ceiling(seq_len(nrow(data)) / geometry$batch_size))
  values <- lapply(indices, function(rows) {
    .msr_response_field_budget(geometry, .analysis_model(model$object$opt_mod), length(rows))
    part <- data[rows, , drop = FALSE]
    part[[geometry$position]] <- glmmTMB::numFactor(as.matrix(part[xy_names]))
    native <- .analysis_model(model$object$opt_mod)
    retained <- levels(native$frame[[geometry$position]])
    added <- setdiff(levels(part[[geometry$position]]), retained)
    if (!isTRUE(model$native_args[["allow.new.levels"]]) && length(added))
      stop("New response-field positions require native_args = list(allow.new.levels = TRUE).", call. = FALSE)
    # Keep original factor coding even in a single-location sub-batch.
    part[[geometry$position]] <- factor(as.character(part[[geometry$position]]),
                                       levels = c(retained, added))
    part[[geometry$group]] <- factor(rep(geometry$group_value, length(rows)),
                                   levels = geometry$group_levels)
    .msr_native_raster_eval(model, part)
  })
  if (is.null(dim(values[[1L]]))) unlist(values, use.names = FALSE) else do.call(rbind, values)
}
