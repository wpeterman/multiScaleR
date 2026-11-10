.msr_gdr_field_binding <- function(object) {
  field <- .msr_gdr_field_guard(object$opt_mod)
  .msr_gdr_hash(list(field$field_id, field$prediction_sha256,
    field$response_sha256, .msr_projection_signature(object)))
}

.msr_gdr_field_geometry <- function(object, coords, crs, batch_size, max_memory_mb) {
  if (!identical(object$backend, "refit"))
    stop("Use a native GDR response-field multiScaleR refit result.", call. = FALSE)
  field <- .msr_gdr_field_guard(object$opt_mod)
  validate_scalar_numeric(batch_size, "batch_size", positive = TRUE, integerish = TRUE)
  validate_scalar_numeric(max_memory_mb, "max_memory_mb", positive = TRUE)
  if (batch_size > .Machine$integer.max || !is.finite(max_memory_mb))
    stop("Use a representable field batch size and finite memory budget.", call. = FALSE)
  declared <- .msr_response_field_crs(crs)
  if (!isTRUE(declared == sf::st_crs(field$crs)) ||
      (!is.matrix(coords) && !is.data.frame(coords)) ||
      !identical(unname(as.matrix(coords)), field$coords))
    stop("Retain the exact fitted GDR field coordinates, row order and projected CRS.", call. = FALSE)
  geometry <- list(schema_version = 2L, provider = "native_GDR_fixed_range",
    crs = declared$wkt, units = declared$units_gdal,
    training_locations = nrow(field$lower), training_rows = nrow(field$coords),
    provider_version = field$provider$TMB_version, batch_size = as.integer(batch_size),
    max_memory_mb = max_memory_mb, range = field$definition$range, range_estimated = FALSE,
    reference_sha256 = .msr_gdr_hash(field$coords), binding_sha256 = .msr_gdr_field_binding(object))
  geometry$definition_sha256 <- .msr_gdr_hash(geometry)
  structure(geometry, class = "multiScaleR_field_geometry")
}

.msr_gdr_field_prediction_budget <- function(geometry, model, rows) {
  field <- model@response_field; p <- length(field$point)
  bytes <- 8 * (32 * as.double(p)^2 + 16 * as.double(nrow(field$coords) + rows) * p)
  if (!is.finite(bytes) || bytes > geometry$max_memory_mb * 2^20)
    stop("GDR response-field prediction exceeds max_memory_mb; reduce batch_size.", call. = FALSE)
  bytes
}

.msr_gdr_field_raster_eval <- function(model, data) {
  geometry <- model$field_geometry; model$field_geometry <- NULL
  internal <- c(".msr_field_x", ".msr_field_y")
  if (!all(internal %in% names(data)) || any(!is.finite(as.matrix(data[internal]))))
    stop("Every GDR field cell needs finite projected field coordinates.", call. = FALSE)
  if (!nrow(data)) return(if (model$se.fit)
    data.frame(fit = numeric(), se.fit = numeric()) else numeric())
  # The coordinate handler bypasses the ordinary raster evaluator. Restore
  # its declared logical types before applying frozen component recipes.
  for (name in intersect(model$logical_names, names(data))) {
    if (any(!is.na(data[[name]]) & !data[[name]] %in% c(0, 1)))
      stop("A logical raster predictor must encode FALSE/TRUE as 0/1.", call. = FALSE)
    data[[name]] <- as.logical(data[[name]])
  }
  values <- lapply(split(seq_len(nrow(data)), ceiling(seq_len(nrow(data))/geometry$batch_size)),
    function(rows) {
      .msr_gdr_field_prediction_budget(geometry, model$object$opt_mod, length(rows))
      part <- data[rows, , drop = FALSE]
      args <- model$native_args
      if (any(c("coords", "crs") %in% names(args)))
        stop("Field geometry supplies raster coords/crs; do not duplicate them in native_args.", call. = FALSE)
      args$coords <- as.matrix(part[internal]); args$crs <- geometry$crs
      value <- .msr_native_prediction_rows(model$object, part, model$type,
        model$se.fit, args, model$unmarked_predictor)
      if (model$se.fit) data.frame(fit = value$Predicted, se.fit = value$SE) else value
    })
  if (model$se.fit) do.call(rbind, values) else unlist(values, use.names = FALSE)
}
