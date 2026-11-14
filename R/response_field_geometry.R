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

# The default refit backend stores no `backend` element on its result. Only
# the opt-in backends name themselves. Treat a missing element as the refit
# backend, as the refit Hessian code does.
.msr_refit_backend_fit <- function(object) {
  is.null(object$backend) || identical(object$backend, "refit")
}

# Bind a field descriptor to the retained fit without retaining a model or tape.
.msr_response_field_binding <- function(object) {
  model <- .analysis_model(object$opt_mod)
  if (methods::is(model, "unmarkedFitGDRField")) return(.msr_gdr_field_binding(object))
  .msr_gdr_hash(list(frame = list(rows = rownames(model$frame),
    columns = lapply(model$frame, identity)),
    fit = model$fit[c("par", "parfull", "objective", "convergence")],
    fields = .msr_glmmTMB_response_fields(model),
    projection = .msr_projection_signature(object)))
}

#' Bind Projected Coordinates to a Fitted Response Field
#'
#' Check the coordinates used to fit a spatial response field and package
#' them, with their CRS, into a small descriptor for raster prediction. A
#' response field is a spatially correlated random effect on the response,
#' separate from the kernel scales of effect. To map such a model, each
#' raster cell needs its own field value, and the descriptor tells
#' [predict.multiScaleR()] how to compute it from cell coordinates.
#'
#' @param object A `multiScaleR` result from [multiScale_optim()] with the
#'   default `backend = "refit"`, whose fitted model is one of:
#'   * a `glmmTMB` model with one shared exponential, Gaussian or Matern
#'     field, such as `exp(pos + 0 | grp)`, and no other random effects. The
#'     family must be Gaussian, Poisson, NB2 or Gamma with a log link, or
#'     binomial with a logit link.
#'   * a distance-removal model with a fixed-range abundance field from
#'     [fit_gdr_response_field()].
#' @param coords Numeric matrix or data frame with two finite projected
#'   coordinate columns (x, y), one row per fitting row in the original row
#'   order. Repeated locations repeat their coordinates. For `glmmTMB`, use the
#'   coordinates passed to `glmmTMB::numFactor()`, in the same units. For a
#'   distance-removal field, use exactly the `coords` passed to
#'   [fit_gdr_response_field()]. Those are block locations, which can be
#'   coarser than the survey points used for the kernel covariates.
#' @param crs The projected CRS of `coords`, as an EPSG code (for example
#'   `32617`) or any input accepted by `sf::st_crs()`. It must have two
#'   horizontal axes and linear units that match `coords`. Geographic,
#'   geocentric, vertical and compound CRSs are rejected. Prediction rasters
#'   must use an equivalent CRS. Nothing is transformed or rescaled, so choose
#'   the projected CRS before fitting.
#' @param batch_size Positive integer, default 128. The largest number of
#'   raster cells passed to each native prediction call. Smaller batches use
#'   less memory but repeat the native setup more often.
#' @param max_memory_mb Positive number, default 512. A conservative budget,
#'   in MiB, for the dense arrays of each prediction batch. A batch whose
#'   estimate exceeds the budget stops before native prediction. Lower
#'   `batch_size` if that happens. This is a planning estimate, not an
#'   operating-system memory cap. For a distance-removal field, the budget set
#'   in [fit_gdr_response_field()] also applies.
#' @return A `multiScaleR_field_geometry` object, a small list that stores no
#'   model, coordinates or predictions. It contains:
#'   * `schema_version`, and for distance-removal fields `provider`
#'     (`"native_GDR_fixed_range"`), `range` and `range_estimated` (`FALSE`).
#'   * `crs` (WKT text) and `units`, the linear units of the CRS.
#'   * For `glmmTMB` fields, `position` and `group`, the names of the position
#'     and group columns, with the fitted `group_value` and `group_levels`.
#'   * `training_locations` (distinct field locations) and `training_rows`.
#'   * `provider_version`, the `glmmTMB` or `TMB` version.
#'   * `batch_size` and `max_memory_mb`, the prediction controls.
#'   * `reference_sha256`, `binding_sha256` and `definition_sha256`, SHA-256
#'     fingerprints of the coordinates, of the fitted context and of the
#'     descriptor itself. Prediction stops if any of them no longer matches,
#'     so recreate the descriptor after refitting or changing its inputs.
#' @details
#' # Map a fitted field
#'
#' Pass the descriptor to [predict.multiScaleR()] through
#' `raster_args = list(field_geometry = geometry)`. Prepare the layers with
#' [kernel_scale.raster()] using `scale_center = TRUE` and
#' `match_buffer = TRUE`. Each cell then supplies its own projected
#' coordinates to the field, and prediction runs on one core, in batches of
#' `batch_size` cells. No covariance among cells is computed or stored.
#' Several rows in `raster_args$const` give one band (or fit and SE pair) per
#' scenario row, in row order.
#'
#' # glmmTMB fields
#'
#' Every raster cell is a new field position, so set
#' `native_args = list(allow.new.levels = TRUE)` in [predict.multiScaleR()].
#' Without it, prediction stops before the first batch. All cells use the
#' single fitted field group, and a scenario that names another group stops.
#' Supply ordinary site and offset scenarios as usual. A position or group
#' scenario passed to [kernel_scale.raster()] is replaced by the cell
#' coordinates during prediction. Do not substitute a constant coordinate.
#' On the response scale, predictions are the expected count for Poisson and
#' NB2 models, the mean for Gaussian and Gamma models, and the per-trial
#' probability for binomial models, including binomial models fitted to
#' success and failure counts.
#'
#' # Distance-removal fields
#'
#' Use `type = "lambda"` for density. Each cell gets the field value predicted
#' at its center from the fitted block values, with the range held fixed.
#' `allow.new.levels` is not used. See [fit_gdr_response_field()] for the
#' full workflow, including scale-of-effect estimation with the field.
#'
#' # What the standard errors include
#'
#' Standard errors are conditional on the selected kernel scales. They include
#' the fitted coefficient and field uncertainty, and for distance-removal
#' fields the extra field variance between blocks. They add no observation
#' noise and no kernel-scale uncertainty. `re.form = NA` gives a zero-field
#' reference, not an average over the field. Checking the coordinates does
#' not establish that the field is well estimated or that its intervals are
#' calibrated.
#' @seealso [predict.multiScaleR()], [kernel_scale.raster()],
#'   [fit_gdr_response_field()], [prediction_contract()],
#'   `vignette("v11-spatial-fields", package = "multiScaleR")`
#' @examples
#' \donttest{
#' if (requireNamespace("glmmTMB", quietly = TRUE)) {
#'   # Give the bundled landscape (in meters) a projected CRS
#'   land <- terra::rast(system.file("extdata", "landscape.tif",
#'                                   package = "multiScaleR"))[["land2"]]
#'   terra::crs(land) <- "EPSG:32617"
#'   data("surv_pts")
#'   pts <- terra::vect(surv_pts)
#'   terra::crs(pts) <- "EPSG:32617"
#'   pts <- pts[seq(1, 100, by = 2)]   # 50 survey points
#'   xy <- terra::crds(pts)
#'
#'   # Counts depend on habitat at a 400 m scale plus a correlated field
#'   set.seed(3)
#'   truth <- kernel_scale.raster(land, sigma = 400, kernel = "gaussian",
#'                                verbose = FALSE)
#'   x2 <- as.numeric(scale(terra::extract(truth, pts, ID = FALSE)[, 1]))
#'   field <- drop(t(chol(0.25 * exp(-as.matrix(stats::dist(xy)) / 800))) %*%
#'                   stats::rnorm(nrow(xy)))
#'   inputs <- kernel_prep(pts = pts, raster_stack = land, max_D = 1200,
#'                         kernel = "gaussian", verbose = FALSE)
#'   df <- data.frame(counts = stats::rpois(nrow(xy), exp(1 + 0.6 * x2 + field)),
#'                    inputs$kernel_dat)
#'   df$pos <- glmmTMB::numFactor(xy[, 1], xy[, 2])   # field positions (m)
#'   df$grp <- factor("all")                          # one shared field
#'   m0 <- glmmTMB::glmmTMB(counts ~ land2 + exp(pos + 0 | grp),
#'                          family = poisson(), data = df,
#'                          start = list(theta = c(log(0.5), log(500))))
#'   fit <- multiScale_optim(m0, inputs, n_cores = 1, verbose = FALSE)
#'
#'   # Bind the fitting coordinates and their CRS to the fitted field
#'   geometry <- response_field_geometry(fit, coords = xy, crs = 32617)
#'   print(geometry$training_locations)   # 50 distinct field locations
#'
#'   # Expected counts with conditional SEs for a 10 x 10 cell window
#'   layers <- kernel_scale.raster(land, multiScaleR = fit,
#'     scale_center = TRUE, match_buffer = TRUE, clamp = TRUE, verbose = FALSE,
#'     site_values = list(pos = as.character(df$pos[1]), grp = "all"))
#'   window <- terra::crop(layers, terra::ext(1500, 1700, 1500, 1700))
#'   attr(window, "multiScaleR_projection") <-
#'     attr(layers, "multiScaleR_projection")
#'   map <- predict(fit, window, se.fit = TRUE, uncertainty = "conditional",
#'                  native_args = list(allow.new.levels = TRUE),
#'                  raster_args = list(field_geometry = geometry))
#'   print(map)   # bands: fit (expected count) and se.fit
#' }
#' }
#' @export
response_field_geometry <- function(object, coords, crs, batch_size = 128L,
                                    max_memory_mb = 512) {
  if (inherits(object, "multiScaleR") && methods::is(object$opt_mod, "unmarkedFitGDRField"))
    return(.msr_gdr_field_geometry(object, coords, crs, batch_size, max_memory_mb))
  if (!inherits(object, "multiScaleR") || !.msr_refit_backend_fit(object) ||
      !inherits(.analysis_model(object$opt_mod), "glmmTMB"))
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
  if (methods::is(model, "unmarkedFitGDRField"))
    return(.msr_gdr_field_prediction_budget(geometry, model, rows))
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
  if (identical(geometry$provider, "native_GDR_fixed_range"))
    return(.msr_gdr_field_raster_eval(model, data))
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
