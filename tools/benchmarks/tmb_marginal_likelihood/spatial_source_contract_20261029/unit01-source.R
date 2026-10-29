#' Check Sources for a Private Spatial Scale Workflow
#'
#' Compare explicit point and raster sources with prepared point IDs and the
#' numeric coordinate columns used by a spatial scale model. This optional
#' preflight describes source agreement without changing a fitted model.
#'
#' @param points An \code{sf} or terra \code{SpatVector} of nonempty projected
#'   two-dimensional points with a defined CRS (coordinate reference system).
#' @param raster One terra \code{SpatRaster} in the same projected CRS.
#' @param inputs Binned \code{multiScaleR_data} inputs with unique point IDs.
#'   These compact summaries do not prove historical source or distance units.
#' @param data Observation data frame. Repeated source-point IDs are allowed.
#' @param coordinates Names of two finite numeric columns in \code{data}, in
#'   the same order and map units as the source geometry coordinates.
#' @param row_id Optional observation ID column. The default uses row names.
#' @param point_id Optional source-point ID column. The default uses row names.
#' @param coordinate_tolerance Nonnegative absolute coordinate tolerance in
#'   source map units. Zero, the default, requires exact agreement. The report
#'   records any requested tolerance; coordinates are never transformed or scaled.
#' @return A private \code{msr_spatial_source_report} list containing:
#'   \describe{
#'     \item{source}{Supplied CRS WKT and linear unit, coordinate names, ID
#'       columns, raster layer names and grid extent/resolution/dimensions. Grid
#'       distances use source map units. Historical source lineage and extraction
#'       unit agreement remain explicitly unverified.}
#'     \item{alignment}{Observation/source-point counts, the recorded absolute
#'       tolerance and maximum coordinate difference in source map units.}
#'     \item{extent}{Nominal buffer radius \code{inputs$max_D}, its unverified
#'       unit assumption and a table in binned point-ID order. The table gives
#'       each point's minimum signed clearance from the raster extent, whether
#'       the point lies inside it and whether its circular buffer fits inside it.
#'       Negative clearance means a point lies outside the supplied extent.
#'       Actual valid-cell and kernel-weight coverage are not evaluated.}
#'     \item{distance_scaling}{The existing internal distance scaling factor
#'       \code{unit_conv}. This is not a conversion between meters and feet.}
#'     \item{fit_eligibility, calibration}{Always false and
#'       \code{"not_established"}, respectively. Source agreement supplies no
#'       inference acceptance or uncertainty calibration.}
#'   }
#' @details Missing/geographic/different CRSs, unsupported geometry, ambiguous
#'   IDs, missing raster source layers or disagreeing coordinate columns stop
#'   this optional report. It does not run automatically during preparation,
#'   fitting or prediction. It reads grid metadata without extracting raster
#'   values, changing prepared summaries or refitting either inferential target.
#'
#'   Read \code{extent$points$buffer_within_extent} as a geometric check only.
#'   A true value does not rule out raster holes, masks or incomplete support
#'   for a particular covariate. A false value flags a potentially clipped
#'   extraction buffer; it does not reject a fitted model. The nominal radius
#'   calculation assumes \code{max_D} uses the supplied CRS linear unit. Compact
#'   inputs cannot certify that historical assumption or raster-value lineage.
#'
#'   The spatial field's covariance range and the ecological kernel scale are
#'   distinct distances. This report estimates neither and changes no priors,
#'   likelihood/posterior target, fitting gates or prediction uncertainty.
#' @seealso \code{\link{kernel_prep}}, \code{\link{kernel_prep_by_group}}
#' @keywords internal
.msr_spatial_source_report <- function(points, raster, inputs, data, coordinates,
    row_id = NULL, point_id = NULL, coordinate_tolerance = 0) {
  validate_scalar_numeric(coordinate_tolerance, "coordinate_tolerance", lower = 0)
  if (!inherits(points, c("sf", "SpatVector")))
    stop("points must be sf or terra SpatVector projected POINT geometry.", call. = FALSE)
  points <- sf::st_as_sf(points)
  if (!nrow(points) || any(sf::st_is_empty(points)) ||
      !all(sf::st_geometry_type(points) == "POINT") || is.na(sf::st_crs(points)) ||
      !isFALSE(sf::st_is_longlat(points)))
    stop("points require nonempty projected POINT geometry and a defined CRS.", call. = FALSE)
  xy <- sf::st_coordinates(points)
  if (ncol(xy) != 2L || any(!is.finite(xy)))
    stop("points require finite two-dimensional coordinates.", call. = FALSE)
  if (!inherits(raster, "SpatRaster") || terra::nlyr(raster) < 1L)
    stop("raster must be one SpatRaster with at least one layer.", call. = FALSE)
  crs <- sf::st_crs(points)
  if (!isTRUE(crs == sf::st_crs(terra::crs(raster))))
    stop("Points and raster must share the same defined projected CRS.", call. = FALSE)
  linear_unit <- crs$units_gdal
  if (length(linear_unit) != 1L || is.na(linear_unit) || !nzchar(linear_unit))
    stop("The source CRS must declare a linear unit.", call. = FALSE)
  if (!inherits(inputs, "multiScaleR_data") || is.null(inputs$binned))
    stop("inputs must be binned multiScaleR_data.", call. = FALSE)
  for (nm in c("max_D", "unit_conv")) validate_scalar_numeric(inputs[[nm]], nm, positive = TRUE)
  if (!is.data.frame(data) || !nrow(data) || anyDuplicated(names(data)))
    stop("data must contain observation rows and distinct columns.", call. = FALSE)
  if (inherits(data, "sf")) {
    if (!isTRUE(sf::st_crs(data) == crs))
      stop("Spatial observation data must share the source CRS.", call. = FALSE)
    data <- sf::st_drop_geometry(data)
  }
  observed_xy <- .msr_spatial_coordinates(data, coordinates)
  # Explicit keys align rows without assuming either source or observation order.
  keys <- function(frame, column, label) {
    if (!is.null(column) && (!is.character(column) || length(column) != 1L ||
        is.na(column) || !column %in% names(frame)))
      stop(label, " must name one ID column.", call. = FALSE)
    value <- if (is.null(column)) rownames(frame) else frame[[column]]
    if (!is.atomic(value) || !is.null(dim(value)) || length(value) != nrow(frame) ||
        anyNA(value) || any(!nzchar(as.character(value))))
      stop(label, " requires one nonmissing nonempty ID per row.", call. = FALSE)
    as.character(value)
  }
  ids <- inputs$binned$point_ids
  if (!is.atomic(ids) || !is.null(dim(ids)) || !length(ids) || anyNA(ids) ||
      any(!nzchar(as.character(ids))) || anyDuplicated(ids) ||
      !length(inputs$binned$vsum) || any(vapply(inputs$binned$vsum,
        function(x) !is.matrix(x) || nrow(x) != length(ids), TRUE)))
    stop("Binned point IDs must uniquely align every prepared matrix row.", call. = FALSE)
  ids <- as.character(ids)
  source_ids <- keys(points, point_id, "point_id")
  if (anyDuplicated(source_ids) || length(source_ids) != length(ids) || !setequal(source_ids, ids))
    stop("Source point IDs must uniquely match all binned point IDs.", call. = FALSE)
  observation_ids <- keys(data, row_id, "row_id")
  rows <- match(observation_ids, ids)
  if (anyNA(rows)) stop("Observation IDs do not match binned point IDs.", call. = FALSE)
  xy <- xy[match(ids, source_ids), , drop = FALSE]
  coordinate_error <- max(abs(observed_xy - xy[rows, , drop = FALSE]))
  if (coordinate_error > coordinate_tolerance)
    stop("Coordinate columns disagree with supplied point geometry in source map units (maximum difference ",
      format(coordinate_error), "; tolerance ", format(coordinate_tolerance), ").", call. = FALSE)
  sources <- unique(as.character(inputs$scale_vars$source))
  if (anyNA(sources) || any(!sources %in% names(raster)))
    stop("The supplied raster is missing prepared source layers.", call. = FALSE)
  bounds <- as.vector(terra::ext(raster))
  clearance <- pmin(xy[, 1L] - bounds[1L], bounds[2L] - xy[, 1L],
    xy[, 2L] - bounds[3L], bounds[4L] - xy[, 2L])
  structure(list(source = list(crs_wkt = crs$wkt, linear_unit = linear_unit,
    coordinates = coordinates, point_id = point_id, row_id = row_id,
    raster_layers = names(raster), prepared_source_layers = sources,
    grid = list(extent = stats::setNames(bounds, c("xmin", "xmax", "ymin", "ymax")),
      resolution = terra::res(raster), dimensions = c(rows = terra::nrow(raster), columns = terra::ncol(raster))),
    historical_raster_lineage = "unverified_from_compact_inputs",
    historical_distance_unit_agreement = "unverified_from_compact_inputs"),
    alignment = list(observations = nrow(data), source_points = length(ids),
      distinct_observed_points = length(unique(observation_ids)),
      coordinate_tolerance = coordinate_tolerance, maximum_coordinate_difference = coordinate_error),
    extent = list(nominal_buffer_radius = inputs$max_D,
      unit_assumption = "max_D is expressed in the supplied projected CRS linear unit; historical agreement is unverified",
      points = data.frame(point_id = ids, minimum_clearance = unname(clearance),
        point_inside_extent = unname(clearance >= 0), buffer_within_extent = unname(clearance >= inputs$max_D)),
      valid_cell_coverage = "not_evaluated", kernel_weight_coverage = "not_evaluated"),
    distance_scaling = list(value = inputs$unit_conv,
      meaning = "internal distance scaling, not a conversion between CRS linear units"),
    fit_eligibility = FALSE, calibration = "not_established"), class = "msr_spatial_source_report")
}
