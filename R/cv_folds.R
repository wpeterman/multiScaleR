# Deterministic balancing keeps indivisible sampling units intact. Stable unit
# labels resolve equal sizes without consuming or changing the user's RNG.
.msr_cv_assign <- function(unit, k) {
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k < 2 || k != floor(k))
    stop("k must be an integer of at least two.", call. = FALSE)
  labels <- sort(unique(unit))
  if (length(labels) < k)
    stop("There must be at least k indivisible groups or spatial units.", call. = FALSE)
  size <- tabulate(match(unit, labels), length(labels))
  allocation <- integer(length(labels)); load <- numeric(k)
  for (i in order(-size, labels)) {
    f <- which.min(load)
    allocation[i] <- f; load[f] <- load[f] + size[i]
  }
  allocation[match(unit, labels)]
}

.msr_fold_ids <- function(ids, n) {
  if (is.null(ids)) ids <- as.character(seq_len(n))
  if (!.msr_valid_observation_ids(ids, n))
    stop("ids must be unique, nonmissing, nonempty character row IDs.", call. = FALSE)
  ids
}

.msr_fold_group <- function(group, ids) {
  if ((!is.atomic(group) && !is.factor(group)) || !is.null(dim(group)) ||
      length(group) != length(ids) || anyNA(group) ||
      (is.numeric(group) && any(!is.finite(group))) ||
      (!is.null(names(group)) && !identical(names(group), ids)))
    stop("group must have one nonmissing label per row, with names aligned to ids.", call. = FALSE)
  group <- as.character(group)
  if (any(!nzchar(group))) stop("group labels must not be empty.", call. = FALSE)
  group
}

.msr_make_folds <- function(assignment, ids, info, excluded = NULL) {
  labels <- paste0("fold", seq_len(max(assignment)))
  folds <- stats::setNames(lapply(seq_along(labels), function(f) {
    omit <- if (is.null(excluded)) character() else excluded[[f]]
    list(analysis = ids[assignment != f & !ids %in% omit], assessment = ids[assignment == f])
  }), labels)
  .msr_cv_folds(folds, ids)
  info$excluded <- if (is.null(excluded)) stats::setNames(rep(list(character()), length(labels)), labels) else
    stats::setNames(excluded, labels)
  attr(folds, "fold_info") <- info
  class(folds) <- c("multiScaleR_folds", "list")
  folds
}

#' Create Grouped Cross-Validation Folds
#'
#' Keep all rows from a sampling group in the same assessment fold. Whole
#' groups are assigned to balance row counts, largest groups first.
#'
#' @param group Vector of sampling-group labels, one per row. Animals, sites or
#'   matched choice sets can define groups. Missing or empty labels are rejected.
#' @param k Number of folds, at least two and no greater than the number of
#'   indivisible groups. The default is five.
#' @param ids Unique character row IDs, in the order of \code{group}. Defaults
#'   to its names, or character row numbers when it is unnamed. Named groups
#'   must align exactly with supplied IDs.
#' @return A named list of class \code{multiScaleR_folds}. Each element contains
#'   disjoint character vectors \code{analysis} and \code{assessment}, accepted
#'   by \code{\link{cross_validate}}. Every row is assessed once. Attribute
#'   \code{fold_info} records the method, row/unit assignments and exclusions.
#' @details Assignment is deterministic and does not use the random-number
#'   generator. Equal-size groups are ordered by their labels. Balancing counts
#'   does not stratify responses. Choose sampling groups and the validation
#'   target before examining assessment outcomes. Fold construction does not
#'   qualify a model class for scoring; matched-set survival/amt scoring still
#'   requires a separate adapter.
#' @examples
#' ids <- paste0("row", 1:12)
#' folds <- grouped_cv_folds(rep(letters[1:6], each = 2), k = 3, ids = ids)
#' folds
#' @seealso \code{\link{spatial_cv_folds}}, \code{\link{cross_validate}},
#'   \code{\link{cv_tab}}
#' @export
grouped_cv_folds <- function(group, k = 5, ids = NULL) {
  if (is.null(ids)) ids <- names(group)
  ids <- .msr_fold_ids(ids, length(group))
  group <- .msr_fold_group(group, ids)
  assignment <- .msr_cv_assign(group, k)
  .msr_make_folds(assignment, ids, list(method = "grouped", k = as.integer(k),
    assignment = data.frame(row_id = ids, group = group, unit = group, fold = paste0("fold", assignment))))
}

# Spatial group/block intersections form connected components. Union-find uses
# storage proportional to rows and blocks rather than a pairwise distance matrix.
.msr_cv_spatial_units <- function(block, group) {
  n <- max(block); parent <- seq_len(n)
  root <- function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i <- parent[i]
    }
    i
  }
  if (!is.null(group)) for (cells in split(block, group)) {
    cells <- unique(cells)
    a <- root(cells[1L])
    for (b in cells[-1L]) {
      b <- root(b)
      if (a != b) {
        parent[max(a, b)] <- min(a, b)
        a <- min(a, b)
      }
    }
  }
  roots <- vapply(seq_len(n), root, integer(1))
  roots[block]
}

.msr_cv_coordinates <- function(coords, crs, ids) {
  if (inherits(coords, "SpatVector")) coords <- sf::st_as_sf(coords)
  if (inherits(coords, "sfc")) coords <- sf::st_sf(geometry = coords)
  if (inherits(coords, "sf")) {
    if (!all(sf::st_geometry_type(coords) == "POINT") || any(sf::st_is_empty(coords)))
      stop("coords must contain nonempty POINT geometries.", call. = FALSE)
    native_crs <- sf::st_crs(coords)
    if (!is.null(crs) && !isTRUE(native_crs == sf::st_crs(crs)))
      stop("crs must match the existing coordinate CRS; transform coordinates explicitly.", call. = FALSE)
    crs <- native_crs
    xy <- sf::st_coordinates(coords)
    rows <- rownames(coords)
  } else {
    if ((!is.matrix(coords) && !is.data.frame(coords)) || ncol(coords) != 2L)
      stop("coords must be projected points or a numeric two-column coordinate table.", call. = FALSE)
    xy <- as.matrix(coords); rows <- rownames(coords)
    if (is.null(crs)) stop("Coordinate tables require an explicit projected crs.", call. = FALSE)
    crs <- sf::st_crs(crs)
  }
  if (is.na(crs) || isTRUE(sf::st_is_longlat(crs)))
    stop("Spatial folds require a known projected CRS, not longitude/latitude.", call. = FALSE)
  if (!is.numeric(xy) || ncol(xy) != 2L || !nrow(xy) || any(!is.finite(xy)))
    stop("coords must contain finite two-dimensional coordinates.", call. = FALSE)
  if (is.null(ids) && !is.null(rows)) ids <- rows
  ids <- .msr_fold_ids(ids, nrow(xy))
  if (!is.null(rows) && !identical(rows, as.character(seq_len(nrow(xy)))) && !identical(rows, ids))
    stop("Coordinate row names must align with ids.", call. = FALSE)
  list(xy = xy, crs = crs, ids = ids)
}

#' Create Spatially Blocked Cross-Validation Folds
#'
#' Assign whole projected grid blocks to assessment folds. Optional sampling
#' groups join intersecting blocks into indivisible units. A distance buffer
#' can exclude nearby training units without dropping assessment rows.
#'
#' @param coords Projected \code{sf} POINTs, a \code{terra::SpatVector} of points,
#'   or a numeric two-column coordinate matrix/data frame. Row order must match
#'   the data. Meaningful row names must match \code{ids}.
#' @param block_size Positive grid-cell width, in coordinate CRS units.
#' @param k Number of folds, at least two. The default is five. At least
#'   \code{k} indivisible spatial units must remain after joining groups.
#' @param ids Unique character row IDs. Defaults to coordinate row names, or
#'   character row numbers. Use the same IDs as the modeling data.
#' @param group Optional sampling-group labels aligned with \code{ids}.
#'   Blocks joined by any shared group are kept together, including transitive
#'   connections. This can reduce the number of available units.
#' @param buffer Nonnegative training-to-assessment separation in coordinate
#'   CRS units. Default zero adds no exclusions. With a positive buffer, an
#'   entire training unit is excluded if any of its points is within or exactly
#'   at this distance from an assessment point. Empty training sets are errors.
#' @param origin Numeric grid origin of length two, in coordinate CRS units.
#'   Defaults to \code{c(0, 0)}. Cells use half-open intervals of width
#'   \code{block_size}, indexed by \code{floor((coordinate - origin)/block_size)}.
#' @param crs Projected CRS declaration for coordinate tables, accepted by
#'   \code{sf::st_crs()}. For spatial objects, an optional declaration must match
#'   their existing CRS. Coordinates are never transformed or relabeled here.
#' @return A \code{multiScaleR_folds} list with analysis/assessment IDs and
#'   \code{fold_info} metadata: CRS, units, grid/buffer controls, compact row,
#'   block, group and fold assignments, and excluded training IDs per fold.
#' @details Units are assigned largest first to balance row counts, with stable
#'   ties and no random draws. Several disjoint blocks can share a fold; this
#'   does not produce contiguous geographic regions. Buffers measure separation
#'   between sampled points, not distances to entire cell boundaries. Whole-unit
#'   exclusions can create a larger gap than requested. Distance calculations
#'   use nearest-point indices and elementwise distances rather than a dense
#'   training-by-assessment distance matrix.
#'
#' Choose block size, groups, buffer and the prediction target from the sampling
#'   design and ecological question. A projected CRS and a distance buffer alone
#'   do not establish independence or validate calibration. All rows are assessed
#'   once; some can be absent from a training set because of the buffer.
#' @examples
#' xy <- cbind(x = rep(c(0, 100, 200), each = 2), y = rep(c(0, 10), 3))
#' spatial_cv_folds(xy, block_size = 50, k = 3, crs = 32617)
#' @seealso \code{\link{grouped_cv_folds}}, \code{\link{cross_validate}},
#'   \code{\link{cv_tab}}
#' @export
spatial_cv_folds <- function(coords, block_size, k = 5, ids = NULL,
    group = NULL, buffer = 0, origin = c(0, 0), crs = NULL) {
  c <- .msr_cv_coordinates(coords, crs, ids); ids <- c$ids
  valid <- function(x) is.numeric(x) && length(x) == 1L && is.finite(x)
  if (!valid(block_size) || block_size <= 0 || !valid(buffer) || buffer < 0 ||
      !is.numeric(origin) || length(origin) != 2L || any(!is.finite(origin)))
    stop("Use positive finite block_size, nonnegative finite buffer and a finite two-coordinate origin.", call. = FALSE)
  grid <- floor(sweep(c$xy, 2L, origin) / block_size)
  if (any(!is.finite(grid) | abs(grid) > 2^52))
    stop("Grid indices exceed exact representable coordinates; revise origin or block_size.", call. = FALSE)
  cells <- unique(as.data.frame(grid))
  cells <- cells[order(cells[[1L]], cells[[2L]]), , drop = FALSE]
  key <- function(x) paste(formatC(x[, 1L], format = "f", digits = 0),
    formatC(x[, 2L], format = "f", digits = 0), sep = ":")
  block <- match(key(grid), key(as.matrix(cells)))
  if (!is.null(group)) group <- .msr_fold_group(group, ids)
  unit <- .msr_cv_spatial_units(block, group)
  assignment <- .msr_cv_assign(unit, k)
  excluded <- rep(list(character()), k)
  if (buffer > 0) {
    points <- sf::st_as_sf(data.frame(x = c$xy[, 1L], y = c$xy[, 2L]), coords = c("x", "y"), crs = c$crs)
    for (f in seq_len(k)) {
      train <- which(assignment != f); test <- which(assignment == f)
      nearest <- sf::st_nearest_feature(points[train, ], points[test, ], longlat = FALSE)
      distances <- as.numeric(sf::st_distance(points[train, ], points[test[nearest], ], by_element = TRUE))
      if (any(!is.finite(distances))) stop("Could not determine all buffer distances.", call. = FALSE)
      excluded[[f]] <- ids[assignment != f & unit %in% unit[train[distances <= buffer]]]
      if (length(excluded[[f]]) == length(train))
        stop("The buffer removes every training unit in fold", f, ". Reduce k/buffer or revise the design.", call. = FALSE)
    }
  }
  table <- data.frame(row_id = ids, x = c$xy[, 1L], y = c$xy[, 2L], block = block,
    unit = unit, fold = paste0("fold", assignment))
  if (!is.null(group)) table$group <- group
  .msr_make_folds(assignment, ids, list(method = "spatial", k = as.integer(k),
    block_size = block_size, buffer = buffer, origin = origin, crs = c$crs$wkt,
    units = c$crs$units_gdal, assignment = table), excluded)
}

#' @rdname grouped_cv_folds
#' @param x A \code{multiScaleR_folds} object.
#' @param ... Reserved for future print controls.
#' @export
print.multiScaleR_folds <- function(x, ...) {
  info <- attr(x, "fold_info", exact = TRUE)
  cat("Cross-validation folds:", info$method, "\n")
  print(data.frame(fold = names(x), analysis = lengths(lapply(x, `[[`, "analysis")),
    assessment = lengths(lapply(x, `[[`, "assessment")), excluded = lengths(info$excluded)), row.names = FALSE)
  invisible(x)
}
