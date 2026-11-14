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
#' Keep all rows from a sampling group in the same assessment fold, so that a
#' held-out row is never predicted with help from its own group. Whole groups
#' are assigned to folds to balance the row counts, largest groups first.
#'
#' @param group Vector of sampling-group labels, one per row, such as site,
#'   animal, year, or choice-set IDs. Missing or empty labels are rejected.
#' @param k Number of folds, at least two and no greater than the number of
#'   groups. The default is five.
#' @param ids Unique character row IDs, in the order of `group`. Use the row
#'   names of the modeling data. The default uses the names of `group`, or the
#'   row numbers as character strings when `group` is unnamed. Named groups
#'   must align exactly with supplied IDs.
#' @return A named list of class `multiScaleR_folds` with one element per fold
#'   (`fold1`, `fold2`, and so on), ready for [cross_validate()] and the other
#'   cross-validation functions. Each element holds two character vectors of
#'   row IDs, `analysis` (rows used to fit) and `assessment` (rows held out).
#'   Every row is assessed exactly once. The attribute `fold_info` holds the
#'   `method` (`"grouped"`), `k`, an `assignment` data frame with the `row_id`,
#'   `group`, `unit`, and assessment `fold` of every row, and `excluded`, the
#'   IDs left out of each fold's training rows (none for grouped folds).
#'   `print()` shows the number of analysis, assessment, and excluded rows in
#'   each fold, and returns `x` invisibly.
#' @details Assignment is deterministic and does not use the random-number
#'   generator, so the same call always gives the same folds. Groups of equal
#'   size are ordered by their labels. Balancing row counts does not stratify
#'   the responses. Choose the sampling groups and the prediction question
#'   before you look at any scores. To keep groups together within spatial
#'   blocks, pass them as `group` to [spatial_cv_folds()]. Building folds does
#'   not decide which models can be scored. For matched choice sets, score
#'   the folds with [cross_validate_clogit()].
#'   `vignette("v08-cross-validation", package = "multiScaleR")` compares
#'   grouped, spatial, and random fold designs.
#' @examples
#' # Twelve records from six sites, two visits each
#' ids <- paste0("row", 1:12)
#' folds <- grouped_cv_folds(rep(letters[1:6], each = 2), k = 3, ids = ids)
#' folds                    # rows fitted and held out in each fold
#' folds$fold1$assessment   # both visits to a site are held out together
#' @seealso [spatial_cv_folds()], [cross_validate()], [cv_tab()]
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
#' Lay a square grid over projected points and assign whole grid blocks to
#' assessment folds, so that the model must predict parts of the landscape it
#' was not fitted to. Optional sampling groups join the blocks they touch into
#' indivisible units. A distance buffer can exclude nearby training units
#' without dropping assessment rows.
#'
#' @param coords Projected `sf` POINTs, a `terra::SpatVector` of points, or a
#'   numeric two-column coordinate matrix or data frame. Row order must match
#'   the data. Meaningful row names must match `ids`.
#' @param block_size Positive width of each grid cell, in the units of the
#'   coordinate CRS (usually meters). Larger blocks put more distance between
#'   held-out points and their nearest training points, which moves the task
#'   closer to predicting a new area.
#' @param k Number of folds, at least two. The default is five. At least `k`
#'   indivisible spatial units must remain after joining groups.
#' @param ids Unique character row IDs. Defaults to the row names of `coords`,
#'   or to the row numbers as character strings. Use the same IDs as the
#'   modeling data.
#' @param group Optional sampling-group labels aligned with `ids`, such as site
#'   or animal IDs. Blocks joined by any shared group are kept together,
#'   including through chains of shared groups. This can reduce the number of
#'   available units.
#' @param buffer Nonnegative training-to-assessment separation, in coordinate
#'   CRS units. The default 0 adds no exclusions. With a positive buffer, an
#'   entire training unit is excluded from a fold if any of its points lies
#'   within, or exactly at, this distance from an assessment point of that
#'   fold. Assessment rows are never dropped. A buffer that leaves a fold with
#'   no training rows is an error.
#' @param origin Numeric grid origin of length two, in coordinate CRS units.
#'   Defaults to `c(0, 0)`. Cells are half-open intervals of width
#'   `block_size`, indexed by `floor((coordinate - origin) / block_size)`.
#' @param crs Projected CRS declaration for coordinate tables, in any form
#'   accepted by `sf::st_crs()`. For spatial objects, an optional declaration
#'   must match their existing CRS. Coordinates are never transformed or
#'   relabeled here, and longitude and latitude are rejected.
#' @return A named list of class `multiScaleR_folds`, structured as described
#'   in [grouped_cv_folds()]. Each fold holds `analysis` and `assessment` row
#'   IDs. Every row is assessed exactly once, although a buffer can leave a
#'   row out of some training sets. The attribute `fold_info` holds the
#'   `method` (`"spatial"`), `k`, `block_size`, `buffer`, `origin`, the `crs`
#'   (as WKT) and its distance `units`, an `assignment` data frame with each
#'   row's `row_id`, coordinates `x` and `y`, grid `block`, joined `unit`,
#'   assessment `fold`, and `group` when supplied, and `excluded`, the
#'   training IDs removed by the buffer in each fold.
#' @details Units are assigned largest first to balance row counts, with
#'   stable tie-breaking and no random draws, so the same call always gives
#'   the same folds. Several disjoint blocks can share a fold, so a fold need
#'   not be one contiguous region. Buffers measure the separation between
#'   sampled points, not the distance to block edges. Because whole units are
#'   excluded, the gap can be larger than requested, and a large buffer can
#'   leave few training points. Distance calculations use nearest-point
#'   searches rather than a dense training-by-assessment distance matrix.
#'
#'   Choose the block size, groups, buffer, and prediction target from the
#'   sampling design and the ecological question. A projected CRS and a
#'   distance buffer alone do not establish independence or validate
#'   calibration. `vignette("v08-cross-validation", package = "multiScaleR")`
#'   maps a set of spatial folds and shows how a buffer changes the training
#'   sets.
#' @examples
#' # Six points in three 50 m blocks (crs = 32617 declares meter units)
#' xy <- cbind(x = rep(c(0, 100, 200), each = 2), y = rep(c(0, 10), 3))
#' spatial_cv_folds(xy, block_size = 50, k = 3, crs = 32617)
#'
#' # The bundled survey points, in 1,000 m blocks. A 320 m buffer removes
#' # training blocks that have a point near any held-out point.
#' data("surv_pts")
#' pts <- crds(vect(surv_pts))
#' buffered <- spatial_cv_folds(pts, block_size = 1000, k = 3, buffer = 320,
#'                              crs = 32617)
#' buffered   # 'excluded' counts rows used for neither fitting nor scoring
#' head(attr(buffered, "fold_info")$assignment)   # block and fold of each row
#' @seealso [grouped_cv_folds()], [cross_validate()], [cv_tab()]
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
#' @param x A `multiScaleR_folds` object, for `print()`.
#' @param ... Not used. Reserved for future print controls.
#' @export
print.multiScaleR_folds <- function(x, ...) {
  info <- attr(x, "fold_info", exact = TRUE)
  cat("Cross-validation folds:", info$method, "\n")
  print(data.frame(fold = names(x), analysis = lengths(lapply(x, `[[`, "analysis")),
    assessment = lengths(lapply(x, `[[`, "assessment")), excluded = lengths(info$excluded)), row.names = FALSE)
  invisible(x)
}
