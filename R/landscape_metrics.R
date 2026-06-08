# Internal helpers for exploratory landscape metric support.

.landscape_adjacency_class_ceiling <- function(metric) {
  switch(
    metric,
    contag = 50L,
    ai = 200L,
    pladj = 200L,
    200L
  )
}


.landscape_validate_categorical_values <- function(values,
                                                   metric,
                                                   max_classes = NULL,
                                                   context = "Landscape metric") {
  finite_values <- as.vector(values)
  finite_values <- finite_values[!is.na(finite_values)]

  if (length(finite_values) == 0) {
    return(invisible(integer(0)))
  }

  if (any(!is.finite(finite_values))) {
    stop(
      sprintf("%s `%s` requires finite categorical class values.", context, metric),
      call. = FALSE
    )
  }

  tolerance <- sqrt(.Machine$double.eps)
  if (any(abs(finite_values - round(finite_values)) > tolerance)) {
    examples <- unique(finite_values[abs(finite_values - round(finite_values)) > tolerance])
    examples <- utils::head(examples, 5)
    stop(
      sprintf(
        "%s `%s` requires a categorical raster encoded with integer-like class values; found non-integer values such as %s.",
        context,
        metric,
        paste(signif(examples, 6), collapse = ", ")
      ),
      call. = FALSE
    )
  }

  classes <- sort(unique(round(finite_values)))
  n_classes <- length(classes)
  if (!is.null(max_classes) && n_classes > max_classes) {
    stop(
      sprintf(
        "%s `%s` has %s classes, which exceeds the current supported ceiling of %s classes. Reclassify or aggregate categories before using this adjacency metric.",
        context,
        metric,
        n_classes,
        max_classes
      ),
      call. = FALSE
    )
  }

  invisible(classes)
}


.landscape_clean_count_matrix <- function(x, tolerance = 1e-8) {
  x[abs(x) < tolerance] <- 0
  rounded <- round(x)
  near_integer <- is.finite(x) & abs(x - rounded) < tolerance
  x[near_integer] <- rounded[near_integer]
  x
}


.landscape_class_totals <- function(values, weights = NULL) {
  .landscape_validate_categorical_values(values, metric = "composition")

  keep <- !is.na(values)
  values <- values[keep]

  if (is.null(weights)) {
    weights <- rep(1, length(values))
  } else {
    weights <- weights[keep]
    valid_weights <- !is.na(weights) & weights > 0
    values <- values[valid_weights]
    weights <- weights[valid_weights]
  }

  if (length(values) == 0 || sum(weights) <= 0) {
    return(numeric(0))
  }

  totals <- rowsum(weights, group = values, reorder = FALSE)[, 1]
  class_names <- names(totals)
  numeric_names <- suppressWarnings(as.numeric(class_names))
  order_idx <- if (all(!is.na(numeric_names))) {
    order(numeric_names)
  } else {
    order(class_names)
  }

  totals[order_idx]
}


.landscape_class_proportions <- function(values, weights = NULL) {
  totals <- .landscape_class_totals(values = values, weights = weights)
  if (length(totals) == 0) {
    return(numeric(0))
  }

  totals / sum(totals)
}


.landscape_composition_metric <- function(values,
                                          metric,
                                          weights = NULL,
                                          base = exp(1),
                                          resolution = NULL,
                                          classes_max = NULL) {
  metric <- match.arg(
    metric,
    c("shdi", "shei", "sidi", "siei", "msidi", "msiei", "pr",
      "prd", "rpr", "ta")
  )
  validate_scalar_numeric(base, "base", positive = TRUE)
  if (base == 1) {
    stop("`base` must not equal 1.", call. = FALSE)
  }

  totals <- .landscape_class_totals(values = values, weights = weights)
  if (length(totals) == 0) {
    return(NA_real_)
  }

  p <- totals / sum(totals)
  richness <- length(p)

  if (metric %in% c("prd", "ta")) {
    if (is.null(resolution)) {
      stop("`resolution` is required for PRD and TA.", call. = FALSE)
    }
    validate_scalar_numeric(resolution, "resolution", positive = TRUE)
    area_total <- sum(totals) * resolution^2 / 10000
  }

  shdi <- -sum(p * log(p, base = base))
  sidi <- 1 - sum(p^2)
  msidi <- -log(sum(p^2))

  switch(
    metric,
    shdi = shdi,
    shei = if (richness == 1) 0 else shdi / log(richness, base = base),
    sidi = sidi,
    siei = if (richness <= 1) NA_real_ else sidi / (1 - (1 / richness)),
    msidi = msidi,
    msiei = msidi / log(richness),
    pr = richness,
    prd = richness / area_total * 100,
    rpr = {
      if (is.null(classes_max)) {
        NA_real_
      } else {
        validate_scalar_numeric(classes_max, "classes_max", positive = TRUE)
        richness / classes_max * 100
      }
    },
    ta = area_total
  )
}


.landscape_shdi <- function(values, weights = NULL, base = exp(1)) {
  .landscape_composition_metric(
    values = values,
    metric = "shdi",
    weights = weights,
    base = base
  )
}


.landscape_composition_by_buffer <- function(d,
                                             r_stack.df,
                                             radius,
                                             metric,
                                             weights = NULL,
                                             base = exp(1),
                                             resolution = NULL,
                                             classes_max = NULL) {
  validate_scalar_numeric(radius, "radius", positive = TRUE)
  validate_scalar_numeric(base, "base", positive = TRUE)

  if (length(d) == 0) {
    stop("`d` must contain at least one distance.", call. = FALSE)
  }

  values <- as.matrix(r_stack.df)
  storage.mode(values) <- "double"
  if (nrow(values) != length(d)) {
    stop("`r_stack.df` must have one row for each distance in `d`.", call. = FALSE)
  }
  .landscape_validate_categorical_values(
    values,
    metric = metric,
    context = "Landscape composition metric"
  )

  if (!is.null(weights) && length(weights) != length(d)) {
    stop("`weights` must have the same length as `d`.", call. = FALSE)
  }

  out <- landscape_composition_metric_cpp(
    d = d,
    values = values,
    radius = radius,
    metric = metric,
    weights_ = weights,
    base = base,
    resolution_ = resolution,
    classes_max_ = classes_max
  )

  stats::setNames(out, colnames(values))
}


.landscape_shdi_by_buffer <- function(d,
                                      r_stack.df,
                                      radius,
                                      weights = NULL,
                                      base = exp(1)) {
  .landscape_composition_by_buffer(
    d = d,
    r_stack.df = r_stack.df,
    radius = radius,
    metric = "shdi",
    weights = weights,
    base = base
  )
}


.landscape_circle_kernel <- function(radius, resolution) {
  validate_scalar_numeric(radius, "radius", positive = TRUE)
  validate_scalar_numeric(resolution, "resolution", positive = TRUE)

  half_width <- ceiling(radius / resolution)
  offsets <- seq(-half_width, half_width)
  d <- sqrt(outer(offsets^2, offsets^2, "+")) * resolution

  matrix(as.numeric(d <= radius),
         nrow = length(offsets),
         ncol = length(offsets))
}


.landscape_cells_to_matrix <- function(values, cells, n_cols) {
  validate_scalar_numeric(n_cols, "n_cols", integerish = TRUE, positive = TRUE)

  if (length(values) != length(cells)) {
    stop("`values` and `cells` must have the same length.", call. = FALSE)
  }

  rows <- ((cells - 1) %/% n_cols) + 1
  cols <- ((cells - 1) %% n_cols) + 1
  row_min <- min(rows)
  col_min <- min(cols)

  out <- matrix(
    NA_real_,
    nrow = max(rows) - row_min + 1,
    ncol = max(cols) - col_min + 1
  )
  out[cbind(rows - row_min + 1, cols - col_min + 1)] <- values
  out
}


.landscape_min_perimeter <- function(n_cells) {
  if (n_cells <= 0) {
    return(NA_real_)
  }

  total_n <- trunc(sqrt(n_cells))
  total_m <- n_cells - total_n^2

  if (total_m == 0) {
    total_n * 4
  } else if (n_cells <= total_n * (1 + total_n)) {
    4 * total_n + 2
  } else {
    4 * total_n + 4
  }
}


.landscape_max_like_adjacencies <- function(n_cells) {
  if (n_cells <= 0) {
    return(NA_real_)
  }

  n <- trunc(sqrt(n_cells))
  m <- n_cells - n^2

  if (m == 0) {
    2 * n * (n - 1)
  } else if (m <= n) {
    2 * n * (n - 1) + 2 * m - 1
  } else {
    2 * n * (n - 1) + 2 * m - 2
  }
}


.landscape_edge_stats <- function(values, resolution) {
  validate_scalar_numeric(resolution, "resolution", positive = TRUE)

  if (!is.matrix(values)) {
    values <- as.matrix(values)
  }

  valid <- !is.na(values)
  n_valid <- sum(valid)
  if (n_valid == 0) {
    return(list(
      n_valid = 0,
      area_ha = NA_real_,
      internal_edge_count = NA_real_,
      boundary_edge_count = NA_real_,
      total_perimeter_count = NA_real_
    ))
  }

  right_diff <- right_valid <- matrix(FALSE, nrow = nrow(values), ncol = 0)
  if (ncol(values) > 1) {
    right_diff <- values[, -ncol(values), drop = FALSE] !=
      values[, -1, drop = FALSE]
    right_valid <- valid[, -ncol(valid), drop = FALSE] &
      valid[, -1, drop = FALSE]
  }

  down_diff <- down_valid <- matrix(FALSE, nrow = 0, ncol = ncol(values))
  if (nrow(values) > 1) {
    down_diff <- values[-nrow(values), , drop = FALSE] !=
      values[-1, , drop = FALSE]
    down_valid <- valid[-nrow(valid), , drop = FALSE] &
      valid[-1, , drop = FALSE]
  }

  internal_edge_count <- sum(right_diff & right_valid, na.rm = TRUE) +
    sum(down_diff & down_valid, na.rm = TRUE)

  padded_valid <- matrix(FALSE, nrow = nrow(valid) + 2, ncol = ncol(valid) + 2)
  padded_valid[seq_len(nrow(valid)) + 1, seq_len(ncol(valid)) + 1] <- valid
  row_idx <- seq_len(nrow(valid)) + 1
  col_idx <- seq_len(ncol(valid)) + 1
  boundary_edge_count <- sum(valid & !padded_valid[row_idx - 1, col_idx]) +
    sum(valid & !padded_valid[row_idx + 1, col_idx]) +
    sum(valid & !padded_valid[row_idx, col_idx - 1]) +
    sum(valid & !padded_valid[row_idx, col_idx + 1])

  list(
    n_valid = n_valid,
    area_ha = n_valid * resolution^2 / 10000,
    internal_edge_count = internal_edge_count,
    boundary_edge_count = boundary_edge_count,
    total_perimeter_count = internal_edge_count + boundary_edge_count
  )
}


.landscape_edge_metric <- function(values, resolution, metric) {
  metric <- match.arg(metric, c("ed", "te", "lsi"))
  stats <- .landscape_edge_stats(values = values, resolution = resolution)

  if (stats$n_valid == 0) {
    return(NA_real_)
  }

  switch(
    metric,
    ed = stats$internal_edge_count * resolution / stats$area_ha,
    te = stats$internal_edge_count * resolution,
    lsi = stats$total_perimeter_count / .landscape_min_perimeter(stats$n_valid)
  )
}


.landscape_edge_density <- function(values, resolution) {
  .landscape_edge_metric(values = values, resolution = resolution, metric = "ed")
}


.landscape_edge_by_buffer <- function(d,
                                      r_stack.df,
                                      cells,
                                      radius,
                                      resolution,
                                      n_cols,
                                      metric) {
  validate_scalar_numeric(radius, "radius", positive = TRUE)
  validate_scalar_numeric(resolution, "resolution", positive = TRUE)
  validate_scalar_numeric(n_cols, "n_cols", integerish = TRUE, positive = TRUE)

  if (length(d) == 0) {
    stop("`d` must contain at least one distance.", call. = FALSE)
  }
  if (length(cells) != length(d)) {
    stop("`cells` must have the same length as `d`.", call. = FALSE)
  }

  values <- as.matrix(r_stack.df)
  storage.mode(values) <- "double"
  if (nrow(values) != length(d)) {
    stop("`r_stack.df` must have one row for each distance in `d`.", call. = FALSE)
  }
  .landscape_validate_categorical_values(
    values,
    metric = metric,
    context = "Landscape edge metric"
  )

  out <- landscape_edge_metric_cpp(
    d = d,
    values = values,
    cells = as.numeric(cells),
    radius = radius,
    resolution = resolution,
    n_cols = n_cols,
    metric = metric
  )

  stats::setNames(out, colnames(values))
}


.landscape_ed_by_buffer <- function(d,
                                    r_stack.df,
                                    cells,
                                    radius,
                                    resolution,
                                    n_cols) {
  .landscape_edge_by_buffer(
    d = d,
    r_stack.df = r_stack.df,
    cells = cells,
    radius = radius,
    resolution = resolution,
    n_cols = n_cols,
    metric = "ed"
  )
}


.landscape_adjacency_table <- function(values) {
  if (!is.matrix(values)) {
    values <- as.matrix(values)
  }

  classes <- sort(unique(stats::na.omit(as.vector(values))))
  if (length(classes) == 0) {
    return(matrix(numeric(0), nrow = 0, ncol = 0))
  }

  out <- matrix(
    0,
    nrow = length(classes),
    ncol = length(classes),
    dimnames = list(as.character(classes), as.character(classes))
  )

  add_pairs <- function(a, b) {
    keep <- !is.na(a) & !is.na(b)
    a <- a[keep]
    b <- b[keep]
    if (length(a) == 0) {
      return(NULL)
    }

    ai <- match(a, classes)
    bi <- match(b, classes)
    for (i in seq_along(ai)) {
      out[ai[[i]], bi[[i]]] <<- out[ai[[i]], bi[[i]]] + 1
      out[bi[[i]], ai[[i]]] <<- out[bi[[i]], ai[[i]]] + 1
    }
    NULL
  }

  if (ncol(values) > 1) {
    add_pairs(
      as.vector(values[, -ncol(values), drop = FALSE]),
      as.vector(values[, -1, drop = FALSE])
    )
  }
  if (nrow(values) > 1) {
    add_pairs(
      as.vector(values[-nrow(values), , drop = FALSE]),
      as.vector(values[-1, , drop = FALSE])
    )
  }

  out
}


.landscape_adjacency_metric <- function(values, metric) {
  metric <- match.arg(metric, c("ai", "pladj", "contag"))

  tab <- .landscape_adjacency_table(values)
  if (length(tab) == 0) {
    return(NA_real_)
  }

  total <- sum(tab)
  if (metric == "ai") {
    class_counts <- .landscape_class_totals(as.vector(values))
    class_counts <- class_counts[match(rownames(tab), names(class_counts))]
    like_adjacencies <- diag(tab) / 2
    max_adjacencies <- vapply(
      class_counts,
      .landscape_max_like_adjacencies,
      numeric(1)
    )
    class_ai <- like_adjacencies / max_adjacencies * 100
    class_ai[!is.finite(class_ai)] <- NA_real_
    return(sum(class_ai * class_counts / sum(class_counts), na.rm = TRUE))
  }

  if (metric == "pladj") {
    if (total == 0) {
      return(0)
    }
    return(sum(diag(tab)) / total * 100)
  }

  if (nrow(tab) < 2 || total == 0) {
    return(NA_real_)
  }

  p <- tab / total
  esum <- sum(p * log(p), na.rm = TRUE)
  (1 + esum / (2 * log(nrow(tab)))) * 100
}


.landscape_adjacency_by_buffer <- function(d,
                                           r_stack.df,
                                           cells,
                                           radius,
                                           n_cols,
                                           metric) {
  validate_scalar_numeric(radius, "radius", positive = TRUE)
  validate_scalar_numeric(n_cols, "n_cols", integerish = TRUE, positive = TRUE)

  if (length(d) == 0) {
    stop("`d` must contain at least one distance.", call. = FALSE)
  }
  if (length(cells) != length(d)) {
    stop("`cells` must have the same length as `d`.", call. = FALSE)
  }

  values <- as.matrix(r_stack.df)
  storage.mode(values) <- "double"
  if (nrow(values) != length(d)) {
    stop("`r_stack.df` must have one row for each distance in `d`.", call. = FALSE)
  }
  .landscape_validate_categorical_values(
    values,
    metric = metric,
    max_classes = .landscape_adjacency_class_ceiling(metric),
    context = "Landscape adjacency metric"
  )

  out <- landscape_adjacency_metric_cpp(
    d = d,
    values = values,
    cells = as.numeric(cells),
    radius = radius,
    n_cols = n_cols,
    metric = metric
  )

  stats::setNames(out, colnames(values))
}


.landscape_validate_single_layer_raster <- function(raster, metric_label) {
  if (!inherits(raster, "SpatRaster")) {
    stop("`raster` must be a terra `SpatRaster`.", call. = FALSE)
  }
  if (terra::nlyr(raster) != 1) {
    stop("`raster` must contain exactly one categorical layer.", call. = FALSE)
  }

  raster_res <- terra::res(raster)
  if (!isTRUE(all.equal(raster_res[[1]], raster_res[[2]]))) {
    stop(sprintf("%s FFT projection currently requires square raster cells.",
                 metric_label),
         call. = FALSE)
  }

  raster_res[[1]]
}


.landscape_class_count_rasters <- function(raster, radius, na.rm = TRUE) {
  validate_scalar_numeric(radius, "radius", positive = TRUE)
  validate_scalar_logical(na.rm, "na.rm")

  resolution <- .landscape_validate_single_layer_raster(raster, "Composition")
  values <- terra::as.matrix(raster, wide = TRUE)
  .landscape_validate_categorical_values(
    values,
    metric = "composition",
    context = "Landscape composition projection"
  )
  classes <- sort(unique(stats::na.omit(as.vector(values))))
  kernel <- .landscape_circle_kernel(radius = radius, resolution = resolution)

  class_counts <- vector("list", length(classes))
  for (i in seq_along(classes)) {
    class_indicator <- matrix(
      as.numeric(values == classes[[i]]),
      nrow = nrow(values),
      ncol = ncol(values)
    )
    class_indicator[is.na(class_indicator)] <- 0
    class_counts[[i]] <- .landscape_clean_count_matrix(fft_convolution(
      x = class_indicator,
      kernel = kernel,
      fun = "sum",
      na.rm = na.rm
    ))
  }

  total <- if (length(class_counts) == 0) {
    matrix(0, nrow = nrow(values), ncol = ncol(values))
  } else {
    .landscape_clean_count_matrix(Reduce("+", class_counts))
  }

  list(
    values = values,
    classes = classes,
    class_counts = class_counts,
    total = total,
    resolution = resolution
  )
}


.landscape_composition_raster_fft <- function(raster,
                                              radius,
                                              metric,
                                              base = exp(1),
                                              classes_max = NULL,
                                              na.rm = TRUE) {
  metric <- match.arg(
    metric,
    c("shdi", "shei", "sidi", "siei", "msidi", "msiei", "pr",
      "prd", "rpr", "ta")
  )
  validate_scalar_numeric(base, "base", positive = TRUE)
  if (base == 1) {
    stop("`base` must not equal 1.", call. = FALSE)
  }
  if (!is.null(classes_max)) {
    validate_scalar_numeric(classes_max, "classes_max", positive = TRUE)
  }

  counts <- .landscape_class_count_rasters(
    raster = raster,
    radius = radius,
    na.rm = na.rm
  )
  out <- terra::rast(raster)

  if (length(counts$classes) == 0) {
    terra::values(out) <- NA_real_
    names(out) <- paste0(names(raster), "_", metric)
    return(out)
  }

  total <- counts$total
  present <- matrix(0, nrow = nrow(total), ncol = ncol(total))
  shdi <- matrix(0, nrow = nrow(total), ncol = ncol(total))
  sum_sq <- matrix(0, nrow = nrow(total), ncol = ncol(total))

  for (count in counts$class_counts) {
    p <- count / total
    active <- is.finite(p) & count > 0 & total > 0
    present <- present + active
    shdi_term <- matrix(0, nrow = nrow(total), ncol = ncol(total))
    shdi_term[active] <- p[active] * log(p[active], base = base)
    shdi <- shdi - shdi_term
    sum_sq[active] <- sum_sq[active] + p[active]^2
  }

  sidi <- 1 - sum_sq
  msidi <- -log(sum_sq)
  area_total <- total * counts$resolution^2 / 10000

  values <- switch(
    metric,
    shdi = shdi,
    shei = ifelse(present <= 1, 0, shdi / log(present, base = base)),
    sidi = sidi,
    siei = ifelse(present <= 1, NA_real_, sidi / (1 - (1 / present))),
    msidi = msidi,
    msiei = msidi / log(present),
    pr = present,
    prd = present / area_total * 100,
    rpr = {
      if (is.null(classes_max)) {
        matrix(NA_real_, nrow = nrow(total), ncol = ncol(total))
      } else {
        present / classes_max * 100
      }
    },
    ta = area_total
  )

  values[!is.finite(total) | total <= 0] <- NA_real_
  terra::values(out) <- as.vector(values)
  names(out) <- paste0(names(raster), "_", metric)

  out
}


.landscape_shdi_raster_fft <- function(raster,
                                       radius,
                                       base = exp(1),
                                       na.rm = TRUE) {
  .landscape_composition_raster_fft(
    raster = raster,
    radius = radius,
    metric = "shdi",
    base = base,
    na.rm = na.rm
  )
}


.landscape_edge_kernel <- function(radius, resolution, orientation) {
  validate_scalar_numeric(radius, "radius", positive = TRUE)
  validate_scalar_numeric(resolution, "resolution", positive = TRUE)
  orientation <- match.arg(orientation, c("right", "down"))

  half_width <- ceiling(radius / resolution)
  offsets <- seq(-half_width, half_width)
  row_offset <- matrix(rep(offsets, times = length(offsets)),
                       nrow = length(offsets))
  col_offset <- matrix(rep(offsets, each = length(offsets)),
                       nrow = length(offsets))

  if (identical(orientation, "right")) {
    neighbour_row <- row_offset
    neighbour_col <- col_offset - 1
  } else {
    neighbour_row <- row_offset - 1
    neighbour_col <- col_offset
  }

  cell_in <- sqrt(row_offset^2 + col_offset^2) * resolution <= radius
  neighbour_in <- sqrt(neighbour_row^2 + neighbour_col^2) * resolution <= radius

  matrix(as.numeric(cell_in & neighbour_in),
         nrow = length(offsets),
         ncol = length(offsets))
}


.landscape_edge_counts_raster_fft <- function(raster, radius, na.rm = TRUE) {
  validate_scalar_numeric(radius, "radius", positive = TRUE)
  validate_scalar_logical(na.rm, "na.rm")
  resolution <- .landscape_validate_single_layer_raster(raster, "Edge")

  values <- terra::as.matrix(raster, wide = TRUE)
  .landscape_validate_categorical_values(
    values,
    metric = "edge",
    context = "Landscape edge projection"
  )
  valid <- !is.na(values)
  right_edges <- matrix(0, nrow = nrow(values), ncol = ncol(values))
  down_edges <- matrix(0, nrow = nrow(values), ncol = ncol(values))
  right_valid_pairs <- matrix(0, nrow = nrow(values), ncol = ncol(values))
  down_valid_pairs <- matrix(0, nrow = nrow(values), ncol = ncol(values))

  if (ncol(values) > 1) {
    right_valid <- valid[, -ncol(valid), drop = FALSE] &
      valid[, -1, drop = FALSE]
    right_valid_pairs[, -ncol(values)] <- as.numeric(right_valid)
    right_edges[, -ncol(values)] <- as.numeric(
      values[, -ncol(values), drop = FALSE] != values[, -1, drop = FALSE] &
        right_valid
    )
  }
  if (nrow(values) > 1) {
    down_valid <- valid[-nrow(valid), , drop = FALSE] &
      valid[-1, , drop = FALSE]
    down_valid_pairs[-nrow(values), ] <- as.numeric(down_valid)
    down_edges[-nrow(values), ] <- as.numeric(
      values[-nrow(values), , drop = FALSE] != values[-1, , drop = FALSE] &
        down_valid
    )
  }

  area_kernel <- .landscape_circle_kernel(radius = radius,
                                          resolution = resolution)
  right_kernel <- .landscape_edge_kernel(radius = radius,
                                         resolution = resolution,
                                         orientation = "right")
  down_kernel <- .landscape_edge_kernel(radius = radius,
                                        resolution = resolution,
                                        orientation = "down")

  area_cells <- .landscape_clean_count_matrix(fft_convolution(valid * 1,
                                                              area_kernel,
                                                              fun = "sum",
                                                              na.rm = na.rm))
  right_count <- .landscape_clean_count_matrix(fft_convolution(right_edges,
                                                               right_kernel,
                                                               fun = "sum",
                                                               na.rm = na.rm))
  down_count <- .landscape_clean_count_matrix(fft_convolution(down_edges,
                                                              down_kernel,
                                                              fun = "sum",
                                                              na.rm = na.rm))
  right_valid_count <- .landscape_clean_count_matrix(fft_convolution(right_valid_pairs,
                                                                     right_kernel,
                                                                     fun = "sum",
                                                                     na.rm = na.rm))
  down_valid_count <- .landscape_clean_count_matrix(fft_convolution(down_valid_pairs,
                                                                    down_kernel,
                                                                    fun = "sum",
                                                                    na.rm = na.rm))

  list(
    values = values,
    area_cells = area_cells,
    internal_edge_count = .landscape_clean_count_matrix(right_count + down_count),
    valid_pair_count = .landscape_clean_count_matrix(right_valid_count + down_valid_count),
    resolution = resolution,
    right_kernel = right_kernel,
    down_kernel = down_kernel
  )
}


.landscape_min_perimeter_matrix <- function(n_cells) {
  n_cells <- round(n_cells)
  out <- matrix(NA_real_, nrow = nrow(n_cells), ncol = ncol(n_cells))
  active <- is.finite(n_cells) & n_cells > 0

  total_n <- trunc(sqrt(n_cells[active]))
  total_m <- n_cells[active] - total_n^2

  out[active] <- ifelse(
    total_m == 0,
    total_n * 4,
    ifelse(n_cells[active] <= total_n * (1 + total_n),
           4 * total_n + 2,
           4 * total_n + 4)
  )

  out
}


.landscape_max_like_adjacencies_matrix <- function(n_cells) {
  n_cells <- round(n_cells)
  out <- matrix(NA_real_, nrow = nrow(n_cells), ncol = ncol(n_cells))
  active <- is.finite(n_cells) & n_cells > 0

  n <- trunc(sqrt(n_cells[active]))
  m <- n_cells[active] - n^2

  out[active] <- ifelse(
    m == 0,
    2 * n * (n - 1),
    ifelse(m <= n,
           2 * n * (n - 1) + 2 * m - 1,
           2 * n * (n - 1) + 2 * m - 2)
  )

  out
}


.landscape_edge_raster_fft <- function(raster,
                                       radius,
                                       metric,
                                       na.rm = TRUE) {
  metric <- match.arg(metric, c("ed", "te", "lsi"))
  counts <- .landscape_edge_counts_raster_fft(
    raster = raster,
    radius = radius,
    na.rm = na.rm
  )

  edge_total <- counts$internal_edge_count * counts$resolution
  area_total <- counts$area_cells * counts$resolution^2 / 10000
  values <- switch(
    metric,
    ed = edge_total / area_total,
    te = edge_total,
    lsi = {
      boundary_edge_count <- 4 * counts$area_cells - 2 * counts$valid_pair_count
      total_perimeter_count <- counts$internal_edge_count + boundary_edge_count
      total_perimeter_count / .landscape_min_perimeter_matrix(counts$area_cells)
    }
  )
  values[!is.finite(values) | counts$area_cells <= 0] <- NA_real_

  out <- terra::rast(raster)
  terra::values(out) <- as.vector(values)
  names(out) <- paste0(names(raster), "_", metric)
  out
}


.landscape_ed_raster_fft <- function(raster,
                                     radius,
                                     na.rm = TRUE) {
  .landscape_edge_raster_fft(
    raster = raster,
    radius = radius,
    metric = "ed",
    na.rm = na.rm
  )
}


.landscape_te_raster_fft <- function(raster,
                                     radius,
                                     na.rm = TRUE) {
  .landscape_edge_raster_fft(
    raster = raster,
    radius = radius,
    metric = "te",
    na.rm = na.rm
  )
}


.landscape_pair_count_raster_fft <- function(values,
                                             from,
                                             to,
                                             right_kernel,
                                             down_kernel,
                                             na.rm) {
  right_pairs <- matrix(0, nrow = nrow(values), ncol = ncol(values))
  down_pairs <- matrix(0, nrow = nrow(values), ncol = ncol(values))

  if (ncol(values) > 1) {
    right_pairs[, -ncol(values)] <- as.numeric(
      values[, -ncol(values), drop = FALSE] == from &
        values[, -1, drop = FALSE] == to
    )
    right_pairs[is.na(right_pairs)] <- 0
  }
  if (nrow(values) > 1) {
    down_pairs[-nrow(values), ] <- as.numeric(
      values[-nrow(values), , drop = FALSE] == from &
        values[-1, , drop = FALSE] == to
    )
    down_pairs[is.na(down_pairs)] <- 0
  }

  .landscape_clean_count_matrix(
    fft_convolution(right_pairs,
                    right_kernel,
                    fun = "sum",
                    na.rm = na.rm) +
      fft_convolution(down_pairs,
                      down_kernel,
                      fun = "sum",
                      na.rm = na.rm)
  )
}


.landscape_adjacency_raster_fft <- function(raster,
                                            radius,
                                            metric,
                                            na.rm = TRUE) {
  metric <- match.arg(metric, c("ai", "pladj", "contag"))
  projection_values <- terra::as.matrix(raster, wide = TRUE)
  .landscape_validate_categorical_values(
    projection_values,
    metric = metric,
    max_classes = .landscape_adjacency_class_ceiling(metric),
    context = "Landscape adjacency projection"
  )

  counts <- .landscape_edge_counts_raster_fft(
    raster = raster,
    radius = radius,
    na.rm = na.rm
  )
  values <- counts$values
  .landscape_validate_categorical_values(
    values,
    metric = metric,
    max_classes = .landscape_adjacency_class_ceiling(metric),
    context = "Landscape adjacency projection"
  )
  classes <- sort(unique(stats::na.omit(as.vector(values))))

  out <- terra::rast(raster)
  if (length(classes) == 0) {
    terra::values(out) <- NA_real_
    names(out) <- paste0(names(raster), "_", metric)
    return(out)
  }

  total_source <- counts$valid_pair_count

  if (metric == "pladj") {
    like_source <- matrix(0, nrow = nrow(values), ncol = ncol(values))
    for (class in classes) {
      like_source <- like_source + .landscape_pair_count_raster_fft(
        values = values,
        from = class,
        to = class,
        right_kernel = counts$right_kernel,
        down_kernel = counts$down_kernel,
        na.rm = na.rm
      )
    }
    result <- like_source / total_source * 100
    result[total_source <= 0 & counts$area_cells > 0] <- 0
    result[counts$area_cells <= 0] <- NA_real_
  } else if (metric == "ai") {
    class_counts <- .landscape_class_count_rasters(
      raster = raster,
      radius = radius,
      na.rm = na.rm
    )
    result <- matrix(0, nrow = nrow(values), ncol = ncol(values))

    for (i in seq_along(class_counts$classes)) {
      class <- class_counts$classes[[i]]
      count <- class_counts$class_counts[[i]]
      like_source <- .landscape_pair_count_raster_fft(
        values = values,
        from = class,
        to = class,
        right_kernel = counts$right_kernel,
        down_kernel = counts$down_kernel,
        na.rm = na.rm
      )
      max_adjacencies <- .landscape_max_like_adjacencies_matrix(count)
      class_ai <- like_source / max_adjacencies * 100
      contribution <- class_ai * count / class_counts$total
      contribution[!is.finite(contribution)] <- 0
      result <- result + contribution
    }

    result[class_counts$total <= 0 | counts$area_cells <= 0] <- NA_real_
  } else {
    class_counts <- .landscape_class_count_rasters(
      raster = raster,
      radius = radius,
      na.rm = na.rm
    )
    present <- matrix(0, nrow = nrow(values), ncol = ncol(values))
    for (count in class_counts$class_counts) {
      present <- present + (count > 0)
    }

    entropy <- matrix(0, nrow = nrow(values), ncol = ncol(values))
    total_ordered <- total_source * 2
    for (from_idx in seq_along(classes)) {
      for (to_idx in seq(from_idx, length(classes))) {
        from <- classes[[from_idx]]
        to <- classes[[to_idx]]
        source_ab <- .landscape_pair_count_raster_fft(
          values = values,
          from = from,
          to = to,
          right_kernel = counts$right_kernel,
          down_kernel = counts$down_kernel,
          na.rm = na.rm
        )
        if (from_idx == to_idx) {
          pair_count <- source_ab * 2
          multiplier <- 1
        } else {
          source_ba <- .landscape_pair_count_raster_fft(
            values = values,
            from = to,
            to = from,
            right_kernel = counts$right_kernel,
            down_kernel = counts$down_kernel,
            na.rm = na.rm
          )
          pair_count <- source_ab + source_ba
          multiplier <- 2
        }
        p <- pair_count / total_ordered
        active <- is.finite(p) & p > 0
        entropy_term <- matrix(0, nrow = nrow(values), ncol = ncol(values))
        entropy_term[active] <- multiplier * p[active] * log(p[active])
        entropy <- entropy + entropy_term
      }
    }

    result <- (1 + entropy / (2 * log(present))) * 100
    result[present < 2 | total_source <= 0 | counts$area_cells <= 0] <- NA_real_
  }

  result[!is.finite(result)] <- NA_real_
  terra::values(out) <- as.vector(result)
  names(out) <- paste0(names(raster), "_", metric)
  out
}
