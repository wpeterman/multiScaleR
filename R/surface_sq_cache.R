# Exact distance thresholds, with stable centered online sample variance.
# Preparation is local to an optimization context, never a cross-fit cache.
.msr_prepare_surface_sq <- function(d_list, cov_df, specs,
                                    budget = 32 * 1024^2) {
  if (is.null(specs) || is.null(d_list) || is.null(cov_df) ||
      length(d_list) != length(cov_df)) return(NULL)
  specs <- specs[specs$type == "surface" & specs$metric == "sq" &
                   !specs$weighted & specs$optimize, , drop = FALSE]
  if (!nrow(specs)) return(NULL)
  out <- vector("list", length(d_list))
  for (i in seq_along(d_list)) {
    d <- d_list[[i]]
    for (source in unique(specs$source)) {
      z <- cov_df[[i]][, source]
      if (!is.numeric(z) || length(z) != length(d)) next
      keep <- is.finite(d) & is.finite(z)
      n <- sum(keep)
      # Bound both persistent storage and the temporary per-source work.
      if (n < 2L || 16 * n > min(budget, 4 * 1024^2)) next
      ord <- order(d[keep])
      distance <- d[keep][ord]
      value <- as.double(z[keep][ord])
      # Shifting first prevents large absolute elevations from degrading the
      # online mean update. Never subtract uncentered sums of squares.
      variance <- msr_surface_sq_prefix_cpp(value)
      if (any(!is.finite(variance[-1L]) | variance[-1L] < 0)) next
      # Include all cells tied at the boundary, matching d <= radius.
      ends <- c(which(diff(distance) != 0), n)
      prepared <- list(metric = "sq", threshold = distance[ends],
                       value = sqrt(variance[ends]), min_distance = min(d[is.finite(d)]))
      cost <- as.numeric(utils::object.size(prepared))
      if (cost > budget) next
      budget <- budget - cost
      for (name in specs$covariate[specs$source == source])
        out[[i]][[name]] <- prepared
    }
  }
  out
}

.msr_prepare_surface_neighbors <- function(d_list, cov_df, specs, resolution, n_cols) {
  neighbors <- .msr_prepare_surface_neighbors_base(d_list, cov_df, specs, resolution, n_cols)
  # Prefix roundoff can alter a nonsmooth radius search trajectory. Keep the
  # established scalar path until end-to-end equivalence has been validated.
  if (!isTRUE(getOption("multiScaleR.experimental_sq_cache", FALSE))) return(neighbors)
  sq <- .msr_prepare_surface_sq(d_list, cov_df, specs)
  if (is.null(neighbors)) return(sq)
  if (is.null(sq)) return(neighbors)
  Map(c, neighbors, sq)
}
