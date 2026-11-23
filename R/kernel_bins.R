# Distance-binned acceleration for kernel-weighted covariates --------------
#
# The optimizer in `multiScale_optim()` repeatedly recomputes, for every point,
# a normalized kernel-weighted mean over all raster cells inside the `max_D`
# buffer. The cell distances are FIXED across all optimizer evaluations; only
# the kernel scale (`sigma`, and `shape` for the exponential power kernel)
# changes. Because the kernel weight is a pure function of distance, the
# weighted mean can be computed from per-point distance-binned summaries that
# are precomputed once in `kernel_prep()`:
#
#   cov.w[i] = sum_b w(d_rep[b]; sigma) * vsum[i, b]
#            / sum_b w(d_rep[b]; sigma) * csum[i, b]
#
# where, for point i and distance bin b,
#   vsum[i, b] = sum of FINITE cell values whose distance falls in bin b,
#   csum[i, b] = count of FINITE cells whose distance falls in bin b,
#   d_rep[b]   = representative (pooled count-weighted mean) distance of bin b.
#
# This turns each per-evaluation weighting from O(n_cells) per point into a
# single O(n_bins) matrix-vector product, and lets a workflow store O(n_bins)
# summaries per point instead of every buffer cell. Normalizing by the
# finite-cell weight sum reproduces the missing-data renormalization used by the
# exact C++ path (`scale_type_sparse`). Only kernel-type covariates are binned;
# landscape-metric covariates require the full per-cell spatial detail and keep
# the exact path.


# Sum `values` into `nbins` bins indexed by integer `bins` (NA bins dropped).
.msr_bin_accumulate <- function(values, bins, nbins) {
  out <- numeric(nbins)
  keep <- !is.na(bins)
  if (!any(keep)) {
    return(out)
  }
  grouped <- rowsum(values[keep], bins[keep], reorder = FALSE)
  out[as.integer(rownames(grouped))] <- grouped[, 1]
  out
}


# Build per-point distance-binned summaries for kernel-type covariates.
#
# d_list      : list (one per point) of scaled distance vectors (as stored in
#               `kernel_inputs$d_list`, i.e. already divided by `unit_conv`).
# value_list  : list (one per point) of objects supporting `[, source]` column
#               extraction of raw raster values (the sparse matrices stored in
#               `kernel_inputs$raw_cov`).
# covariates  : character vector of kernel-covariate names to bin.
# sources     : character vector (same length) of source raster-layer names.
# nbins       : number of equal-width distance bins.
#
# Returns NULL when there is nothing to bin, otherwise a list with `bin_rep`
# (length nbins), `vsum`/`csum` (named lists of n_points x nbins matrices), and
# bookkeeping fields.
.msr_build_kernel_bins <- function(d_list,
                                   value_list,
                                   covariates,
                                   sources,
                                   nbins = 256L,
                                   point_ids = NULL) {
  n_pts <- length(d_list)
  if (n_pts == 0 || length(covariates) == 0) {
    return(NULL)
  }
  nbins <- as.integer(nbins)
  if (length(nbins) != 1 || is.na(nbins) || nbins < 2) {
    stop("`nbins` must be a single integer >= 2.", call. = FALSE)
  }

  max_d <- max(vapply(d_list,
                      function(x) if (length(x)) max(x, na.rm = TRUE) else 0,
                      numeric(1)))
  if (!is.finite(max_d) || max_d <= 0) {
    return(NULL)
  }

  # Equal-width breaks over [0, max_d]; nudge the top edge so the farthest cell
  # is included by `.bincode(..., include.lowest = TRUE)`.
  breaks <- seq(0, max_d * (1 + 1e-9), length.out = nbins + 1L)
  mids <- (breaks[-1L] + breaks[-length(breaks)]) / 2

  # Cache each point's bin codes once; accumulate the pooled representative
  # distance per bin (count-weighted mean of the actual cell distances).
  bcode_list <- vector("list", n_pts)
  d_sum <- numeric(nbins)
  d_cnt <- numeric(nbins)
  for (i in seq_len(n_pts)) {
    di <- d_list[[i]]
    b <- .bincode(di, breaks, include.lowest = TRUE)
    bcode_list[[i]] <- b
    d_sum <- d_sum + .msr_bin_accumulate(di, b, nbins)
    keep <- !is.na(b)
    if (any(keep)) {
      d_cnt <- d_cnt + tabulate(b[keep], nbins)
    }
  }
  bin_rep <- d_sum / d_cnt
  empty <- !is.finite(bin_rep)
  bin_rep[empty] <- mids[empty]

  vsum <- stats::setNames(vector("list", length(covariates)), covariates)
  csum <- stats::setNames(vector("list", length(covariates)), covariates)

  for (k in seq_along(covariates)) {
    cov_name <- covariates[[k]]
    src <- sources[[k]]
    V <- matrix(0, nrow = n_pts, ncol = nbins)
    C <- matrix(0, nrow = n_pts, ncol = nbins)
    for (i in seq_len(n_pts)) {
      vals <- as.numeric(value_list[[i]][, src])
      fin <- is.finite(vals)
      if (!any(fin)) {
        next
      }
      b <- bcode_list[[i]]
      V[i, ] <- .msr_bin_accumulate(vals[fin], b[fin], nbins)
      C[i, ] <- tabulate(b[fin][!is.na(b[fin])], nbins)
    }
    vsum[[cov_name]] <- V
    csum[[cov_name]] <- C
  }

  if (is.null(point_ids)) {
    point_ids <- names(d_list)
  }

  list(
    bin_rep = bin_rep,
    vsum = vsum,
    csum = csum,
    covariates = covariates,
    sources = stats::setNames(sources, covariates),
    nbins = nbins,
    n_points = n_pts,
    point_ids = point_ids
  )
}


# Distance-only kernel weight (the normalization constant cancels in the
# numerator/denominator ratio, so only the shape matters). Mirrors the kernels
# implemented in `src/scale_type_sparse.cpp`.
.msr_kernel_shape_weight <- function(d, sigma, shape, kernel) {
  switch(
    kernel,
    gaussian = exp(-(d * d) / (2 * sigma * sigma)),
    exp      = exp(-d / sigma),
    fixed    = as.numeric(d < sigma),
    expow    = exp(-(d / sigma)^shape),
    stop("Invalid kernel type.", call. = FALSE)
  )
}


# Compute binned kernel-weighted means for `covs`, vectorized across points.
#
# binned          : structure from `.msr_build_kernel_bins()`.
# covs            : kernel covariates to evaluate (subset of binned$covariates).
# param_covariates: optimized-covariate ordering used to index `sigma`/`shape`
#                   (matches `.msr_eval_scale_vars()`).
# sigma, shape    : current parameter vectors (scaled units).
# kernel          : kernel name.
# eval_idx        : optional integer vector selecting/ordering the point rows
#                   (the `complete_idx` subset used by `kernel_scale_fn`).
#
# Returns an n_eval x length(covs) matrix (column names = covs).
.msr_binned_kernel_means <- function(binned,
                                     covs,
                                     param_covariates,
                                     sigma,
                                     shape,
                                     kernel,
                                     eval_idx = NULL) {
  d_rep <- binned$bin_rep
  n_eval <- if (is.null(eval_idx)) binned$n_points else length(eval_idx)
  out <- matrix(NA_real_, nrow = n_eval, ncol = length(covs),
                dimnames = list(NULL, covs))

  for (cov_name in covs) {
    param_idx <- match(cov_name, param_covariates)
    sig <- sigma[[param_idx]]
    shp <- if (identical(kernel, "expow")) shape[[param_idx]] else NULL
    w <- .msr_kernel_shape_weight(d_rep, sig, shp, kernel)

    V <- binned$vsum[[cov_name]]
    C <- binned$csum[[cov_name]]
    # Complete-case indices frequently select every row. Reuse these read-only
    # matrices instead of copying both for every optimizer/model refit call.
    if (!is.null(eval_idx) && !identical(eval_idx, seq_len(nrow(V)))) {
      V <- V[eval_idx, , drop = FALSE]
      C <- C[eval_idx, , drop = FALSE]
    }

    num <- as.numeric(V %*% w)
    den <- as.numeric(C %*% w)
    res <- num / den
    res[!is.finite(den) | den <= 0] <- NA_real_
    out[, cov_name] <- res
  }

  out
}


# Restrict a binned structure to a subset of covariates (used to keep only the
# covariates that appear in the fitted model).
.msr_subset_binned <- function(binned, covariates) {
  if (is.null(binned)) {
    return(NULL)
  }
  keep <- intersect(binned$covariates, covariates)
  if (length(keep) == 0) {
    return(NULL)
  }
  binned$covariates <- keep
  binned$vsum <- binned$vsum[keep]
  binned$csum <- binned$csum[keep]
  binned$sources <- binned$sources[keep]
  binned
}


# Distance-binned acceleration for kernel-weighted surface moments ----------
#
# The weighted surface metrics sq, ssk and sku (`surface_var(weighted = TRUE)`)
# are functions of kernel-weighted power sums of each buffer's cell values,
#
#   A_k[i] = sum_c w(d_ic; sigma) * (v_ic - vbar_i)^k,   k = 0, ..., K,
#
# where vbar_i is the unweighted mean of the finite cells around point i.
# Centering on vbar_i keeps the conversion from power sums to central moments
# stable. As for kernel means, only the weights change with sigma, so the power
# sums can be summarized once per point in `kernel_prep()`.
#
# A single representative distance per bin is too coarse for second and higher
# moments: it leaves weighted SQ about ten times less accurate than a kernel
# mean at the same bin count. Instead each cell is split between the two
# nearest nodes of an equal-width distance grid, in proportion to its position
# between them (linear binning). Evaluating the kernel at the nodes then
# reproduces the exact power sums with an error of second order in the node
# spacing. At the default 256 bins this is comparable to, or smaller than, the
# kernel-mean binning error. `sa` (mean absolute deviation) puts the
# sigma-dependent weighted mean inside an absolute value and has no such
# decomposition, and the step weights of `kernel = "fixed"` do not interpolate,
# so both keep the exact per-cell path.

.msr_surface_bin_metrics <- function() c("sq", "ssk", "sku")

# Weighted surface specifications eligible for linear binning.
.msr_surface_bin_specs <- function(scale_vars, kernel) {
  if (is.null(scale_vars) || !nrow(scale_vars) || identical(kernel, "fixed") ||
      !all(c("type", "metric", "weighted") %in% names(scale_vars)))
    return(NULL)
  keep <- scale_vars$type == "surface" & scale_vars$weighted %in% TRUE &
    scale_vars$metric %in% .msr_surface_bin_metrics()
  if ("optimize" %in% names(scale_vars)) keep <- keep & scale_vars$optimize %in% TRUE
  if (!any(keep)) return(NULL)
  scale_vars[keep, , drop = FALSE]
}

# Build per-point, per-node centered power sums for weighted surface metrics.
#
# d_list, value_list: per-point distances and source-cell values, as stored in
#   `kernel_inputs$d_list` and `kernel_inputs$raw_cov`.
# specs: rows of `scale_vars` from `.msr_surface_bin_specs()`.
# nodes: optional node grid from a previous build, so that cross-validation
#   folds reuse the prepared geometry. Otherwise `nbins` equal intervals span
#   0 to the largest prepared distance.
#
# Returns NULL when there is nothing to bin, otherwise a list with `nodes`,
# `nbins`, the covariate-to-source and covariate-to-metric maps, `sums` (one
# list of n_points x (nbins + 1) matrices per source, for powers 0 to K), and
# `point_ids`.
.msr_build_surface_bins <- function(d_list, value_list, specs, nbins = 256L,
                                    point_ids = NULL, nodes = NULL) {
  n_pts <- length(d_list)
  if (n_pts == 0 || is.null(specs) || !nrow(specs)) return(NULL)
  if (is.null(nodes)) {
    nbins <- as.integer(nbins)
    if (length(nbins) != 1 || is.na(nbins) || nbins < 2)
      stop("`nbins` must be a single integer >= 2.", call. = FALSE)
    max_d <- max(vapply(d_list, function(x)
      if (length(x)) max(x, na.rm = TRUE) else 0, numeric(1)))
    if (!is.finite(max_d) || max_d <= 0) return(NULL)
    nodes <- seq(0, max_d * (1 + 1e-9), length.out = nbins + 1L)
  }
  nbins <- length(nodes) - 1L
  h <- nodes[[2L]]
  sources <- unique(specs$source)
  # Powers up to 2 suffice for SQ; skewness and kurtosis need 4.
  max_power <- vapply(sources, function(s)
    if (any(specs$metric[specs$source == s] %in% c("ssk", "sku"))) 4L else 2L,
    integer(1))
  names(max_power) <- sources
  sums <- stats::setNames(vector("list", length(sources)), sources)
  for (s in sources) {
    n_pow <- max_power[[s]] + 1L
    mats <- replicate(n_pow, matrix(0, n_pts, nbins + 1L), simplify = FALSE)
    for (i in seq_len(n_pts)) {
      d <- d_list[[i]]
      v <- as.numeric(value_list[[i]][, s])
      keep <- is.finite(v) & is.finite(d)
      if (!any(keep)) next
      x <- d[keep] / h
      if (any(x > nbins + 1e-6))
        stop("Prepared distances exceed the surface bin grid.", call. = FALSE)
      lo <- pmin(floor(x), nbins - 1)
      frac <- x - lo
      v <- v[keep] - mean(v[keep])
      powers <- outer(v, seq_len(n_pow) - 1L, "^")
      # Split every cell's powers between its two neighboring nodes.
      lower <- rowsum(powers * (1 - frac), lo + 1L, reorder = FALSE)
      upper <- rowsum(powers * frac, lo + 2L, reorder = FALSE)
      node <- matrix(0, nbins + 1L, n_pow)
      node[as.integer(rownames(lower)), ] <- lower
      ui <- as.integer(rownames(upper))
      node[ui, ] <- node[ui, , drop = FALSE] + upper
      for (k in seq_len(n_pow)) mats[[k]][i, ] <- node[, k]
    }
    sums[[s]] <- mats
  }
  if (is.null(point_ids)) point_ids <- names(d_list)
  list(nodes = nodes, nbins = nbins,
       covariates = specs$covariate,
       sources = stats::setNames(specs$source, specs$covariate),
       metrics = stats::setNames(specs$metric, specs$covariate),
       sums = sums, n_points = n_pts, point_ids = point_ids)
}

# Keep only the binned surface covariates used by a model.
.msr_subset_surface_bins <- function(surface_binned, covariates) {
  if (is.null(surface_binned)) return(NULL)
  keep <- intersect(surface_binned$covariates, covariates)
  if (!length(keep)) return(NULL)
  surface_binned$covariates <- keep
  surface_binned$sources <- surface_binned$sources[keep]
  surface_binned$metrics <- surface_binned$metrics[keep]
  surface_binned$sums <- surface_binned$sums[unique(surface_binned$sources)]
  surface_binned
}

# Rows of the binned sums for the evaluated points. Points are matched by ID so
# that any subset or reordering of the cell lists stays aligned. NULL means the
# rows cannot be matched and the caller must use the exact path.
.msr_surface_bin_rows <- function(surface_binned, ids, n) {
  if (is.null(surface_binned) || !n) return(NULL)
  if (is.null(ids))
    return(if (n == surface_binned$n_points) seq_len(n) else NULL)
  if (is.null(surface_binned$point_ids) || length(ids) != n) return(NULL)
  rows <- match(ids, surface_binned$point_ids)
  if (anyNA(rows)) NULL else rows
}

# Weighted surface metrics from binned power sums, one column per covariate.
# Mirrors `.surface_weighted_metric()`: SQ is zero for zero spread, skewness
# and kurtosis are missing there, and a point without positive weight is NA.
.msr_binned_surface_metrics <- function(surface_binned, covs, param_covariates,
                                        sigma, shape, kernel, rows) {
  out <- matrix(NA_real_, nrow = length(rows), ncol = length(covs),
                dimnames = list(NULL, covs))
  all_rows <- identical(rows, seq_len(surface_binned$n_points))
  for (cov_name in covs) {
    param_idx <- match(cov_name, param_covariates)
    w <- .msr_kernel_shape_weight(surface_binned$nodes, sigma[[param_idx]],
      if (identical(kernel, "expow")) shape[[param_idx]] else NULL, kernel)
    sums <- surface_binned$sums[[surface_binned$sources[[cov_name]]]]
    metric <- surface_binned$metrics[[cov_name]]
    power <- if (identical(metric, "sq")) 2L else 4L
    a <- vapply(seq_len(power + 1L), function(k) {
      s <- sums[[k]]
      as.numeric((if (all_rows) s else s[rows, , drop = FALSE]) %*% w)
    }, numeric(length(rows)))
    if (is.null(dim(a))) a <- matrix(a, nrow = length(rows))
    W <- a[, 1L]
    mu <- a[, 2L] / W
    m2 <- a[, 3L] / W - mu^2
    value <- switch(metric,
      sq = ifelse(m2 <= 0, 0, sqrt(pmax(m2, 0))),
      ssk = {
        m3 <- a[, 4L] / W - 3 * mu * a[, 3L] / W + 2 * mu^3
        ifelse(m2 <= 0, NA_real_, m3 / pmax(m2, .Machine$double.xmin)^1.5)
      },
      sku = {
        m4 <- a[, 5L] / W - 4 * mu * a[, 4L] / W + 6 * mu^2 * a[, 3L] / W - 3 * mu^4
        ifelse(m2 <= 0, NA_real_, m4 / pmax(m2, .Machine$double.xmin)^2 - 3)
      })
    value[!is.finite(W) | W <= 0] <- NA_real_
    out[, cov_name] <- value
  }
  out
}
