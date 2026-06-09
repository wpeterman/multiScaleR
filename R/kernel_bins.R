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
    if (!is.null(eval_idx)) {
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
