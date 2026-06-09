#' @title Scale Function
#' @description Scaling function to be applied to rasters
#' @param d Vector of distances
#' @param kernel Kernel function to be used ('gaussian', 'exp', 'fixed', 'expow'; Default: 'gaussian')
#' @param sigma Scaling parameter
#' @param shape Shape parameter if using exponential power kernel
#' @param r_stack.df Dataframe values extracted from rasters
#' @param output If NULL, a vector of weights is returned, otherwise a weighted raster values are returned, Default: NULL
#' @return A vector of weights or vector of weighted raster values
#' @details DETAILS
#' @examples
#' ### TO BE COMPLETED ###
#' @rdname scale_type
#' @keywords internal

scale_type <- function(d,
                       kernel = c('gaussian', 'exp', 'expow', 'fixed'),
                       sigma,
                       shape = NULL,
                       r_stack.df = NULL,
                       output = NULL) {
  kernel <- match.arg(kernel)

  # Dense raster values are weighted in vectorized R, which avoids the
  # element-by-element sparse random access in `scale_type_sparse()` (an
  # Armadillo `sp_mat(i, j)` binary search per cell) and is materially faster
  # for the one-time `kernel_prep()` pass and any exact per-cell evaluation.
  # Sparse inputs (e.g. legacy `multiScaleR_data` objects) and the weight-output
  # mode keep the original compiled path.
  if (is.null(output) && !is.null(r_stack.df) &&
      !methods::is(r_stack.df, "sparseMatrix")) {
    return(.msr_scale_type_dense(d = d,
                                 kernel = kernel,
                                 sigma = sigma,
                                 shape = shape,
                                 r_stack.df = r_stack.df))
  }

  scale_type_sparse(d = d,
                    kernel = kernel,
                    sigma_ = sigma,
                    shape_ = shape,
                    r_stack_df = r_stack.df,
                    output = output)
}


# Vectorized dense weighted mean reproducing `scale_type_sparse()` exactly:
# result[j] = sum_finite v_i w(d_i) / sum_finite w(d_i). The kernel
# normalization constant cancels in the ratio, so only the distance-dependent
# shape (`.msr_kernel_shape_weight()`) is needed. Missing cells are excluded
# from both the numerator and denominator, matching the C++ NA renormalization.
.msr_scale_type_dense <- function(d, kernel, sigma, shape = NULL, r_stack.df) {
  m <- as.matrix(r_stack.df)
  storage.mode(m) <- "double"
  nc <- ncol(m)

  if (length(sigma) == 1L && nc > 1L) {
    sigma <- rep(sigma, nc)
  }
  if (identical(kernel, "expow")) {
    if (is.null(shape)) {
      stop("shape must be provided when kernel is 'expow'.", call. = FALSE)
    }
    if (length(shape) == 1L && nc > 1L) {
      shape <- rep(shape, nc)
    }
  }

  out <- numeric(nc)
  for (j in seq_len(nc)) {
    w <- .msr_kernel_shape_weight(d, sigma[[j]],
                                  if (identical(kernel, "expow")) shape[[j]] else NULL,
                                  kernel)
    v <- m[, j]
    fin <- is.finite(v)
    den <- sum(w[fin])
    out[[j]] <- if (is.finite(den) && den > 0) sum(w[fin] * v[fin]) / den else NA_real_
  }
  out
}
