#' @title Fast Fourier Transformation
#' @description Kernel smoothing with fft
#' @param x A matrix of the raster layer
#' @param kernel A weight matrix for the smoothing
#' @param fun Either "mean" or "sum"
#' @param na.rm Logical; should NA values be removed from mean calculation?
#'
#' @return A matrix
#' @rdname fft_convolution
#' @keywords internal
#' @importFrom stats fft

fft_convolution <- function(x, kernel, fun = "mean", na.rm = TRUE) {
  x <- x_mat <- as.matrix(x)

  if (!is.matrix(x) && !is.data.frame(x)) {
    stop("'x' must be a matrix or data frame")
  }
  if (!is.matrix(kernel)) {
    stop("'kernel' must be a matrix")
  }
  if (!is.numeric(x) || !is.numeric(kernel)) {
    stop("'x' and 'kernel' must be numeric")
  }
  if (!(fun %in% c("mean", "sum"))) {
    stop("'fun' must be either 'mean' or 'sum'")
  }

  nr_x <- nrow(x_mat)
  nc_x <- ncol(x_mat)
  nr_k <- nrow(kernel)
  nc_k <- ncol(kernel)

  row_off <- floor(nr_k / 2)
  col_off <- floor(nc_k / 2)

  i1 <- row_off + 1
  i2 <- row_off + nr_x
  j1 <- col_off + 1
  j2 <- col_off + nc_x

  nr <- nr_x + nr_k - 1
  nc <- nc_x + nc_k - 1

  x_na_mask <- is.na(x_mat)
  has_na <- any(x_na_mask)

  # numerator data
  if (has_na) {
    x_mat[x_na_mask] <- 0
  }

  # pad x and kernel
  pad_x <- matrix(0, nr, nc)
  pad_k <- matrix(0, nr, nc)

  pad_x[1:nr_x, 1:nc_x] <- x_mat
  pad_k[1:nr_k, 1:nc_k] <- kernel

  # FFTs
  fft_x <- fft(pad_x)
  fft_k <- fft(pad_k)

  # numerator convolution
  conv_full <- Re(fft(fft_x * fft_k, inverse = TRUE)) / (nr * nc)
  conv_crop <- conv_full[i1:i2, j1:j2, drop = FALSE]

  # Raster edges have missing support even when the mapped cells contain no
  # NA. Always normalize means over the actual finite support when requested.
  if (fun == "mean" || !na.rm) {
    pad_valid <- matrix(0, nr, nc)
    pad_valid[1:nr_x, 1:nc_x] <- !x_na_mask
    fft_valid <- fft(pad_valid)
  }
  if (fun == "mean") {
    k_sum <- sum(kernel)
    if (k_sum == 0) stop("Sum of kernel is zero; cannot normalize by mean.")
    if (na.rm) {
      weight_full <- Re(fft(fft_valid * fft_k, inverse = TRUE)) / (nr * nc)
      weight_crop <- weight_full[i1:i2, j1:j2, drop = FALSE]
      conv_crop <- conv_crop / weight_crop
      conv_crop[weight_crop <= 0] <- NA_real_
    } else conv_crop <- conv_crop / k_sum
  }
  if (!na.rm) {
    support <- matrix(0, nr, nc)
    support[1:nr_k, 1:nc_k] <- kernel != 0
    count_full <- Re(fft(fft_valid * fft(support), inverse = TRUE)) / (nr * nc)
    # This convolution counts cells, so its exact values are integers.
    complete <- round(count_full[i1:i2, j1:j2, drop = FALSE]) == sum(kernel != 0)
    conv_crop[!complete] <- NA_real_
  }

  conv_crop[x_na_mask] <- NA_real_

  t(conv_crop)
}

# Direct sum counterpart with the same transposed output and zero padding as
# fft_convolution(). Reverse the stencil because focal uses correlation, while
# FFT uses convolution. This matters for directional adjacency/edge kernels.
.msr_focal_convolution <- function(x, kernel, fun = "sum", na.rm = TRUE) {
  if (!identical(fun, "sum"))
    stop("Direct metric convolution supports sums only.", call. = FALSE)
  x <- as.matrix(x)
  if (!any(kernel != 0)) {
    result <- matrix(0, nrow(x), ncol(x))
  } else {
    weights <- kernel[nrow(kernel):1L, ncol(kernel):1L, drop = FALSE]
    weights[weights == 0] <- NA_real_
    result <- terra::as.matrix(terra::focal(terra::rast(x), w = weights,
      fun = "sum", na.rm = na.rm, fillvalue = if (na.rm) 0 else NA_real_,
      expand = FALSE), wide = TRUE)
  }
  result[is.na(x)] <- NA_real_
  t(result)
}
