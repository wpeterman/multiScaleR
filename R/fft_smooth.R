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

  if (fun == "mean") {

    # fast path: no missing data, or na.rm = FALSE
    if (!na.rm || !has_na) {
      k_sum <- sum(kernel)
      if (k_sum == 0) stop("Sum of kernel is zero; cannot normalize by mean.")
      conv_crop <- conv_crop / k_sum

    } else {
      # NA-aware denominator only when needed
      valid_mask <- matrix(1, nr_x, nc_x)
      valid_mask[x_na_mask] <- 0

      pad_valid <- matrix(0, nr, nc)
      pad_valid[1:nr_x, 1:nc_x] <- valid_mask

      fft_valid <- fft(pad_valid)
      weight_full <- Re(fft(fft_valid * fft_k, inverse = TRUE)) / (nr * nc)
      weight_crop <- weight_full[i1:i2, j1:j2, drop = FALSE]

      conv_crop <- conv_crop / weight_crop
      conv_crop[weight_crop == 0] <- NA_real_
    }
  }

  conv_crop[x_na_mask] <- NA_real_

  t(conv_crop)
}
