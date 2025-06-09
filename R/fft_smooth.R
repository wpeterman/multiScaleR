#' @title Fast Fourier Transformation
#' @description Kernel smoothing with fft
#' @param x A matrix of the raster layer
#' @param kernel A weight matrix for the smoothing
#' @param fun Only calculates mean
#' @param na.rm Always remove NA values from mean calculation

#' @return A matrix
#' @rdname fft_convolution
#' @keywords internal
#' @importFrom stats fft

fft_convolution <- function(x, kernel, fun = "mean", na.rm = TRUE) {
  # Check inputs
    # Convert input to matrix
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

  # Store NA positions and replace NAs with 0
  x_na_mask <- is.na(x_mat)
  x_mat[x_na_mask] <- 0

  # Determine padded size for convolution
  nr <- nrow(x_mat) + nrow(kernel) - 1
  nc <- ncol(x_mat) + ncol(kernel) - 1

  # Zero-pad x and kernel
  pad_x <- matrix(0, nr, nc)
  pad_k <- matrix(0, nr, nc)
  pad_x[1:nrow(x_mat), 1:ncol(x_mat)] <- x_mat
  pad_k[1:nrow(kernel), 1:ncol(kernel)] <- kernel

  # Perform FFT-based convolution
  fft_x <- fft(pad_x)
  fft_k <- fft(pad_k)
  conv_full <- Re(fft(fft_x * fft_k, inverse = TRUE)) / (nr * nc)

  # Determine crop boundaries to match original size
  row_off <- floor(nrow(kernel) / 2)
  col_off <- floor(ncol(kernel) / 2)
  i1 <- row_off + 1
  i2 <- row_off + nrow(x_mat)
  j1 <- col_off + 1
  j2 <- col_off + ncol(x_mat)

  # Crop to original dimensions
  conv_crop <- conv_full[i1:i2, j1:j2]

  # Normalize if 'mean' is requested
  if (fun == "mean") {
    k_sum <- sum(kernel, na.rm = na.rm)
    if (k_sum == 0) stop("Sum of kernel is zero; cannot normalize by mean.")
    conv_crop <- conv_crop / k_sum
  }

  # Restore NA values
  conv_crop[x_na_mask] <- NA

  # Return result with transposition to match original orientation
  return(t(conv_crop))
}
