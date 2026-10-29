#' Prepare a Proper Gaussian Projection Block
#'
#' Internal building blocks for projected Gaussian scale and response effects.
#' These do not add public random-effect formula syntax.
#'
#' @param design Finite numeric observation-by-coefficient projection matrix.
#' @param covariance Proper coefficient covariance matrix. The default is an
#'   identity matrix; no estimated covariance parameters are introduced.
#' @param label Nonempty character label for the block.
#' @return A compact list containing the projection, lower Cholesky factor,
#'   fixed precision and log determinant of the Cholesky factor. Coefficients
#'   are represented as \code{exp(log_sd) * chol * z}, with independent
#'   standard-normal \code{z}. The centered density includes its normalizer.
#' @details Fixed covariance must be symmetric and positive definite. Singular
#'   or ill-conditioned matrices are rejected without covariance repair. This
#'   interface represents a proper finite Gaussian block, not an intrinsic field.
#' @keywords internal
.msr_gaussian_block <- function(design, covariance = diag(ncol(design)), label = "Gaussian block") {
  if (!is.matrix(design) || !is.numeric(design) || !nrow(design) || !ncol(design) ||
      any(!is.finite(design))) stop("Gaussian block design must be a finite nonempty numeric matrix.", call. = FALSE)
  q <- ncol(design)
  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      !identical(dim(covariance), c(q, q)) || any(!is.finite(covariance)) ||
      !isSymmetric(covariance, tol = 1e-12))
    stop("Gaussian block covariance must be finite, symmetric and match the coefficient dimension.", call. = FALSE)
  L <- tryCatch(t(chol(covariance)), error = function(e) NULL)
  if (is.null(L) || rcond(covariance) < 1e-12)
    stop("Gaussian block covariance must be positive definite and well conditioned; no covariance was repaired.", call. = FALSE)
  if (!is.character(label) || length(label) != 1L || is.na(label) || !nzchar(label))
    stop("Supply a nonempty Gaussian block label.", call. = FALSE)
  list(design = design, chol = L, precision = chol2inv(t(L)),
    log_chol_determinant = sum(log(diag(L))), dimension = q, label = label)
}

# Whitened latent coordinates have a standard-normal density already included
# by the shared posterior/likelihood objectives. Do not add a second normalizer.
.msr_gaussian_effect <- function(block, z, log_sd) {
  if (length(z) != block$dimension || length(log_sd) != 1L)
    stop("Gaussian block parameters do not match the prepared dimension.", call. = FALSE)
  exp(log_sd) * drop(block$design %*% (block$chol %*% z))
}

# Independent centered-coordinate reference, also valid for automatic
# differentiation. The covariance is fixed; only b and log_sd are differentiated.
.msr_gaussian_centered_nll <- function(block, b, log_sd) {
  if (length(b) != block$dimension || length(log_sd) != 1L)
    stop("Centered Gaussian parameters do not match the prepared dimension.", call. = FALSE)
  .5 * exp(-2 * log_sd) * sum(b * drop(block$precision %*% b)) +
    block$dimension * (log_sd + .5 * log(2 * pi)) + block$log_chol_determinant
}

.msr_has_gaussian_blocks <- function(spec) {
  terms <- c(unlist(lapply(spec$scales, function(sm) sm$random), recursive = FALSE), spec$response$random)
  any(vapply(terms, function(term) !is.null(term$gaussian_block), TRUE))
}
