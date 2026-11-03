# The hurdle response is one overall expected count from two predictors.
# Keep its native mean; propagate signed covariance across all coefficients.
.msr_gam_ziplss_model <- function(object) {
  if (!identical(object$family$family, "ziplss") || !is.list(object$pterms) ||
      inherits(object$pterms, "terms") || length(object$pterms) != 2L ||
      !requireNamespace("mgcv", quietly = TRUE)) return(FALSE)
  # A family name does not qualify a distinct provider's uncertainty method.
  native <- utils::getS3method("predict", "gam", envir = asNamespace("mgcv"))
  for (kind in class(object)) {
    method <- utils::getS3method("predict", kind, optional = TRUE)
    if (!is.null(method)) return(identical(method, native))
  }
  FALSE
}

.msr_gam_ziplss_response_metadata <- function(object, args) {
  if (!.msr_gam_ziplss_model(object) || !identical(args$type, "response") ||
      !isTRUE(args$se.fit)) return(object)
  original <- object$family$predict
  reference <- mgcv::ziplss()$predict
  if (!is.function(original) || !identical(formals(original), formals(reference)) ||
      !identical(body(original), body(reference)))
    stop("Cannot verify response uncertainty for this custom ziplss predictor.", call. = FALSE)
  object$family$predict <- .msr_gam_ziplss_predictor(original)
  object
}

.msr_gam_ziplss_predictor <- function(native_prediction) {
  function(family, se = FALSE, eta = NULL, y = NULL, X = NULL,
           beta = NULL, off = NULL, Vb = NULL) {
    point <- native_prediction(family, se = FALSE, eta = eta, y = y,
      X = X, beta = beta, off = off, Vb = Vb)
    if (!isTRUE(se)) return(point)
    if (!is.null(eta) || !is.matrix(X) || !is.numeric(X) ||
        any(!is.finite(X)) || !is.numeric(beta) ||
        length(beta) != ncol(X) || any(!is.finite(beta)) ||
        !is.matrix(Vb) || !identical(dim(Vb), c(ncol(X), ncol(X))) ||
        any(!is.finite(Vb)))
      stop("Cannot verify the native ziplss coefficient covariance and design.")
    indices <- attr(X, "lpi"); n <- nrow(X)
    if (!is.list(indices) || length(indices) != 2L ||
        any(!vapply(indices, function(i) is.numeric(i) && length(i) > 0L &&
          all(is.finite(i)) && all(i == as.integer(i)) && !anyDuplicated(i) &&
          all(i >= 1L & i <= ncol(X)), logical(1))))
      stop("Cannot align the native ziplss coefficient spans.")
    if (is.null(off)) off <- list(0, 0)
    off[[3L]] <- 0
    for (j in 1:2) {
      if (is.null(off[[j]])) off[[j]] <- 0
      if (!is.numeric(off[[j]]) || !length(off[[j]]) ||
          !length(off[[j]]) %in% c(1L, n) || anyNA(off[[j]]))
        stop("Cannot align the native ziplss prediction offsets.")
    }
    gamma <- drop(X[, indices[[1]], drop = FALSE] %*% beta[indices[[1]]]) + off[[1]]
    occurrence <- drop(X[, indices[[2]], drop = FALSE] %*% beta[indices[[2]]]) + off[[2]]
    lambda <- exp(gamma); intensity <- exp(occurrence)
    probability <- -expm1(-intensity)
    mu <- rep(1, n); dmu <- numeric(n)
    active <- gamma >= log(.Machine$double.eps) / 2
    mu[active] <- lambda[active] / -expm1(-lambda[active])
    small <- active & lambda < .001
    a <- lambda[small]
    dmu[small] <- a / 2 + a^2 / 6 - a^4 / 180 + a^6 / 5040
    large <- active & !small
    dmu[large] <- mu[large] * (1 - lambda[large] / expm1(lambda[large]))
    dprobability <- numeric(n)
    finite <- is.finite(occurrence)
    dprobability[finite] <- exp(occurrence[finite] - intensity[finite])
    fit <- point[[1L]]
    if (!is.numeric(fit) || length(fit) != n || any(!is.finite(fit)) ||
        any(!is.finite(probability * mu)) ||
        any(abs(as.numeric(fit) - probability * mu) > 1e-8 * pmax(1, abs(fit))))
      stop("Cannot verify the native ziplss response mean at these inputs.")
    # Add both derivatives for any coefficient shared by the two predictors.
    gradient <- matrix(0, n, ncol(X), dimnames = dimnames(X))
    derivatives <- list(probability * dmu, mu * dprobability)
    for (j in 1:2) gradient[, indices[[j]]] <-
      gradient[, indices[[j]], drop = FALSE] +
      sweep(X[, indices[[j]], drop = FALSE], 1L, derivatives[[j]], `*`)
    variance <- rowSums((gradient %*% Vb) * gradient)
    if (any(!is.finite(variance)) || any(variance < 0))
      stop("Native ziplss response uncertainty requires finite nonnegative variance.")
    standard_error <- sqrt(variance)
    attributes(standard_error) <- attributes(fit)
    list(fit = fit, se.fit = standard_error)
  }
}
