# Internal derivative adapter. Public backend support is gated separately by
# formula, row, prediction and native-likelihood checks.
.msr_prepare_surface_engine <- function(inputs, covariates, rows,
                                        standardize = TRUE, implementation = "reference") {
  specs <- inputs$scale_vars
  used <- specs[match(covariates, specs$covariate), , drop = FALSE]
  surface <- used$type == "surface" & used$metric == "sq" & used$weighted & used$optimize
  kernel <- used$type == "kernel" & used$optimize
  if (anyNA(surface | kernel) || any(!(surface | kernel)) ||
      !inputs$kernel %in% c("gaussian", "exp"))
    stop("Surface derivatives require optimized weighted SQ and smooth kernel means.", call. = FALSE)
  if (!any(surface) || anyDuplicated(covariates) || !length(rows) ||
      anyNA(rows) || any(rows != as.integer(rows) | rows < 1) ||
      is.null(inputs$d_list) || is.null(inputs$raw_cov) ||
      length(inputs$d_list) != length(inputs$raw_cov) || any(rows > length(inputs$d_list)))
    stop("Surface derivatives require matching raw buffers and retained rows.", call. = FALSE)
  means <- if (any(kernel)) .msr_prepare_kernel_engine(inputs$binned,
    covariates[kernel], rows, inputs$kernel, FALSE, implementation) else NULL
  buffers <- stats::setNames(vector("list", sum(surface)), covariates[surface])
  for (j in which(surface)) {
    buffers[[covariates[j]]] <- lapply(rows, function(i) {
      d <- inputs$d_list[[i]]
      v <- inputs$raw_cov[[i]][, used$source[j]]
      if (!is.numeric(v) || length(v) != length(d)) stop("Invalid surface buffer.", call. = FALSE)
      keep <- is.finite(d) & is.finite(v)
      if (!any(keep)) stop("Every surface row requires finite-cell support.", call. = FALSE)
      # Fully finite buffers can share the input vectors under R's copy-on-write
      # semantics. Subset only when the source-specific mask removes cells.
      if (all(keep)) list(d = d, v = as.double(v)) else
        list(d = d[keep], v = as.double(v[keep]))
    })
  }
  structure(list(covariates = covariates, p = length(covariates), n = length(rows),
    rows = rows, kernel = inputs$kernel, standardize = standardize,
    implementation = implementation, means = means, buffers = buffers,
    surface_specs = used[surface, , drop = FALSE],
    parameter_names = paste0("log_scale:", covariates)),
    class = c("msr_surface_engine", "msr_kernel_engine"))
}

.msr_eval_surface_engine <- function(engine, par) {
  if (!is.numeric(par) || length(par) != engine$p || any(!is.finite(par)) ||
      any(!is.finite(exp(par)) | exp(par) <= 0) ||
      (!is.null(names(par)) && !identical(names(par), engine$parameter_names)))
    stop("Supply finite log parameters in the surface engine order.", call. = FALSE)
  raw <- derivative <- matrix(0, engine$n, engine$p,
                              dimnames = list(NULL, engine$covariates))
  if (!is.null(engine$means)) {
    j <- match(engine$means$covariates, engine$covariates)
    z <- .msr_eval_kernel_engine(engine$means, unname(par[j]))
    raw[, j] <- z$raw; derivative[, j] <- z$d_scale
  }
  for (name in names(engine$buffers)) {
    j <- match(name, engine$covariates)
    sigma <- exp(par[j])
    if (identical(engine$implementation, "stream")) {
      z <- msr_surface_sq_derivative_cpp(engine$buffers[[name]], sigma, engine$kernel)
      raw[, j] <- z[, 1L]; derivative[, j] <- z[, 2L]
      next
    }
    for (i in seq_len(engine$n)) {
      buffer <- engine$buffers[[name]][[i]]
      ratio <- buffer$d / sigma
      w <- if (engine$kernel == "gaussian") exp(-ratio^2 / 2) else exp(-ratio)
      keep <- is.finite(w) & w > 0
      v <- buffer$v[keep]; w <- w[keep]; ratio <- ratio[keep]
      total <- sum(w)
      if (!length(v) || !is.finite(total) || total <= 0)
        stop("No effective surface support at this scale.", call. = FALSE)
      dev <- v - sum(w * v) / total
      variance <- sum(w * dev^2) / total
      raw[i, j] <- sqrt(variance)
      score <- if (engine$kernel == "gaussian") ratio^2 else ratio
      derivative[i, j] <- if (variance == 0) 0 else
        sum(w * score * (dev^2 - variance)) / (2 * total * sqrt(variance))
    }
  }
  if (any(!is.finite(raw)) || any(!is.finite(derivative)))
    stop("Nonfinite surface values or derivatives.", call. = FALSE)
  center <- if (engine$standardize) colMeans(raw) else rep(0, engine$p)
  spread <- if (engine$standardize) apply(raw, 2L, stats::sd) else rep(1, engine$p)
  if (any(!is.finite(spread) | spread <= 0))
    stop("Spatial predictors must have positive finite spread.", call. = FALSE)
  value <- sweep(sweep(raw, 2L, center), 2L, spread, "/")
  if (engine$standardize) {
    dc <- colMeans(derivative)
    centered_derivative <- sweep(derivative, 2L, dc)
    ds <- colSums(sweep(raw, 2L, center) * centered_derivative) /
      ((engine$n - 1) * spread)
    derivative <- sweep(centered_derivative, 2L, spread, "/") -
      sweep(value, 2L, ds / spread, "*")
  }
  list(raw = raw, value = value, d_scale = derivative,
       d_shape = matrix(0, engine$n, engine$p, dimnames = dimnames(raw)),
       center = center, spread = spread)
}
