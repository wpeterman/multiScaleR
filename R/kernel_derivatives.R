# Internal, experimental smooth-kernel engine. Preparation fixes the fitting
# rows and parameter order once. Counts refer to finite cells, including zeros.
# These helpers do not change the public optimization backend.
.msr_prepare_kernel_engine <- function(binned, covariates = names(binned$vsum),
                                       rows = seq_len(nrow(binned$vsum[[1L]])),
                                       kernel = "gaussian", standardize = TRUE,
                                       implementation = "reference") {
  implementation <- match.arg(implementation, c("reference", "stream"))
  if (!kernel %in% c("gaussian", "exp", "expow") || length(kernel) != 1L)
    stop("The derivative engine requires a smooth kernel.", call. = FALSE)
  if (!is.logical(standardize) || length(standardize) != 1L || is.na(standardize))
    stop("standardize must be TRUE or FALSE.", call. = FALSE)
  d <- binned$bin_rep
  if (!is.numeric(d) || !length(d) || any(!is.finite(d) | d < 0))
    stop("Bin distances must be finite and nonnegative.", call. = FALSE)
  if (!length(covariates) || anyNA(covariates) || anyDuplicated(covariates) ||
      !all(covariates %in% names(binned$vsum)) ||
      !all(covariates %in% names(binned$csum)))
    stop("Covariates must name distinct sum and count matrices.", call. = FALSE)
  nr <- nrow(binned$vsum[[covariates[1L]]])
  if (!is.numeric(rows) || !length(rows) || anyNA(rows) ||
      any(rows != as.integer(rows) | rows < 1 | rows > nr))
    stop("rows must contain valid positive row indices.", call. = FALSE)
  if (standardize && length(rows) < 2L)
    stop("Standardization requires at least two rows.", call. = FALSE)
  V <- C <- stats::setNames(vector("list", length(covariates)), covariates)
  for (nm in covariates) {
    v <- binned$vsum[[nm]]
    c <- binned$csum[[nm]]
    if (!is.matrix(v) || !is.numeric(v) || !is.matrix(c) || !is.numeric(c) ||
        !identical(dim(v), dim(c)) || nrow(v) != nr || ncol(v) != length(d))
      stop("Invalid finite-cell sum/count matrices for ", nm, ".", call. = FALSE)
    if (implementation == "stream") {
      checked <- msr_kernel_validate_cpp(v, c, as.integer(rows))
      if (!checked[1L])
        stop("Invalid finite-cell sum/count matrices for ", nm, ".", call. = FALSE)
      if (!checked[2L])
        stop("Every selected row must have finite-cell support.", call. = FALSE)
    } else if (any(!is.finite(v)) || any(!is.finite(c) | c < 0) ||
               any(v[c == 0] != 0))
      stop("Invalid finite-cell sum/count matrices for ", nm, ".", call. = FALSE)
    if (implementation == "stream") {
      # Share unchanged double matrices. R copy-on-modify protects the caller;
      # selected/reordered rows still get their own correctly ordered matrix.
      all_rows <- identical(as.integer(rows), seq_len(nr))
      V[[nm]] <- if (all_rows) v else v[rows, , drop = FALSE]
      C[[nm]] <- if (all_rows) c else c[rows, , drop = FALSE]
      if (!is.double(V[[nm]])) storage.mode(V[[nm]]) <- "double"
      if (!is.double(C[[nm]])) storage.mode(C[[nm]]) <- "double"
    } else {
      V[[nm]] <- v[rows, , drop = FALSE] * 1.0
      C[[nm]] <- c[rows, , drop = FALSE] * 1.0
    }
    if (implementation == "reference" && any(rowSums(C[[nm]]) <= 0))
      stop("Every selected row must have finite-cell support.", call. = FALSE)
  }
  structure(list(V = V, C = C, distance = as.double(d),
                 covariates = covariates, rows = rows, n = length(rows),
                 p = length(covariates), kernel = kernel,
                 standardize = standardize, implementation = implementation,
                 parameter_names = c(paste0("log_scale:", covariates),
                   if (kernel == "expow") paste0("log_shape:", covariates))),
            class = "msr_kernel_engine")
}

.msr_eval_kernel_engine <- function(engine, par) {
  if (inherits(engine, "msr_surface_engine")) return(.msr_eval_surface_engine(engine, par))
  if (!inherits(engine, "msr_kernel_engine"))
    stop("Use a prepared kernel engine.", call. = FALSE)
  if (!is.numeric(par) || length(par) != length(engine$parameter_names) ||
      any(!is.finite(par)) || any(!is.finite(exp(par)) | exp(par) == 0))
    stop("Supply finite log parameters in the engine's parameter order.", call. = FALSE)
  if (!is.null(names(par)) && !identical(names(par), engine$parameter_names))
    stop("Parameter names must match the prepared engine order.", call. = FALSE)
  p <- engine$p
  evaluator <- if (identical(engine$implementation, "stream"))
    msr_kernel_stream_cpp else msr_kernel_eval_cpp
  ans <- evaluator(
    engine$V, engine$C, engine$distance, par[seq_len(p)],
    if (engine$kernel == "expow") par[p + seq_len(p)] else rep(0, p),
    match(engine$kernel, c("gaussian", "exp", "expow")) - 1L,
    engine$standardize)
  for (nm in c("raw", "value", "d_scale", "d_shape"))
    colnames(ans[[nm]]) <- engine$covariates
  ans
}
