validate_scalar_logical <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", arg), call. = FALSE)
  }
}

validate_scalar_numeric <- function(x,
                                    arg,
                                    lower = -Inf,
                                    upper = Inf,
                                    inclusive_lower = TRUE,
                                    inclusive_upper = TRUE,
                                    integerish = FALSE,
                                    positive = FALSE) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x)) {
    stop(sprintf("`%s` must be a single finite numeric value.", arg), call. = FALSE)
  }

  if (integerish && x != as.integer(x)) {
    stop(sprintf("`%s` must be a whole number.", arg), call. = FALSE)
  }

  if (positive && x <= 0) {
    stop(sprintf("`%s` must be > 0.", arg), call. = FALSE)
  }

  lower_ok <- if (inclusive_lower) x >= lower else x > lower
  upper_ok <- if (inclusive_upper) x <= upper else x < upper

  if (!lower_ok || !upper_ok) {
    if (is.finite(lower) && is.finite(upper)) {
      comparator <- paste0(
        if (inclusive_lower) "[" else "(",
        lower, ", ", upper,
        if (inclusive_upper) "]" else ")"
      )
      stop(sprintf("`%s` must be within %s.", arg, comparator), call. = FALSE)
    } else if (is.finite(lower)) {
      comparator <- if (inclusive_lower) paste0(">= ", lower) else paste0("> ", lower)
      stop(sprintf("`%s` must be %s.", arg, comparator), call. = FALSE)
    } else {
      comparator <- if (inclusive_upper) paste0("<= ", upper) else paste0("< ", upper)
      stop(sprintf("`%s` must be %s.", arg, comparator), call. = FALSE)
    }
  }
}

validate_numeric_vector <- function(x,
                                    arg,
                                    length_ = NULL,
                                    positive = FALSE) {
  if (!is.numeric(x) || length(x) == 0 || anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must be a numeric vector with finite values.", arg), call. = FALSE)
  }

  if (!is.null(length_) && length(x) != length_) {
    stop(sprintf("`%s` must have length %s.", arg, length_), call. = FALSE)
  }

  if (positive && any(x <= 0)) {
    stop(sprintf("`%s` must contain only values > 0.", arg), call. = FALSE)
  }
}

validate_character_scalar <- function(x, arg) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    stop(sprintf("`%s` must be a non-empty character string.", arg), call. = FALSE)
  }
}

validate_multiScaleR_input <- function(x, arg = "multiScaleR") {
  if (!inherits(x, c("multiScaleR", "multiScaleR_data"))) {
    stop(sprintf("`%s` must be a `multiScaleR` or `multiScaleR_data` object.", arg), call. = FALSE)
  }
}

validate_covariates_before_scale <- function(x, context = "calculated covariates") {
  x <- as.matrix(x)
  if (ncol(x) == 0) {
    return(invisible(TRUE))
  }

  for (j in seq_len(ncol(x))) {
    covariate <- colnames(x)[[j]]
    if (is.null(covariate) || is.na(covariate) || !nzchar(covariate)) {
      covariate <- paste0("column ", j)
    }

    values <- x[, j]
    non_missing <- !is.na(values)
    finite <- is.finite(values)

    if (!any(non_missing)) {
      stop(
        sprintf("Cannot scale %s: covariate `%s` is all NA.", context, covariate),
        call. = FALSE
      )
    }
    if (any(non_missing & !finite)) {
      stop(
        sprintf("Cannot scale %s: covariate `%s` contains non-finite values.", context, covariate),
        call. = FALSE
      )
    }
    if (!any(finite)) {
      stop(
        sprintf("Cannot scale %s: covariate `%s` has no finite values.", context, covariate),
        call. = FALSE
      )
    }

    finite_values <- values[finite]
    if (length(unique(finite_values)) < 2 || isTRUE(stats::sd(finite_values) == 0)) {
      stop(
        sprintf("Cannot scale %s: covariate `%s` has zero variance.", context, covariate),
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}
